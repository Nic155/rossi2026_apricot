# =============================================================================
# Figure 2C: circos plot of selective sweeps and enriched deleterious-load terms
#
# Draws a circular genome plot (chr1-8) showing the selective-sweep regions
# together with the genomic locations of the W2-enriched GO/KEGG terms whose
# deleterious-variant burden is concentrated in the chromosome-4 sweep hotspot.
#
# Input:
#   circos.ontology.RData, providing
#     - go_variants, kegg_variants : positions (CHR, POS) of deleterious
#       variants annotated to each significantly enriched term (from 04.R)
#     - regions_tbl                : selective-sweep regions (chr, start, end)
#
# Output (to figures_final/):
#   circos_enriched_terms_vs_sweeps.png - Fig 2C
#   circos_track_map.csv                - track-number -> term/colour legend
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(circlize)
library(scales)

# NOTE [repo path]: replace with a repo-relative path, e.g. "data/05.RData"
load("~/Paper parmeniaca/scripts/final/final final/rossi2026_apricot/05.RData")
# =========================================================
# SETTINGS
# =========================================================
# NOTE [repo path]: replace with a repo-relative output dir, e.g. "figures"
out_dir <- "C:/Users/nirossi/Documents/Paper parmeniaca/figures_final"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(out_dir, "circos_enriched_terms_vs_sweeps.png")
out_map <- file.path(out_dir, "circos_track_map.csv")

min_variants <- 10L        # keep enriched terms with >= 10 variants genome-wide
window_size  <- 10000L     # 10 kb bins for term-occupancy blocks
chr_levels   <- paste0("chr", 1:8)

# SAME style for every data track
track_height_all <- 0.055
track_bg_border  <- "grey85"
track_alpha      <- 0.35
track_lwd        <- 0.5

region_col <- "#2CA25F"

# chromosome lengths (apricot Marouch v3.1)
chr_sizes <- tibble(
  chr = paste0("chr", 1:8),
  chr_len = c(44417667, 26986824, 22808846, 25153484,
              16905379, 26898547, 20744237, 20012128)
)

# =========================================================
# HELPERS
# =========================================================
make_windows <- function(chr_sizes, window_size) {
  bind_rows(lapply(seq_len(nrow(chr_sizes)), function(i) {
    starts <- seq.int(1L, chr_sizes$chr_len[i], by = window_size)
    ends   <- pmin(starts + window_size - 1L, chr_sizes$chr_len[i])
    tibble(chr = chr_sizes$chr[i], bin = seq_along(starts),
           start = starts, end = ends)
  }))
}

draw_rect_track <- function(track_df, track_colour) {
  track_df <- as.data.frame(track_df)
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = track_height_all,
    bg.border = track_bg_border,
    panel.fun = function(x, y) {
      sector_chr <- CELL_META$sector.index
      d <- track_df[track_df$chr == sector_chr, , drop = FALSE]
      if (nrow(d) > 0) {
        d <- d[order(d$start), , drop = FALSE]
        circos.rect(
          xleft = d$start, ybottom = 0,
          xright = d$end,  ytop = 1,
          col = alpha(track_colour, track_alpha),
          border = track_colour, lwd = track_lwd
        )
      }
    }
  )
}

# =========================================================
# 1) BUILD CATEGORY POSITION TABLE from go_variants + kegg_variants
#    Each enriched term becomes one "category"
# =========================================================
classified_pos <- bind_rows(
  go_variants   %>% transmute(chr = CHR, pos = POS, pathway_category = description),
  kegg_variants %>% transmute(chr = CHR, pos = POS, pathway_category = description)
) %>%
  filter(!is.na(pathway_category), !is.na(chr), !is.na(pos)) %>%
  distinct(chr, pos, pathway_category)

category_summary <- classified_pos %>%
  count(pathway_category, name = "n_variants", sort = TRUE)

cat("\nVariants per enriched term:\n")
print(category_summary)

keep_categories <- category_summary %>%
  filter(n_variants >= min_variants) %>%
  pull(pathway_category)

if (length(keep_categories) == 0) {
  stop("No enriched terms with n_variants >= ", min_variants, ".")
}

classified_pos <- classified_pos %>%
  filter(pathway_category %in% keep_categories)

# =========================================================
# 2) RESTRICT TO chr1-8 AND PREP chr_sizes / regions
# =========================================================
chr_sizes <- chr_sizes %>%
  filter(chr %in% chr_levels) %>%
  mutate(chr = factor(chr, levels = chr_levels)) %>%
  arrange(chr) %>%
  mutate(chr = as.character(chr))

classified_pos <- classified_pos %>% filter(chr %in% chr_sizes$chr)
regions_tbl    <- regions_tbl    %>% filter(chr %in% chr_sizes$chr)

xlim_mat <- cbind(rep(1, nrow(chr_sizes)), chr_sizes$chr_len)

# =========================================================
# 3) CATEGORY BLOCKS FROM WINDOW OCCUPANCY
# =========================================================
# For each term, mark the 10 kb bins it occupies, then merge runs of
# consecutive occupied bins into contiguous blocks (block_id increments
# wherever the bin index is non-consecutive). Each block becomes one rectangle.
windows <- make_windows(chr_sizes, window_size)

category_blocks <- classified_pos %>%
  mutate(bin = ((pos - 1L) %/% window_size) + 1L) %>%
  count(pathway_category, chr, bin, name = "n_hits") %>%
  filter(n_hits > 0) %>%
  inner_join(windows, by = c("chr", "bin")) %>%
  arrange(pathway_category, chr, bin) %>%
  group_by(pathway_category, chr) %>%
  mutate(block_id = cumsum(c(TRUE, diff(bin) != 1L))) %>%
  group_by(pathway_category, chr, block_id) %>%
  summarise(start = min(start), end = max(end), .groups = "drop")

active_categories <- category_blocks %>% distinct(pathway_category) %>% pull(pathway_category)
keep_categories <- keep_categories[keep_categories %in% active_categories]
category_blocks <- category_blocks %>% filter(pathway_category %in% keep_categories)

# =========================================================
# 4) COLOURS
# =========================================================
cat_cols <- setNames(
  hcl.colors(length(keep_categories), palette = "Set 2"),
  keep_categories
)

# =========================================================
# 5) TRACK LIST (Selected sweep regions first, then each term)
# =========================================================
track_list <- c(
  list(list(
    name   = "Selective sweep regions",
    data   = regions_tbl %>% transmute(chr, start, end),
    colour = region_col
  )),
  lapply(keep_categories, function(cat_name) {
    list(
      name   = cat_name,
      data   = category_blocks %>%
        filter(pathway_category == cat_name) %>%
        transmute(chr, start, end),
      colour = unname(cat_cols[cat_name])
    )
  })
)

track_map <- tibble(
  track_no   = 0:length(track_list),
  track_name = c("Chromosome labels / axis",
                 vapply(track_list, function(x) x$name, character(1))),
  colour     = c(NA_character_, vapply(track_list, function(x) x$colour, character(1))),
  n_blocks   = c(NA_integer_, vapply(track_list, function(x) nrow(x$data), integer(1)))
)
write.csv(track_map, out_map, row.names = FALSE)

cat("\n================ TRACK MAP ================\n")
print(track_map, n = Inf)

# =========================================================
# 6) PLOT
# =========================================================
png(
  filename = out_png,
  width = 12,
  height = 12,
  units = "in",
  res = 600,
  bg = "transparent"      # was "white"
)
par(mar = c(0.5, 0.5, 0.5, 0.5))

circos.clear()
circos.par(
  start.degree = 90, gap.degree = 3,
  cell.padding = c(0, 0, 0, 0),
  track.margin = c(0.0015, 0.0015),
  points.overflow.warning = FALSE
)

circos.initialize(factors = chr_sizes$chr, xlim = xlim_mat)

# Axis track
circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.06, bg.border = NA,
  panel.fun = function(x, y) {
    xlim <- CELL_META$xlim
    sector_chr <- CELL_META$sector.index
    circos.axis(h = "top", labels.cex = 0.4,
                major.at = pretty(xlim, n = 4),
                minor.ticks = 0, major.tick.length = mm_y(1.2))
    circos.text(x = mean(xlim), y = 0.2, labels = sector_chr,
                cex = 0.6, facing = "bending.inside", niceFacing = TRUE)
  }
)

# data tracks, identical style, colour differs
for (tr in track_list) {
  message("Drawing: ", tr$name)
  draw_rect_track(tr$data, tr$colour)
}

dev.off()
circos.clear()

cat("\nSaved PNG to:\n", out_png, "\n")
cat("Saved track map to:\n", out_map, "\n")