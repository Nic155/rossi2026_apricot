# =============================================================================
# Figure 1 (D, E): predicted severity of homozygous deleterious variants
#
# Compares the predicted functional severity of homozygous derived deleterious
# variants between the cultivated (C1) and wild (W2) gene pools.
#
# IMPORTANT: the "enriched regions" defined here (top 5% of 100 kb windows by
# homozygous-genotype frequency, per group) are a visualisation device for
# Fig 1D/E ONLY. They are NOT the selective-sweep regions of Fig 2 (which are
# defined upstream from the pi_W2/pi_C1 diversity ratio). Do not conflate them.
#
# Input:
#   03.RData, providing data frame `df` with, per SNP:
#     - CHR, POS
#     - PROVEAN_SCORE          : lower = stronger predicted deleterious effect
#     - DEL_CLASS, IMPACT      : effect-class annotations (see Fig 1E below)
#     - one column per individual named "<sample>__C1" or "<sample>__W2",
#       genotype coded as dosage of the derived allele (0/1/2; NA = missing)
#
# Output (PNG, to figures_final/):
#   p_density.png  - Fig 1D: normalised PROVEAN-score density, C1 vs W2
#   Proportion of start stop codon-disrupting homozygous variants.png
#                  - Fig 1E: per-individual proportion of homozygous loss-of-
#                    function (start/stop-disrupting) variants
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# NOTE [repo path]: replace with a repo-relative path, e.g. "data/03.RData"
load("~/Paper parmeniaca/scripts/final/final final/rossi2026_apricot/03.RData")

# NOTE [repo path]: replace with a repo-relative output dir, e.g. "figures"
out_dir <- "C:/Users/nirossi/Documents/Paper parmeniaca/figures_final"

# -----------------------------
# Identify sample columns
# -----------------------------
# Per-individual genotype columns are suffixed by gene pool (__C1 / __W2)

c1_cols <- grep("__C1$", names(df), value = TRUE)
w2_cols <- grep("__W2$", names(df), value = TRUE)

# -----------------------------
# Define group colours
# -----------------------------
# C1 = cultivated (orange), W2 = wild (blue)

my_cols <- c(
  C1 = "#D55E00",
  W2 = "#0072B2"
)

# ============================================================
# p_density : Normalized PROVEAN density in
#             homozygous-enriched regions
# ============================================================

# -----------------------------
# Settings for enriched regions
# -----------------------------
# Fig 1D restricts the PROVEAN density to windows where homozygous variants are
# concentrated, so the comparison reflects realised (recessive) load. These are
# the top 5% of windows by homozygous frequency, computed separately per group.

window_size <- 100000      # 100 kb windows
enrichment_quantile <- 0.95   # keep the top 5% of windows (95th-percentile cutoff)
min_hom_variants <- 10        # require >=10 homozygous-variant rows per window

# -----------------------------
# Add genomic windows
# -----------------------------
# Assign each SNP to a non-overlapping 100 kb window by position

df_win <- df %>%
  mutate(
    window_start = floor((POS - 1) / window_size) * window_size + 1,
    window_end = window_start + window_size - 1
  )

# -----------------------------
# Identify regions enriched with homozygous mutations
# -----------------------------
# Per group and window: hom_freq = homozygous derived calls / callable genotypes.
# row_hom_calls counts dosage==2 genotypes; row_callable counts non-missing
# genotypes. The 95th-percentile cutoff is taken within each group separately.

region_enrichment <- bind_rows(
  df_win %>%
    transmute(
      group = "C1",
      CHR,
      window_start,
      window_end,
      row_hom_calls = rowSums(across(all_of(c1_cols), ~ !is.na(.x) & .x == 2)),
      row_callable = rowSums(across(all_of(c1_cols), ~ !is.na(.x)))
    ),

  df_win %>%
    transmute(
      group = "W2",
      CHR,
      window_start,
      window_end,
      row_hom_calls = rowSums(across(all_of(w2_cols), ~ !is.na(.x) & .x == 2)),
      row_callable = rowSums(across(all_of(w2_cols), ~ !is.na(.x)))
    )
) %>%
  group_by(group, CHR, window_start, window_end) %>%
  summarise(
    hom_calls = sum(row_hom_calls),
    callable = sum(row_callable),
    n_hom_variant_rows = sum(row_hom_calls > 0),
    hom_freq = hom_calls / callable,
    .groups = "drop"
  ) %>%
  filter(
    callable > 0,
    n_hom_variant_rows >= min_hom_variants
  ) %>%
  group_by(group) %>%
  mutate(
    enrichment_cutoff = quantile(hom_freq, enrichment_quantile, na.rm = TRUE),
    enriched = hom_freq >= enrichment_cutoff
  ) %>%
  ungroup()

enriched_regions <- region_enrichment %>%
  filter(enriched) %>%
  select(group, CHR, window_start, window_end, hom_freq, hom_calls, callable)

# -----------------------------
# Extract PROVEAN variants in enriched homozygous regions
# -----------------------------
# Keep a variant for a group's density if it is present in that group (dosage>0
# in >=1 individual). `present` flags any-state presence; `hom_present` flags
# homozygous presence. The inner_join restricts to the enriched windows above.

provean_in_enriched <- bind_rows(
  df_win %>%
    transmute(
      group = "C1",
      CHR,
      window_start,
      window_end,
      PROVEAN_SCORE,
      present = rowSums(across(all_of(c1_cols), ~ !is.na(.x) & .x > 0)) > 0,
      hom_present = rowSums(across(all_of(c1_cols), ~ !is.na(.x) & .x == 2)) > 0
    ),

  df_win %>%
    transmute(
      group = "W2",
      CHR,
      window_start,
      window_end,
      PROVEAN_SCORE,
      present = rowSums(across(all_of(w2_cols), ~ !is.na(.x) & .x > 0)) > 0,
      hom_present = rowSums(across(all_of(w2_cols), ~ !is.na(.x) & .x == 2)) > 0
    )
) %>%
  inner_join(
    enriched_regions %>%
      select(group, CHR, window_start, window_end),
    by = c("group", "CHR", "window_start", "window_end")
  ) %>%
  filter(
    present,
    !is.na(PROVEAN_SCORE)
  ) %>%
  mutate(
    group = factor(group, levels = c("C1", "W2"))
  )

# Split groups so W2 is drawn first and C1 overlaid on top (C1 is the focal
# comparison; its curve should sit above W2's where they overlap)
provean_w2_enriched <- provean_in_enriched %>%
  filter(group == "W2")

provean_c1_enriched <- provean_in_enriched %>%
  filter(group == "C1")

# -----------------------------
# p_density: normalized density
# -----------------------------

p_density <- ggplot() +
  geom_density(
    data = provean_w2_enriched,
    aes(
      x = PROVEAN_SCORE,
      color = group
    ),
    fill = NA,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_density(
    data = provean_c1_enriched,
    aes(
      x = PROVEAN_SCORE,
      color = group
    ),
    fill = NA,
    linewidth = 0.9,
    show.legend = FALSE
  ) +
  scale_color_manual(values = my_cols, breaks = c("C1", "W2")) +
  theme_bw() +
  labs(
    x = "PROVEAN score",
    y = "Normalized Density"
  ) +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1.2
    )
  )

p_density

ggsave(
  # NOTE [repo path]: out_dir is an absolute path; set it to "figures" for the repo
  filename = file.path(out_dir, "p_density.png"),
  plot = p_density,
  bg = "transparent",
  width = 12,
  height = 8,
  dpi = 300
)

# ============================================================
# p_hom_delclass_box (Fig 1E):
# Per individual, proportion of homozygous loss-of-function variants =
#   (homozygous start/stop-disrupting variants)
#   ---------------------------------------------------
#   (all homozygous deleterious variants: HIGH-impact + PROVEAN-flagged)
# A lower proportion in C1 indicates preferential removal of large-effect
# (start/stop-disrupting) recessive variants.
#   Numerator   = DEL_CLASS in {HIGH, STOP_START_LOSS}
#   Denominator = IMPACT    in {HIGH, PROVEAN}
# Both counted in the homozygous state only (dosage == 2).
# ============================================================

sample_cols <- c(c1_cols, w2_cols)

# -----------------------------
# Numerator: DEL_CLASS HIGH + STOP_START_LOSS
# -----------------------------

numerator_rows <- !is.na(df$DEL_CLASS) &
  df$DEL_CLASS %in% c(
    "HIGH",
    "STOP_START_LOSS"
  )

# -----------------------------
# Denominator: IMPACT HIGH + PROVEAN
# -----------------------------

denominator_rows <- !is.na(df$IMPACT) &
  df$IMPACT %in% c(
    "HIGH",
    "PROVEAN"
  )

numerator_idx <- which(numerator_rows)
denominator_idx <- which(denominator_rows)

# -----------------------------
# Homozygous-only counting function (dosage == 2)
# -----------------------------
# For each individual column, count rows (at the given indices) where the
# genotype is homozygous derived; NA genotypes are ignored.

count_hom_per_sample <- function(idx, cols) {
  vapply(cols, function(s) {
    g <- df[[s]][idx]
    sum(g == 2, na.rm = TRUE)
  }, numeric(1))
}

# -----------------------------
# Per-sample counts
# -----------------------------

hom_numerator <- count_hom_per_sample(
  idx = numerator_idx,
  cols = sample_cols
)

hom_denominator <- count_hom_per_sample(
  idx = denominator_idx,
  cols = sample_cols
)

# -----------------------------
# Build plotting dataframe
# -----------------------------

hom_delclass_prop <- data.frame(
  sample = sample_cols,
  numerator = as.numeric(hom_numerator),
  denominator = as.numeric(hom_denominator)
) %>%
  mutate(
    group = case_when(
      grepl("__C1$", sample) ~ "C1",
      grepl("__W2$", sample) ~ "W2",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = c("C1", "W2")),
    proportion = if_else(
      denominator > 0,
      numerator / denominator,
      NA_real_
    )
  ) %>%
  filter(!is.na(group))

# -----------------------------
# Y-axis limits
# -----------------------------

y_max <- max(hom_delclass_prop$proportion, na.rm = TRUE)

if (!is.finite(y_max) || y_max == 0) {
  y_max <- 0.01
}

# -----------------------------
# Boxplot theme
# -----------------------------

box_theme <- theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

# -----------------------------
# Boxplot
# -----------------------------

p_hom_delclass_box <- ggplot(
  hom_delclass_prop,
  aes(x = group, y = proportion, fill = group)
) +
  geom_boxplot(
    width         = 0.65,
    alpha         = 0.85,
    outlier.shape = NA,
    color         = "black",
    coef          = Inf,
    na.rm         = TRUE
  ) +
  geom_jitter(
    color       = "grey20",
    width       = 0.12,
    height      = 0,
    size        = 1.8,
    alpha       = 0.85,
    show.legend = FALSE,
    na.rm       = TRUE
  ) +
  scale_fill_manual(values = my_cols, drop = FALSE) +
  scale_y_continuous(
    # fixed lower bound (0.25) zooms in on the data range; adjust if values fall
    # below 25%, otherwise points would be clipped
    labels = scales::percent_format(accuracy = 0.1),
    limits = c(0.25, y_max * 1.18)
  ) +
  labs(
    x = NULL,
    y = "Proportion of start/stop codon-disrupting homozygous variants"
  ) +
  box_theme +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1.1
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.8
    ),
    axis.ticks.y = element_line(
      colour = "black",
      linewidth = 0.8
    )
  )

p_hom_delclass_box

ggsave(
  # NOTE [repo path]: out_dir is an absolute path; set it to "figures" for the repo
  filename = file.path(out_dir, "Proportion of start stop codon-disrupting homozygous variants.png"),
  plot = p_hom_delclass_box,
  width = 6,
  height = 5,
  dpi = 600
)
