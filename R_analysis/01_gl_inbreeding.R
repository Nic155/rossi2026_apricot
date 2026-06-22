# =============================================================================
# Figure 1 (A, B, C): genetic load and inbreeding in cultivated vs wild apricot
#
# Compares deleterious load between the European cultivated (C1) and Central
# Asian wild (W2) gene pools, and relates load to individual inbreeding.
#
# Input:
#   01.RData, providing
#     - plot_counts : per-individual counts of derived deleterious variants,
#                     long format, split by class (Heterozygous / Homozygous)
#     - plot_ratio  : per-individual deleterious-to-synonymous (DEL/SYN) ratio,
#                     same long format and classes
#     - inbreeding      : per-individual inbreeding, incl. froh_800kb (F_ROH from
#                     runs of homozygosity >= 800 kb; see Methods for threshold)
#
# Output (PNG, to figures_final/):
#   gl.png          - Fig 1B: deleterious load (counts + DEL/SYN), Het vs Hom
#   froh.png        - Fig 1A: individual inbreeding F_ROH (ROH >= 800 kb)
#   correlation.png - Fig 1C: homozygous load vs F_ROH, counts and DEL/SYN
#
# Note: counts and the DEL/SYN ratio are complementary load measures. The ratio
# normalises for total variant number, so a lower C1 ratio reflects a smaller
# deleterious *fraction* rather than simply fewer variants overall.
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(purrr)
library(patchwork)


# NOTE [repo path]: absolute home-directory path. For the shared repo, replace
# with a repo-relative path, e.g. load("data/01.RData").
load("~/Paper parmeniaca/scripts/final/final final/rossi2026_apricot/01.RData")

# ── colour coding ─────────────────────────────────────────────────────────────
# C1 = cultivated (orange), W2 = wild (blue); kept consistent across all figures
group_cols <- c("C1" = "#D55E00", "W2" = "#0072B2")

# ── shared theme ──────────────────────────────────────────────────────────────
box_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position   = "none",
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background  = element_rect(fill = "white", color = "black"),
    strip.text        = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank()
  )

scatter_theme <- theme_bw(base_size = 12) +
  theme(
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background  = element_rect(fill = "white", color = "black"),
    strip.text        = element_text(size = 11),
    legend.position   = "bottom",
    legend.title      = element_blank()
  )

# ═══════════════════════════════════════════════════════════════════════════════
# i)  STACKED BOXPLOTS — counts and ratio for Het & Hom
# ═══════════════════════════════════════════════════════════════════════════════

# filter to Het and Hom only (drop Total/third level if present)
counts_hw <- plot_counts %>%
  filter(group %in% c("C1", "W2"),
         class %in% c("Heterozygous", "Homozygous")) %>%
  mutate(group = factor(group, levels = c("C1", "W2")),
         class = factor(class, levels = c("Heterozygous", "Homozygous")))

ratio_hw <- plot_ratio %>%
  filter(group %in% c("C1", "W2"),
         class %in% c("Heterozygous", "Homozygous")) %>%
  mutate(group = factor(group, levels = c("C1", "W2")),
         class = factor(class, levels = c("Heterozygous", "Homozygous")))

make_box <- function(df, ylab) {
  ggplot(df, aes(x = group, y = value, fill = group)) +
    # coef = Inf -> whiskers span the full data range (no points flagged as
    # outliers); raw points are overlaid by the jitter layer below
    geom_boxplot(width = 0.65, alpha = 0.85, outlier.shape = NA,
                 coef = Inf, na.rm = TRUE) +
    geom_jitter(
      color  = "grey20",    # fixed dark grey instead of aes(color = group)
      width  = 0.12,
      height = 0,
      size   = 1.8,
      alpha  = 0.85,        # slightly more opaque too
      show.legend = FALSE,
      na.rm  = TRUE
    ) +
    facet_wrap(~ class, nrow = 1, scales = "free_y") +
    scale_fill_manual(values  = group_cols, drop = FALSE) +
    scale_color_manual(values = group_cols, drop = FALSE) +
    scale_y_continuous(labels = comma) +
    labs(x = NULL, y = ylab) +
    box_theme
}

p_counts <- make_box(counts_hw, "Number of deleterious variants")
p_ratio  <- make_box(ratio_hw,  "DEL/SYN ratio")

# Fig 1B: counts (top row) over DEL/SYN ratio (bottom row), Het | Hom facets
p_load <- p_counts / p_ratio 

p_load

# NOTE [repo path]: replace absolute output path with e.g. "figures/gl.png"
ggsave("C:/Users/nirossi/Documents/Paper parmeniaca/figures_final/gl.png", p_load, width = 8, height = 6, dpi = 300)



# ── Fig 1A: individual inbreeding ─────────────────────────────────────────────
# FROH_800kb = fraction of the autosomal genome in ROH >= 800 kb. The 800 kb
# cutoff (~2 cM at ~2.5 cM/Mb) targets long ROH from recent inbreeding.
plot_froh <- inbreeding %>%
  filter(population %in% c("C1", "W2")) %>%
  transmute(
    sample,
    group      = factor(population, levels = c("C1", "W2")),
    FROH_800kb = froh_800kb
  ) %>%
  filter(!is.na(group), !is.na(FROH_800kb))

p_froh <- ggplot(plot_froh, aes(x = group, y = FROH_800kb, fill = group)) +
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
  scale_fill_manual(values = group_cols, drop = FALSE) +
  scale_y_continuous(labels = comma) +
  labs(
    x = NULL,
    y = "FROH_800kb"
  ) +
  box_theme

p_froh

# NOTE [repo path]: replace absolute output path with e.g. "figures/froh.png"
ggsave("C:/Users/nirossi/Documents/Paper parmeniaca/figures_final/froh.png", p_froh, width = 8, height = 6, dpi = 300)


# ═══════════════════════════════════════════════════════════════════════════════
# Homozygous load vs FROH_800kb only — counts and ratio stacked
# ═══════════════════════════════════════════════════════════════════════════════

# Pivot load to wide (one column per class) so the Homozygous column can be
# regressed against F_ROH per individual; then attach froh_800kb
# Fig 1C uses homozygous load only (recessive variants exposed by inbreeding).
wide_counts <- plot_counts %>%
  filter(group %in% c("C1", "W2")) %>%
  pivot_wider(names_from = class, values_from = value) %>%
  rename(sample = genotype, population = group) %>%
  left_join(
    inbreeding %>% select(sample, froh_800kb),
    by = "sample"
  ) %>%
  mutate(
    population = factor(population, levels = c("C1", "W2")),
    dataset = "counts"
  )

wide_ratio <- plot_ratio %>%
  filter(group %in% c("C1", "W2")) %>%
  pivot_wider(names_from = class, values_from = value) %>%
  rename(sample = genotype, population = group) %>%
  left_join(
    inbreeding %>% select(sample, froh_800kb),
    by = "sample"
  ) %>%
  mutate(
    population = factor(population, levels = c("C1", "W2")),
    dataset = "ratio"
  )

make_scatter <- function(df, xlab) {
  ggplot(
    df,
    aes(
      x = Homozygous,
      y = froh_800kb,
      color = population,
      fill = population
    )
  ) +
    geom_point(size = 1.8, alpha = 0.7, show.legend = FALSE) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.8,
      alpha = 0.15,
      show.legend = FALSE
    ) +
    scale_color_manual(values = group_cols, drop = FALSE) +
    scale_fill_manual(values  = group_cols, drop = FALSE) +
    scale_y_continuous(labels = comma) +
    labs(
      x = xlab,
      y = "FROH_800kb"
    ) +
    scatter_theme +
    theme(
      plot.title = element_blank(),
      legend.position = "none"
    )
}

p_scatter_counts <- make_scatter(
  wide_counts,
  "Number of homozygous deleterious variants"
)

p_scatter_ratio <- make_scatter(
  wide_ratio,
  "DEL/SYN ratio"
)

p_scatter_hom_froh <- p_scatter_counts / p_scatter_ratio +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "none",
    plot.title = element_blank()
  )

p_scatter_hom_froh



# NOTE [repo path]: replace absolute output path with e.g. "figures/correlation.png"
ggsave("C:/Users/nirossi/Documents/Paper parmeniaca/figures_final/correlation.png", p_scatter_hom_froh, width = 8, height = 9, dpi = 300)
