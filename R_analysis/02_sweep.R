# =============================================================================
# Figure 2A: deleterious burden in selective-sweep vs background coding regions
#
# Contrasts per-individual deleterious dosage in swept coding sequence against
# all other (background/control) coding sequence, separately for the cultivated
# (C1) and wild (W2) gene pools. Tests whether sweeps carry reduced load
# (linked purifying selection) rather than inflated load (hitchhiking).
#
# Input:
#   sweep.RData, providing
#     - burden_long : long table of per-individual deleterious dosage, already
#                     expressed per 100 kb of CDS, with columns
#                     group (C1/W2) and region (Control / Sweep).
#                     Sweep windows are the top-5% (95th-percentile) windows of
#                     the wild/cultivated diversity ratio pi_W2/pi_C1 from
#                     Groppi et al. 2021 (computed upstream, not in this script).
#
# Output (PNG, to figures_final/):
#   sweep.png - Fig 2A: Control vs Sweep dosage, faceted by group
# =============================================================================

library(tidyverse)
library(ggplot2)
library(scales)


# =============================================================================
# Load data
# =============================================================================
# NOTE [repo path]: replace with a repo-relative path, e.g. "data/02.RData"
load("~/Paper parmeniaca/scripts/final/final final/rossi2026_apricot/02.RData")

# =============================================================================
# Colours and levels
# =============================================================================
matched_group_levels <- c("C1", "W2")
group_cols <- c("C1" = "#D55E00", "W2" = "#0072B2")

# =============================================================================
# Shared theme
# =============================================================================
shared_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position    = "none",
    panel.spacing.x    = grid::unit(1, "lines"),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background   = element_rect(fill = "white", color = "black"),
    strip.text         = element_text(size = 13),
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    panel.grid.major.x = element_blank()
  )

# =============================================================================
# 1. Sweep vs Control
# =============================================================================
# Order region so Control plots before Sweep within each group facet
burden_long <- burden_long %>%
  filter(group %in% matched_group_levels) %>%
  mutate(
    group  = factor(group,  levels = matched_group_levels),
    region = factor(region, levels = c("Control", "Sweep"))
  )

# ── boxplot ───────────────────────────────────────────────────────────────────
p_sweep_box <- ggplot(
  burden_long,
  aes(x = region, y = burden_per_100kb_CDS, fill = group)
) +
  geom_boxplot(
    width        = 0.55,
    alpha        = 0.85,
    outlier.shape = NA,
    color        = "black",       # black frame
    coef         = Inf,           # whiskers span full range; points overlaid below
    na.rm        = TRUE
  ) +
  geom_jitter(
    color  = "grey20",    # fixed dark grey instead of aes(color = group)
    width  = 0.12,
    height = 0,
    size   = 1.8,
    alpha  = 0.85,        # slightly more opaque too
    show.legend = FALSE,
    na.rm  = TRUE
  ) +
  facet_wrap(~ group, scales = "fixed") +
  scale_fill_manual(values  = group_cols, drop = FALSE) +
  scale_color_manual(values = group_cols, drop = FALSE) +
  scale_y_continuous(labels = comma) +
  labs(
    x = NULL,
    y = "Deleterious dosage per 100 kb CDS"
  ) +
  shared_theme

p_sweep_box


# NOTE [repo path]: replace absolute output path with e.g. "figures/sweep.png"
ggsave("C:/Users/nirossi/Documents/Paper parmeniaca/figures_final/sweep.png", p_sweep_box, width = 8, height = 6, dpi = 300)

