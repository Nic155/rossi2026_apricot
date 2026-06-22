# =============================================================================
# Chromosome-4 word cloud of W2-enriched deleterious-load gene symbols
#
# Visualises which genes drive the wild (W2) deleterious-burden enrichment within
# the major chromosome-4 selective-sweep hotspot, sized by the number of
# deleterious variants they carry.
#
# Input:
#   wordcloud.RData, providing
#     - go_variants       : W2-enriched GO-term deleterious variants with
#                           CHR, POS, gene_id (from 04.R)
#     - apricot_with_tair : gene_id -> SYMBOL mapping (TAIR/Arabidopsis symbols)
#
# Output (PNG, to figures_final/):
#   wordcloud_chr4.png
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggwordcloud)


# NOTE [repo path]: replace with a repo-relative path, e.g. "data/06.RData"
load("~/Paper parmeniaca/data_final/06.RData")

go_with_symbol <- go_variants %>%
  left_join(apricot_with_tair %>% select(gene_id, SYMBOL),
            by = "gene_id")

# See result
go_with_symbol %>% select(gene_id, term_id, SYMBOL) %>% head()


# subset chr4, split multi-symbol cells, drop ones starting with AT, count unique variants
# (SYMBOL can hold several ";"-separated symbols; AT-prefixed = raw Arabidopsis
#  locus IDs, dropped in favour of named gene symbols. n_distinct(POS) counts
#  unique variant positions per symbol, which sets word size.)
chr4_symbols <- go_with_symbol %>%
  filter(CHR == "chr4", !is.na(SYMBOL), nzchar(SYMBOL)) %>%
  separate_rows(SYMBOL, sep = ";") %>%
  mutate(SYMBOL = str_squish(SYMBOL)) %>%
  filter(nzchar(SYMBOL),
         !str_starts(toupper(SYMBOL), "AT")) %>%
  group_by(SYMBOL) %>%
  summarise(n_variants = n_distinct(POS), .groups = "drop") %>%
  arrange(desc(n_variants))

cat("Unique non-AT symbols on chr4:", nrow(chr4_symbols), "\n")
print(head(chr4_symbols, 20))

# distinct colours per word
palette_cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e",
                  "#e6ab02","#a6761d","#666666","#1f78b4","#33a02c",
                  "#fb9a99","#fdbf6f","#cab2d6","#b15928","#8da0cb")

chr4_symbols <- chr4_symbols %>%
  mutate(colour = rep_len(palette_cols, n()),
         # floor word size at 10 so low-count symbols stay legible in the cloud
         n_variants_plot = pmax(n_variants, 10))

set.seed(42)   # fixes the word-cloud layout (reproducible placement)
p=ggplot(chr4_symbols,
       aes(label = SYMBOL, size = n_variants_plot, color = colour)) +
  geom_text_wordcloud_area(eccentricity = 0.55, rm_outside = TRUE) +
  scale_size_area(max_size = 28) +
  scale_color_identity() +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) 
p

# NOTE [repo path]: replace absolute output path with e.g. "figures/wordcloud_chr4.png"
ggsave("C:/Users/nirossi/Documents/Paper parmeniaca/figures_final/wordcloud_chr4.png", p, width = 8, height = 6, dpi = 300)



