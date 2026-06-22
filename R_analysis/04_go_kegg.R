# =============================================================================
# GO / KEGG over-representation analysis of differentially loaded genes
# in wild (W2) vs cultivated (C1) Prunus armeniaca
#
# Produces the gene-level functional enrichment behind Figure 2B, plus the
# variant-position tables (go_variants / kegg_variants) consumed by the circos
# script (05.R).
#
# Input:
#   go_kegg_new.RData, providing
#     - df1           : deleterious variants with annotations (IMPACT, EFFECT,
#                       DEL_CLASS, PROVEAN_SCORE, geneID) and one dosage column
#                       per individual, suffixed _C1 / _W2
#     - genes2Go      : gene_id -> GO term/description annotation table
#     - KEGG_pathways : gene_id -> KEGG pathway annotation table
#     - gene_coords   : gene_id -> chr/start/end (for the tandem-collapse check)
#
# Output:
#   W2_enrichment_GO_KEGG.png - Fig 2B (W2 GO + KEGG fold-enrichment barplots)
#   go_variants / kegg_variants - in-memory position tables for 05.R
#
# Workflow:
#   - background = all genes carrying a candidate deleterious variant
#   - per-gene log2FC = log2(mean_W2_dosage + 0.001) - log2(mean_C1_dosage + 0.001)
#   - per-gene Wilcoxon test of per-individual dosage, BH-adjusted across genes
#   - selection  = Wilcoxon p_adj < 0.05 AND |log2FC| > 1
#                  (W2-enriched: log2FC > 1; C1-enriched: log2FC < -1)
#   - ORA via clusterProfiler::enricher (hypergeometric, BH), term size 5-500
#   - KEGG pathways filtered to plant-relevant categories via KEGG BRITE hierarchy
#   - PARALLEL ROBUSTNESS CHECK: enrichment re-run on tandem-collapsed loci
#     (genes within 2 kb on the same chromosome merged into one locus; see
#     `gap` below) to control for non-independence of clustered gene families
#     (terpene synthases, NB-LRRs).
#   - Visualisation: horizontal barplots (fold enrichment), printed to R.
# =============================================================================

# NOTE [repo path]: replace with a repo-relative path, e.g. "data/go_kegg_new.RData"
load("~/Paper parmeniaca/scripts/final/final final/rossi2026_apricot/04.RData")   # df1, genes2Go, KEGG_pathways, gene_coords

library(tidyverse)
library(clusterProfiler)
library(forcats)
library(patchwork)

# -----------------------------------------------------------------------------
# Run ONCE on a machine with internet to cache the KEGG BRITE hierarchy:
# hierarchy_lines <- readLines("https://rest.kegg.jp/get/br:ko00001", warn = FALSE)
# saveRDS(hierarchy_lines, "~/Paper parmeniaca/data_final/kegg_ko00001_hierarchy.rds")
# -----------------------------------------------------------------------------

# =============================================================================
# 1. Deleterious variants (PROVEAN / HIGH) with a geneID
# =============================================================================
# Candidate deleterious variants = below-threshold PROVEAN substitutions plus
# HIGH-impact SnpEff variants. The genes carrying these define the ENRICHMENT
# BACKGROUND (universe) used throughout.
del_variants <- df1 %>%
  filter(IMPACT %in% c("PROVEAN", "HIGH"), !is.na(geneID))

cat("Total PROVEAN/HIGH variants with geneID:", nrow(del_variants), "\n")
cat("Unique genes:", n_distinct(del_variants$geneID), "\n")

# Per-individual dosage columns are suffixed by gene pool (_C1 / _W2)
c1_cols <- grep("_C1$", colnames(del_variants), value = TRUE)
w2_cols <- grep("_W2$", colnames(del_variants), value = TRUE)
cat("C1 samples:", length(c1_cols), "   W2 samples:", length(w2_cols), "\n")

# =============================================================================
# 2. Per-gene dosage and log2FC
# =============================================================================
# log2FC summarises direction/magnitude of the W2-vs-C1 dosage difference.
# The +0.001 pseudocount avoids log(0) when a group's mean dosage is zero.
gene_burden <- del_variants %>%
  group_by(geneID) %>%
  summarise(
    n_variants = n(),
    mean_dosage_C1 = mean(rowMeans(across(all_of(c1_cols)), na.rm = TRUE), na.rm = TRUE),
    mean_dosage_W2 = mean(rowMeans(across(all_of(w2_cols)), na.rm = TRUE), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    log2FC = log2(mean_dosage_W2 + 0.001) - log2(mean_dosage_C1 + 0.001),
    # NB: this ±0.5 `direction` label is descriptive only and is NOT the
    # selection cutoff. Gene-set selection downstream uses |log2FC| > 1 (which
    # matches the manuscript). See FLAG in the accompanying notes.
    direction = case_when(
      log2FC >  0.5 ~ "W2_enriched",
      log2FC < -0.5 ~ "C1_enriched",
      TRUE          ~ "similar"
    )
  )

# =============================================================================
# 3. Per-gene Wilcoxon test (per-individual total dosage, group comparison)
# =============================================================================
# Per gene: sum each individual's dosage across that gene's variants, then
# Wilcoxon rank-sum test of W2 vs C1 individuals. p-values are BH-adjusted
# across all tested genes (p_adj) to control the false discovery rate.
gene_wilcox <- del_variants %>%
  group_by(geneID) %>%
  summarise(across(all_of(c(c1_cols, w2_cols)), ~ sum(.x, na.rm = TRUE)),
            n_variants = n(), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    p_wilcox = tryCatch(
      wilcox.test(c_across(all_of(w2_cols)), c_across(all_of(c1_cols)),
                  exact = FALSE)$p.value,
      error = function(e) NA_real_)
  ) %>%
  ungroup() %>%
  select(geneID, n_variants, p_wilcox) %>%
  left_join(gene_burden %>% select(geneID, mean_dosage_C1, mean_dosage_W2, log2FC, direction),
            by = "geneID") %>%
  mutate(p_adj = p.adjust(p_wilcox, method = "BH"))

cat("\nGenes significant (p_adj < 0.05):", sum(gene_wilcox$p_adj < 0.05, na.rm = TRUE), "\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

normalize_gene_id <- function(gene_id) {
  gene_id <- as.character(gene_id); gene_id <- trimws(gene_id)
  gene_id[gene_id %in% c("", "NA", "N/A", "NULL", "NaN")] <- NA_character_
  gene_id <- sub("\\.t[0-9]+.*$", "", gene_id, perl = TRUE)
  gene_id <- sub("\\.p[0-9]+.*$", "", gene_id, perl = TRUE)
  gene_id <- sub("\\.mRNA[0-9A-Za-z_.-]*$", "", gene_id, ignore.case = TRUE, perl = TRUE)
  gene_id <- sub("_v[0-9]+(\\.[0-9]+)*$", "", gene_id, perl = TRUE)
  gene_id
}

clean_character <- function(x) {
  x <- as.character(x); x <- trimws(x)
  x[toupper(x) %in% c("", "NA", "N/A", "NULL", "NONE", "NAN")] <- NA_character_
  x
}

first_non_missing <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

clean_annotation_table <- function(annotation_table, table_name = "annotation table") {
  annotation_table <- as.data.frame(annotation_table, stringsAsFactors = FALSE)
  if (ncol(annotation_table) < 3) stop(table_name, " must have >= 3 columns.", call. = FALSE)
  
  cleaned <- annotation_table[, 1:3, drop = FALSE]
  names(cleaned) <- c("query_raw", "term_id", "description")
  cleaned$query_raw   <- clean_character(cleaned$query_raw)
  cleaned$term_id     <- clean_character(cleaned$term_id)
  cleaned$description <- clean_character(cleaned$description)
  cleaned$gene_id     <- normalize_gene_id(cleaned$query_raw)
  
  all_empty <- is.na(cleaned$query_raw) & is.na(cleaned$term_id) & is.na(cleaned$description)
  
  query_header_values <- c("query","query id","query_id","gene","gene id","gene_id",
                           "protein","protein id","protein_id","transcript","transcript id")
  term_header_values  <- c("match","term","term id","term_id","go","go id","go_id",
                           "pathway","pathway id","pathway_id","kegg pathway",
                           "ortholog","ortholog id","ipr","ipr id")
  description_header_values <- c("description","name","annotation","term description")
  
  combined_text <- paste(cleaned$query_raw, cleaned$term_id, cleaned$description)
  header_like <- tolower(cleaned$query_raw) %in% query_header_values |
    tolower(cleaned$term_id) %in% term_header_values |
    tolower(cleaned$description) %in% description_header_values |
    grepl("may contain multiple worksheets", combined_text, ignore.case = TRUE)
  
  keep_rows <- !all_empty & !header_like &
    !is.na(cleaned$gene_id) & nzchar(cleaned$gene_id) &
    !is.na(cleaned$term_id) & nzchar(cleaned$term_id)
  
  cleaned <- cleaned[keep_rows, , drop = FALSE]
  cleaned$description[is.na(cleaned$description) | !nzchar(cleaned$description)] <-
    cleaned$term_id[is.na(cleaned$description) | !nzchar(cleaned$description)]
  
  pair_key <- paste(cleaned$gene_id, cleaned$term_id, sep = "\r")
  cleaned <- cleaned[!duplicated(pair_key), , drop = FALSE]
  rownames(cleaned) <- NULL
  cleaned
}

parse_go_descriptions <- function(go_annotation) {
  go_annotation$ontology <- NA_character_
  go_annotation$clean_description <- go_annotation$description
  ontology_prefixes <- c("Biological Process", "Molecular Function", "Cellular Component")
  for (prefix in ontology_prefixes) {
    pattern <- paste0("^", prefix, "\\s*:\\s*")
    matched <- grepl(pattern, go_annotation$description, ignore.case = TRUE)
    go_annotation$ontology[matched] <- prefix
    go_annotation$clean_description[matched] <-
      sub(pattern, "", go_annotation$description[matched], ignore.case = TRUE)
  }
  go_annotation$description <- go_annotation$clean_description
  go_annotation
}

fetch_kegg_pathway_hierarchy <- function() {
  # Prefer a locally cached BRITE hierarchy (see the cache-once block near the
  # top); fall back to a live KEGG REST call, returning NULL if offline.
  # NOTE [repo path]: cache file is an absolute path; make it repo-relative,
  # e.g. "data/kegg_ko00001_hierarchy.rds".
  local_path <- path.expand("~/Paper parmeniaca/data_final/kegg_ko00001_hierarchy.rds")
  if (file.exists(local_path)) return(readRDS(local_path))
  tryCatch(readLines("https://rest.kegg.jp/get/br:ko00001", warn = FALSE),
           error = function(e) NULL)
}

parse_kegg_pathway_hierarchy <- function(hierarchy_lines) {
  if (is.null(hierarchy_lines) || length(hierarchy_lines) == 0) {
    return(data.frame(term_id=character(), kegg_top_category=character(),
                      kegg_subcategory=character(), stringsAsFactors=FALSE))
  }
  current_top_category <- NA_character_; current_subcategory <- NA_character_
  records <- list()
  strip_kegg_markup <- function(x) trimws(gsub("<[^>]+>", "", x))
  for (line in hierarchy_lines) {
    line <- strip_kegg_markup(line)
    if (grepl("^A", line)) { current_top_category <- trimws(sub("^A[0-9]*\\s*", "", line)); current_subcategory <- NA_character_; next }
    if (grepl("^B", line)) { current_subcategory <- trimws(sub("^B[0-9]*\\s*", "", line)); next }
    if (grepl("^C", line)) {
      pathway_id <- NA_character_
      path_match <- regmatches(line, regexpr("\\[PATH:(ko|map)[0-9]{5}\\]", line))
      if (length(path_match) > 0 && nzchar(path_match)) {
        pathway_id <- sub("^\\[PATH:(ko|map)([0-9]{5})\\]$", "ko\\2", path_match)
      } else {
        numeric_match <- regmatches(line, regexpr("[0-9]{5}", line))
        if (length(numeric_match) > 0 && nzchar(numeric_match)) pathway_id <- paste0("ko", numeric_match)
      }
      if (!is.na(pathway_id)) {
        records[[length(records)+1]] <- data.frame(term_id=pathway_id,
                                                   kegg_top_category=current_top_category, kegg_subcategory=current_subcategory,
                                                   stringsAsFactors=FALSE)
      }
    }
  }
  if (length(records) == 0) {
    return(data.frame(term_id=character(), kegg_top_category=character(),
                      kegg_subcategory=character(), stringsAsFactors=FALSE))
  }
  hierarchy <- do.call(rbind, records)
  hierarchy <- hierarchy[!duplicated(hierarchy$term_id), , drop = FALSE]
  rownames(hierarchy) <- NULL
  hierarchy
}

filter_kegg_pathway_annotations <- function(kegg_annotations) {
  # Restrict KEGG to plant-relevant pathways. KEGG's reference pathways include
  # many animal/disease terms that are meaningless for apricot; these are dropped
  # in two ways: (1) by BRITE top-level category (Human Diseases, Drug
  # Development, Organismal Systems) when the hierarchy is available, and
  # (2) by a conservative keyword regex on the description as a fallback/backup.
  pathway_descriptions <- stats::aggregate(description ~ term_id, data = kegg_annotations,
                                           FUN = first_non_missing)
  hierarchy <- parse_kegg_pathway_hierarchy(fetch_kegg_pathway_hierarchy())
  if (nrow(hierarchy) > 0) {
    pathway_descriptions <- merge(pathway_descriptions, hierarchy, by = "term_id", all.x = TRUE)
    filter_method <- "KEGG pathway hierarchy plus conservative description regex"
  } else {
    pathway_descriptions$kegg_top_category <- NA_character_
    pathway_descriptions$kegg_subcategory  <- NA_character_
    filter_method <- "conservative description regex only; KEGG hierarchy unavailable"
  }
  
  excluded_top_categories <- c("Human Diseases", "Drug Development", "Organismal Systems")
  disease_or_nonplant_regex <- paste(c(
    "Measles","Legionellosis","Toxoplasmosis","Alzheimer","Parkinson","Huntington",
    "cancer","carcinoma","leukemia","lymphoma","addiction","alcoholism","amphetamine",
    "cocaine","morphine","nicotine","tuberculosis","influenza","hepatitis","HIV",
    "malaria","salmonella","shigellosis","pertussis","leishmaniasis","Chagas",
    "rheumatoid","asthma","diabetes","cardiomyopathy","viral myocarditis",
    "human immunodeficiency virus","herpes","papillomavirus","coronavirus",
    "Epstein-Barr","Kaposi","prion","amyotrophic","multiple sclerosis",
    "inflammatory bowel disease","systemic lupus","allograft","graft-versus-host",
    "Staphylococcus aureus infection","Pathogenic Escherichia coli infection",
    "Yersinia infection","Vibrio cholerae infection",
    "Epithelial cell signaling in Helicobacter","African trypanosomiasis","amoebiasis",
    "estrogen signaling","GnRH signaling","oxytocin signaling","relaxin signaling",
    "prolactin signaling","thyroid hormone","insulin secretion","glucagon signaling",
    "adipocytokine","PPAR signaling","renin","aldosterone","adrenergic","dopaminergic",
    "serotonergic","cholinergic","GABAergic","glutamatergic","olfactory",
    "taste transduction","salivary secretion","pancreatic secretion","gastric acid",
    "bile secretion","platelet activation","complement and coagulation","Fc gamma",
    "Fc epsilon","T cell receptor","B cell receptor","natural killer","Th17","IL-17",
    "NOD-like receptor","Toll-like receptor","C-type lectin receptor","chemokine signaling",
    "leukocyte","hematopoietic","osteoclast","cardiac muscle","vascular smooth muscle",
    "long-term potentiation","long-term depression","circadian entrainment"
  ), collapse = "|")
  
  removed_by_category <- !is.na(pathway_descriptions$kegg_top_category) &
    pathway_descriptions$kegg_top_category %in% excluded_top_categories
  removed_by_regex <- grepl(disease_or_nonplant_regex, pathway_descriptions$description,
                            ignore.case = TRUE)
  pathway_descriptions$remove_from_main_kegg <- removed_by_category | removed_by_regex
  
  kept_term_ids <- pathway_descriptions$term_id[!pathway_descriptions$remove_from_main_kegg]
  filtered_annotations <- kegg_annotations[kegg_annotations$term_id %in% kept_term_ids, , drop = FALSE]
  
  list(annotations = filtered_annotations, filter_method = filter_method)
}

count_significant <- function(result, fdr_threshold) {
  sum(!is.na(result$adjusted_p_value) & result$adjusted_p_value < fdr_threshold)
}

# clusterProfiler::enricher wrapper -> tidy data.frame (carries ontology + fold enrichment)
# Run with pvalueCutoff/qvalueCutoff = 1 so NO term is dropped inside enricher;
# all terms are returned and filtered later at adjusted_p_value < fdr_threshold.
# Term sizes restricted to 5-500 genes (minGSSize/maxGSSize). fold_enrichment is
# the gene ratio (selected) divided by the background ratio.
run_enricher <- function(annotation, selected_genes, background_genes, group_label) {
  t2g <- annotation %>% select(term_id, gene_id) %>% distinct()
  t2n <- annotation %>% distinct(term_id, description)
  ont_lookup <- if ("ontology" %in% names(annotation)) {
    annotation %>% distinct(term_id, ontology)
  } else tibble(term_id = unique(annotation$term_id), ontology = NA_character_)
  
  res <- enricher(
    gene = unique(selected_genes), universe = unique(background_genes),
    TERM2GENE = t2g, TERM2NAME = t2n,
    pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
    minGSSize = 5, maxGSSize = 500
  )
  
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    return(data.frame(term_id=character(), description=character(), ontology=character(),
                      selected_gene_count=integer(), background_gene_count=integer(),
                      gene_ratio=numeric(), fold_enrichment=numeric(),
                      p_value=numeric(), adjusted_p_value=numeric(),
                      group=character(), stringsAsFactors=FALSE))
  }
  
  as.data.frame(res) %>%
    separate(GeneRatio, into = c("gr_n","gr_d"), sep = "/", convert = TRUE, remove = FALSE) %>%
    separate(BgRatio,   into = c("bg_n","bg_d"), sep = "/", convert = TRUE, remove = FALSE) %>%
    transmute(
      term_id               = ID,
      description           = Description,
      selected_gene_count   = Count,
      background_gene_count  = bg_n,
      gene_ratio            = gr_n / gr_d,
      fold_enrichment       = (gr_n / gr_d) / (bg_n / bg_d),
      p_value               = pvalue,
      adjusted_p_value      = p.adjust,
      group                 = group_label
    ) %>%
    left_join(ont_lookup, by = "term_id")
}

cat("Helper functions defined.\n")

# settings
fdr_threshold <- 0.05

# =============================================================================
# Selected gene sets (significance + effect size) and background
# =============================================================================
# Foreground = significant (BH p_adj < 0.05) genes with a large directional
# effect (|log2FC| > 1). Split by direction: W2-enriched (log2FC > 1) and
# C1-enriched (log2FC < -1). normalize_gene_id strips transcript/isoform suffixes
# so variant geneIDs match the annotation/coordinate tables.
results_final_W2_enriched <- gene_wilcox %>%
  filter(p_adj < 0.05, log2FC >  1) %>%
  mutate(gene_id = normalize_gene_id(geneID)) %>%
  filter(!is.na(gene_id), nzchar(gene_id))

results_final_C1_enriched <- gene_wilcox %>%
  filter(p_adj < 0.05, log2FC < -1) %>%
  mutate(gene_id = normalize_gene_id(geneID)) %>%
  filter(!is.na(gene_id), nzchar(gene_id))

# Background universe = every gene carrying a candidate deleterious variant
tested_genes <- del_variants %>%
  filter(!is.na(geneID)) %>%
  group_by(gene_id = normalize_gene_id(geneID)) %>%
  summarise(n_deleterious_variants = n(), .groups = "drop") %>%
  filter(!is.na(gene_id), nzchar(gene_id))

# =============================================================================
# Clean annotations
# =============================================================================
go_annotations           <- parse_go_descriptions(clean_annotation_table(genes2Go, "genes2Go"))
kegg_pathway_unfiltered  <- clean_annotation_table(KEGG_pathways, "KEGG_pathways")
kegg_pathway_filter      <- filter_kegg_pathway_annotations(kegg_pathway_unfiltered)
kegg_pathway_annotations <- kegg_pathway_filter$annotations
cat("KEGG filter method:", kegg_pathway_filter$filter_method, "\n")

# =============================================================================
# Gene-level gene sets / background vectors
# =============================================================================
c1_genes <- unique(na.omit(normalize_gene_id(results_final_C1_enriched$gene_id)))
w2_genes <- unique(na.omit(normalize_gene_id(results_final_W2_enriched$gene_id)))
c1_genes <- c1_genes[nzchar(c1_genes)]
w2_genes <- w2_genes[nzchar(w2_genes)]
background_genes <- unique(na.omit(normalize_gene_id(tested_genes$gene_id)))
background_genes <- background_genes[nzchar(background_genes)]

cat("\n[GENE LEVEL] C1:", length(c1_genes), " W2:", length(w2_genes),
    " background:", length(background_genes), "\n")


# =============================================================================
# FOREGROUND VARIANT-CLASS SWITCH
#   Choose which variant class defines the between-group foreground.
#   BACKGROUND (background_genes, ~22k) was built above from
#   IMPACT %in% c("PROVEAN","HIGH") and is NOT modified here.
#     "provean_high" -> original foreground (IMPACT PROVEAN/HIGH contrast)
#     "severe"       -> foreground from DEL_CLASS %in% severe_classes
# =============================================================================

# per-gene Wilcoxon + log2FC on whatever variant subset you pass (mirrors the
# original gene_wilcox/gene_burden logic exactly, just on fewer rows)
compute_gene_wilcox <- function(variant_df) {
  c1c <- grep("_C1$", colnames(variant_df), value = TRUE)
  w2c <- grep("_W2$", colnames(variant_df), value = TRUE)
  
  burden <- variant_df %>%
    group_by(geneID) %>%
    summarise(
      mean_dosage_C1 = mean(rowMeans(across(all_of(c1c)), na.rm = TRUE), na.rm = TRUE),
      mean_dosage_W2 = mean(rowMeans(across(all_of(w2c)), na.rm = TRUE), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(log2FC = log2(mean_dosage_W2 + 0.001) - log2(mean_dosage_C1 + 0.001))
  
  variant_df %>%
    group_by(geneID) %>%
    summarise(across(all_of(c(c1c, w2c)), ~ sum(.x, na.rm = TRUE)),
              n_variants = n(), .groups = "drop") %>%
    rowwise() %>%
    mutate(p_wilcox = tryCatch(
      wilcox.test(c_across(all_of(w2c)), c_across(all_of(c1c)), exact = FALSE)$p.value,
      error = function(e) NA_real_)) %>%
    ungroup() %>%
    select(geneID, n_variants, p_wilcox) %>%
    left_join(burden, by = "geneID") %>%
    mutate(p_adj = p.adjust(p_wilcox, method = "BH"))
}

foreground_variant_set <- "severe"                      # "provean_high" or "severe"
severe_classes        <- c("HIGH", "STOP_START_LOSS")   # <- match your DEL_CLASS labels

if (foreground_variant_set == "severe") {
  cat("\nDEL_CLASS values in df1 (confirm these match severe_classes):\n")
  print(table(df1$DEL_CLASS, useNA = "ifany"))
  
  del_variants_fg <- df1 %>%
    filter(toupper(trimws(DEL_CLASS)) %in% toupper(severe_classes), !is.na(geneID))
  
  cat("\n[FOREGROUND = severe] variants:", nrow(del_variants_fg),
      "  genes:", n_distinct(del_variants_fg$geneID), "\n")
  stopifnot(nrow(del_variants_fg) > 0)
  
  fg_wilcox <- compute_gene_wilcox(del_variants_fg)
} else {
  fg_wilcox <- gene_wilcox                              # original PROVEAN/HIGH contrast
}

# Reselect W2/C1 foreground genes from the chosen variant class (same
# p_adj < 0.05 & |log2FC| > 1 cutoffs). NOTE: w2_genes / c1_genes are
# overwritten AGAIN by the PROVEAN_SCORE < -9 block below, which is what
# actually feeds the enrichment. The values set here are intermediate.
w2_genes <- fg_wilcox %>%
  filter(p_adj < 0.05, log2FC >  1) %>%
  transmute(gene_id = normalize_gene_id(geneID)) %>%
  filter(!is.na(gene_id), nzchar(gene_id)) %>%
  pull(gene_id) %>% unique()

c1_genes <- fg_wilcox %>%
  filter(p_adj < 0.05, log2FC < -1) %>%
  transmute(gene_id = normalize_gene_id(geneID)) %>%
  filter(!is.na(gene_id), nzchar(gene_id)) %>%
  pull(gene_id) %>% unique()

cat("\n[FOREGROUND =", foreground_variant_set, "]  W2:", length(w2_genes),
    "  C1:", length(c1_genes),
    "   (background unchanged:", length(background_genes), "genes)\n")
cat("  foreground inside the 22k background  ->  W2:",
    sum(w2_genes %in% background_genes), "/", length(w2_genes),
    "   C1:", sum(c1_genes %in% background_genes), "/", length(c1_genes), "\n")
cat("  (a foreground gene NOT in the background is dropped by enricher, by design)\n")




# =============================================================================
# FOREGROUND = variants with PROVEAN_SCORE < -9
# This overwrites w2_genes and c1_genes for the downstream enrichment.
# BACKGROUND remains unchanged: background_genes from IMPACT %in% c("PROVEAN","HIGH")
#
# FLAG [differs from Methods]: the manuscript describes the per-gene burden test
# on PROVEAN/HIGH deleterious variants selected at |log2FC| > 1; it does not
# mention a PROVEAN_SCORE < -9 (strong-effect) foreground. This block, run last,
# is what actually defines the W2/C1 gene sets fed to the enrichment below.
# Reconcile with the Methods or document the threshold before release.
# =============================================================================

provean_score_threshold <- -9

del_variants_fg_provean9 <- df1 %>%
  filter(
    !is.na(geneID),
    !is.na(PROVEAN_SCORE),
    PROVEAN_SCORE < provean_score_threshold
  )

cat("\n[FOREGROUND = PROVEAN_SCORE <", provean_score_threshold, "]\n")
cat("  variants selected:", nrow(del_variants_fg_provean9), "\n")
cat("  genes before Wilcoxon/group selection:",
    n_distinct(normalize_gene_id(del_variants_fg_provean9$geneID)), "\n")

stopifnot(nrow(del_variants_fg_provean9) > 0)

fg_wilcox_provean9 <- compute_gene_wilcox(del_variants_fg_provean9) %>%
  mutate(gene_id = normalize_gene_id(geneID)) %>%
  filter(!is.na(gene_id), nzchar(gene_id))

cat("  genes after grouping variants by gene:",
    n_distinct(fg_wilcox_provean9$gene_id), "\n")

w2_genes <- fg_wilcox_provean9 %>%
  filter(p_adj < 0.05, log2FC > 1) %>%
  pull(gene_id) %>%
  unique()

c1_genes <- fg_wilcox_provean9 %>%
  filter(p_adj < 0.05, log2FC < -1) %>%
  pull(gene_id) %>%
  unique()

foreground_background_sizes_provean9 <- data.frame(
  set = c(
    "W2 foreground genes",
    "C1 foreground genes",
    "All foreground genes",
    "Background genes"
  ),
  n_genes = c(
    length(w2_genes),
    length(c1_genes),
    length(unique(c(w2_genes, c1_genes))),
    length(background_genes)
  ),
  n_inside_background = c(
    sum(w2_genes %in% background_genes),
    sum(c1_genes %in% background_genes),
    sum(unique(c(w2_genes, c1_genes)) %in% background_genes),
    length(background_genes)
  )
)

print(foreground_background_sizes_provean9)

cat("\nGene-level counts after grouping variants by gene, PROVEAN_SCORE <",
    provean_score_threshold, ":\n")
cat("  W2 foreground genes :", length(w2_genes), "\n")
cat("  C1 foreground genes :", length(c1_genes), "\n")
cat("  All foreground genes:", length(unique(c(w2_genes, c1_genes))), "\n")
cat("  Background genes    :", length(background_genes), "\n")
cat("  W2 inside background:", sum(w2_genes %in% background_genes), "/",
    length(w2_genes), "\n")
cat("  C1 inside background:", sum(c1_genes %in% background_genes), "/",
    length(c1_genes), "\n")





# =============================================================================
# (A) GENE-LEVEL ENRICHMENT  --  headline result
# =============================================================================
# ORA of each foreground set (C1/W2 enriched genes) against the deleterious-gene
# background, for GO and KEGG. W2 GO + KEGG are the panels shown in Fig 2B.
c1_go_enrichment   <- run_enricher(go_annotations,          c1_genes, background_genes, "C1")
w2_go_enrichment   <- run_enricher(go_annotations,          w2_genes, background_genes, "W2")
c1_kegg_enrichment <- run_enricher(kegg_pathway_annotations, c1_genes, background_genes, "C1")
w2_kegg_enrichment <- run_enricher(kegg_pathway_annotations, w2_genes, background_genes, "W2")

cat("\n=== GENE-LEVEL enrichment (FDR <", fdr_threshold, ") ===\n")
cat("  C1 GO  :", count_significant(c1_go_enrichment,  fdr_threshold), "\n")
cat("  W2 GO  :", count_significant(w2_go_enrichment,  fdr_threshold), "\n")
cat("  C1 KEGG:", count_significant(c1_kegg_enrichment, fdr_threshold), "\n")
cat("  W2 KEGG:", count_significant(w2_kegg_enrichment, fdr_threshold), "\n")

# =============================================================================
# (B) LOCUS-LEVEL ENRICHMENT  --  robustness check for tandem clustering
#     gene_coords (chr,start,end,gene_id) must be in the loaded RData
# =============================================================================
# Tandemly duplicated families (terpene synthases, NB-LRRs) sit in physical
# clusters, so counting each paralogue as independent could inflate enrichment.
# As a control, collapse genes within `gap` bp on the same chromosome into one
# locus and re-run the ORA on loci. Enrichment surviving this is robust to
# tandem-array structure.
gene_coords <- gene_coords %>%
  mutate(gene_id = normalize_gene_id(gene_id),
         start = as.integer(start), end = as.integer(end)) %>%
  filter(!is.na(gene_id), nzchar(gene_id)) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  arrange(chr, start)

gap <- 2000  # 2 kb collapse window (matches Methods); try larger values to confirm stability

# Within each chromosome, start a new cluster whenever the gap to the previous
# gene's start exceeds `gap`; genes in the same run share a locus_id.
cluster_map <- gene_coords %>%
  group_by(chr) %>%
  mutate(new_cluster = cumsum(c(TRUE, diff(start) > gap))) %>%
  ungroup() %>%
  mutate(locus_id = paste(chr, new_cluster, sep = "_")) %>%
  select(gene_id, locus_id)

g2l <- setNames(cluster_map$locus_id, cluster_map$gene_id)  # gene_id -> locus_id lookup
to_loci <- function(genes) unique(na.omit(g2l[genes]))      # map a gene set to its loci

cat("\n[LOCUS LEVEL] genes:", nrow(cluster_map),
    " -> loci:", n_distinct(cluster_map$locus_id), " (gap =", gap, "bp)\n")

# remap annotations gene -> locus (a locus carries a term if any of its genes does)
remap_annotation_to_loci <- function(annotation) {
  annotation %>%
    mutate(locus_id = g2l[gene_id]) %>%
    filter(!is.na(locus_id)) %>%
    transmute(term_id, gene_id = locus_id, description,
              ontology = if ("ontology" %in% names(.)) ontology else NA_character_) %>%
    distinct()
}

go_annotations_loci   <- remap_annotation_to_loci(go_annotations)
kegg_annotations_loci <- remap_annotation_to_loci(kegg_pathway_annotations)

c1_loci         <- to_loci(c1_genes)
w2_loci         <- to_loci(w2_genes)
background_loci <- to_loci(background_genes)

c1_go_loci   <- run_enricher(go_annotations_loci,   c1_loci, background_loci, "C1")
w2_go_loci   <- run_enricher(go_annotations_loci,   w2_loci, background_loci, "W2")
c1_kegg_loci <- run_enricher(kegg_annotations_loci, c1_loci, background_loci, "C1")
w2_kegg_loci <- run_enricher(kegg_annotations_loci, w2_loci, background_loci, "W2")

cat("\n=== LOCUS-LEVEL enrichment after tandem collapse (FDR <", fdr_threshold, ") ===\n")
cat("  C1 GO  :", count_significant(c1_go_loci,  fdr_threshold), "\n")
cat("  W2 GO  :", count_significant(w2_go_loci,  fdr_threshold), "\n")
cat("  C1 KEGG:", count_significant(c1_kegg_loci, fdr_threshold), "\n")
cat("  W2 KEGG:", count_significant(w2_kegg_loci, fdr_threshold), "\n")

cat("\nW2 GO terms surviving collapse:\n")
print(w2_go_loci %>% filter(adjusted_p_value < fdr_threshold) %>%
        select(description, ontology, selected_gene_count, fold_enrichment, adjusted_p_value))
cat("\nW2 KEGG terms surviving collapse:\n")
print(w2_kegg_loci %>% filter(adjusted_p_value < fdr_threshold) %>%
        select(description, selected_gene_count, fold_enrichment, adjusted_p_value))


### PLOT ###



plot_enrichment_bar <- function(enrichment_df, plot_title,
                                fdr_threshold = 0.05,
                                label_col = "selected_gene_count") {
  
  df <- enrichment_df %>%
    filter(!is.na(adjusted_p_value), adjusted_p_value < fdr_threshold) %>%
    mutate(description = fct_reorder(description, fold_enrichment))  # largest at top
  
  if (nrow(df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0,
                 label = paste0("No terms with FDR < ", fdr_threshold)) +
        theme_void() + ggtitle(plot_title)
    )
  }
  
  ggplot(df, aes(x = description, y = fold_enrichment, fill = adjusted_p_value)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = .data[[label_col]]),
              hjust = -0.25, size = 3.2, colour = "grey20") +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +     # room for labels
    scale_fill_gradient(low = "#b2182b", high = "#2166ac",
                        name = "adj. p", labels = scales::label_scientific()) +
    labs(title = plot_title, x = NULL, y = "Fold enrichment") +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# build the two plots (Fig 2B = W2 GO over W2 KEGG)
p_w2_go   <- plot_enrichment_bar(w2_go_enrichment,   "W2 GO enrichment (gene-level) - Molecular Function (MF)")
p_w2_kegg <- plot_enrichment_bar(w2_kegg_enrichment, "W2 KEGG enrichment (gene-level)")

# shared fill limits so both panels use one identical adj-p colour legend
all_padj <- c(
  w2_go_enrichment$adjusted_p_value[w2_go_enrichment$adjusted_p_value < 0.05],
  w2_kegg_enrichment$adjusted_p_value[w2_kegg_enrichment$adjusted_p_value < 0.05]
)
fill_limits <- range(all_padj, na.rm = TRUE)

shared_fill <- scale_fill_gradient(
  low = "#b2182b", high = "#2166ac",
  name = "adj. p", labels = scales::scientific,
  limits = fill_limits
)

p_w2_go   <- p_w2_go   + shared_fill
p_w2_kegg <- p_w2_kegg + shared_fill

# stack with a single collected legend (GO panel 3x the KEGG panel height)
p_stacked <- p_w2_go / p_w2_kegg +
  plot_layout(guides = "collect", heights = c(3, 1)) +
  plot_annotation(tag_levels = "A")

print(p_stacked)

# NOTE [repo path]: replace absolute output path with e.g. "figures/W2_enrichment_GO_KEGG.png"
ggsave("C:/Users/nirossi/Documents/Paper parmeniaca/figures_final/W2_enrichment_GO_KEGG.png", p_stacked, width = 8, height = 6, dpi = 300)

# =============================================================================
# NOMINAL VS ANNOTATED FOREGROUND / BACKGROUND GENE COUNTS
# =============================================================================

count_annotated_gene_sets <- function(annotation, database_name) {
  
  annotated_genes <- unique(annotation$gene_id)
  
  data.frame(
    database = database_name,
    set = c(
      "W2 foreground genes",
      "C1 foreground genes",
      "All foreground genes",
      "Background genes"
    ),
    nominal_n_genes = c(
      length(unique(w2_genes)),
      length(unique(c1_genes)),
      length(unique(c(w2_genes, c1_genes))),
      length(unique(background_genes))
    ),
    annotated_n_genes = c(
      length(intersect(unique(w2_genes), annotated_genes)),
      length(intersect(unique(c1_genes), annotated_genes)),
      length(intersect(unique(c(w2_genes, c1_genes)), annotated_genes)),
      length(intersect(unique(background_genes), annotated_genes))
    )
  )
}

foreground_background_annotation_sizes <- bind_rows(
  count_annotated_gene_sets(go_annotations, "GO"),
  count_annotated_gene_sets(kegg_pathway_annotations, "KEGG pathways")
)

print(foreground_background_annotation_sizes)




# =============================================================================
# CREATE DF FOR PLOTTING THE POSITION OF GO-ENRICHED VARIANTS
#
# Builds go_variants / kegg_variants: the genomic positions (CHR, POS) of the
# deleterious variants in W2-enriched-term genes. These tables are the input to
# the circos plot (05.R, Fig 2C). NB: del_pos is still the PROVEAN/HIGH
# del_variants set, so positions cover all candidate deleterious variants in the
# selected genes (not only the PROVEAN < -9 foreground used for gene selection).
# =============================================================================


library(tidyverse)

# =============================================================================
# 1. Get the genes belonging to each significant W2 term
# =============================================================================
# significant GO terms
sig_go_terms <- w2_go_enrichment %>%
  filter(adjusted_p_value < 0.05) %>%
  pull(term_id)

# significant KEGG terms
sig_kegg_terms <- w2_kegg_enrichment %>%
  filter(adjusted_p_value < 0.05) %>%
  pull(term_id)

# genes annotated to those terms (from the cleaned annotation tables)
go_term_genes <- go_annotations %>%
  filter(term_id %in% sig_go_terms) %>%
  select(term_id, description, gene_id)

kegg_term_genes <- kegg_pathway_annotations %>%
  filter(term_id %in% sig_kegg_terms) %>%
  select(term_id, description, gene_id)

# =============================================================================
# 2. Restrict to genes that are actually in the W2-enriched selected set
#    (the term contains many genes; you want the ones driving your signal)
# =============================================================================
go_hits   <- go_term_genes   %>% filter(gene_id %in% w2_genes)
kegg_hits <- kegg_pathway_annotations %>%
  filter(term_id %in% sig_kegg_terms, gene_id %in% w2_genes) %>%
  select(term_id, description, gene_id)

# =============================================================================
# 3. Attach chromosome + position of the deleterious variants in those genes
# =============================================================================
del_pos <- del_variants %>%
  mutate(gene_id = normalize_gene_id(geneID)) %>%
  select(gene_id, CHR, POS, IMPACT, EFFECT, PROVEAN_SCORE)

go_variants <- go_hits %>%
  left_join(del_pos, by = "gene_id") %>%
  arrange(term_id, CHR, POS)

kegg_variants <- kegg_hits %>%
  left_join(del_pos, by = "gene_id") %>%
  arrange(term_id, CHR, POS)

cat("GO term variants:\n"); print(head(go_variants, 20))
cat("\nKEGG term variants:\n"); print(head(kegg_variants, 20))

cat("\nUnique genes (GO):", n_distinct(go_variants$gene_id), "\n")
cat("Unique variants (GO):", nrow(distinct(go_variants, CHR, POS)), "\n")




go_variants %>% distinct(CHR, POS, gene_id) %>% nrow()


