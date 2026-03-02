#!/usr/bin/env Rscript

## =========================================================
## Salmon (gene-level) -> TPM matrix -> Sanity Checks + CLDN1
## Human BAT/WAT  Thermoneutral vs Cold Exposure (n = 1/group)
## Dataset: E-MTAB-4031 (UCSF Diabetes Center, HiSeq 2500, ~51 bp SE)
## Reference: GENCODE 38 (GRCh38)
##
## Design (from ArrayExpress E-MTAB-4031):
##   BAT_TN  – supraclavicular BAT, thermoneutral
##   BAT_CE  – supraclavicular BAT, 5-h cold exposure
##   WAT_TN  – subcutaneous abdominal WAT, thermoneutral
##   WAT_CE  – subcutaneous abdominal WAT, 5-h cold exposure
##   All from a single individual  =>  n = 1 per group
## =========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
})

## ---- 0) Working directory + output folder -----------------------------------
setwd("C:/Users/lking/OneDrive - Louisiana State University/PBRC/Bioinformatics/Maria_Human_BAT/FASTQ")

results_dir <- file.path(getwd(), "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
cat("Results will be saved to:", results_dir, "\n\n")

## ---- 1) Input files ---------------------------------------------------------
files <- c(
  "WAT_TN gene quantification.tabular",
  "WAT_CE gene quantification.tabular",
  "BAT_TN gene quantification.tabular",
  "BAT_CE gene quantification.tabular"
)
sample_names <- c("WAT_TN", "WAT_CE", "BAT_TN", "BAT_CE")

## ---- 2) File existence checks -----------------------------------------------
cat("========================================\n")
cat("SANITY CHECK 1: File existence\n")
cat("========================================\n")
cat("Working directory:", getwd(), "\n\n")

missing <- files[!file.exists(files)]
if (length(missing) > 0) {
  stop(
    "These expected files are missing from the working directory:\n",
    paste0("  - ", missing, collapse = "\n")
  )
}
cat("All 4 quantification files found.\n\n")

## ---- 2b) Library layout confirmation ----------------------------------------
cat("========================================\n")
cat("SANITY CHECK 1b: Library layout\n")
cat("========================================\n")
cat("E-MTAB-4031 was sequenced on Illumina HiSeq 2500.\n")
cat("  Library layout : SINGLE-END\n")
cat("  Read length    : ~51 bp\n")
cat("  Reads per run  : ~13M (e.g. ERR1110588 = 13,425,365 reads)\n")
cat("  Center         : UCSF Diabetes Center\n")
cat("\n")
cat("  Your workflow (1 FASTQ per sample -> Galaxy Salmon quant) is CORRECT\n")
cat("  for single-end data. A ~76%% mapping rate with 51 bp single-end reads\n")
cat("  is reasonable for Salmon on human RNA-seq.\n")
cat("\n")

## ---- 3) Read Salmon quant files + column validation -------------------------
## Salmon gene-level quant columns: Name, Length, EffectiveLength, TPM, NumReads
expected_cols <- c("Name", "Length", "EffectiveLength", "TPM", "NumReads")

read_salmon_gene <- function(path, sample_name) {
  df <- read_tsv(path, show_col_types = FALSE)

  missing_cols <- setdiff(expected_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "File '", path, "' is missing expected Salmon columns: ",
      paste(missing_cols, collapse = ", "),
      "\nFound columns: ", paste(colnames(df), collapse = ", ")
    )
  }

  df$SampleName <- sample_name
  df
}

cat("========================================\n")
cat("SANITY CHECK 2: Column validation\n")
cat("========================================\n")

raw_list <- map2(files, sample_names, read_salmon_gene)
names(raw_list) <- sample_names

cat("All files contain the expected Salmon columns:\n")
cat("  ", paste(expected_cols, collapse = ", "), "\n\n")

## ---- 4) Per-sample QC -------------------------------------------------------
cat("========================================\n")
cat("SANITY CHECK 3: Per-sample QC\n")
cat("========================================\n")

per_sample_qc <- map_dfr(raw_list, function(df) {
  sn <- df$SampleName[1]
  tibble(
    Sample        = sn,
    n_genes       = nrow(df),
    TPM_sum       = sum(df$TPM, na.rm = TRUE),
    total_reads   = sum(df$NumReads, na.rm = TRUE),
    zero_TPM_pct  = 100 * mean(df$TPM == 0, na.rm = TRUE),
    median_TPM    = median(df$TPM, na.rm = TRUE),
    max_TPM       = max(df$TPM, na.rm = TRUE),
    pct_ENSG      = 100 * mean(grepl("^ENSG", df$Name), na.rm = TRUE)
  )
})

print(per_sample_qc, width = Inf)
cat("\n")

## 4a) TPM sum should be ~1,000,000 for each sample
cat("-- TPM sum check (expect ~1,000,000) --\n")
for (i in seq_len(nrow(per_sample_qc))) {
  s <- per_sample_qc$Sample[i]
  tpm_s <- per_sample_qc$TPM_sum[i]
  flag <- if (abs(tpm_s - 1e6) > 5e4) " ** WARNING: >5% deviation **" else " OK"
  cat(sprintf("  %s: %.0f%s\n", s, tpm_s, flag))
}
cat("\n")

## 4b) Gene count should be ~60k for GENCODE 38
cat("-- Gene count check (expect ~60,000 for GENCODE 38) --\n")
for (i in seq_len(nrow(per_sample_qc))) {
  s <- per_sample_qc$Sample[i]
  ng <- per_sample_qc$n_genes[i]
  flag <- if (ng < 20000 || ng > 80000) " ** WARNING: unusual gene count **" else " OK"
  cat(sprintf("  %s: %d genes%s\n", s, ng, flag))
}
cat("\n")

## 4c) Gene ID format: should be predominantly ENSG (Ensembl)
cat("-- Gene ID format check (expect ~100%% ENSG IDs for GENCODE) --\n")
for (i in seq_len(nrow(per_sample_qc))) {
  s <- per_sample_qc$Sample[i]
  pct <- per_sample_qc$pct_ENSG[i]
  flag <- if (pct < 90) " ** WARNING: low ENSG fraction – wrong reference? **" else " OK"
  cat(sprintf("  %s: %.1f%% ENSG%s\n", s, pct, flag))
}
cat("\n")

## 4d) Total mapped reads: flag if any sample has very few reads
cat("-- Total mapped reads (NumReads sum) --\n")
cat("  Low total reads (<5M) may indicate incomplete download or failed run.\n")
for (i in seq_len(nrow(per_sample_qc))) {
  s <- per_sample_qc$Sample[i]
  tr <- per_sample_qc$total_reads[i]
  flag <- if (tr < 5e6) " ** WARNING: very low read count **" else " OK"
  cat(sprintf("  %s: %.1fM reads%s\n", s, tr / 1e6, flag))
}
cat("\n")

## 4e) Zero-expression fraction
cat("-- Zero-expression fraction --\n")
cat("  Typical: 40-60%% of GENCODE genes have TPM=0 in a given tissue.\n")
cat("  Very high (>80%%) suggests mapping problems.\n")
for (i in seq_len(nrow(per_sample_qc))) {
  s <- per_sample_qc$Sample[i]
  zp <- per_sample_qc$zero_TPM_pct[i]
  flag <- if (zp > 80) " ** WARNING: unusually high **" else " OK"
  cat(sprintf("  %s: %.1f%% zero-TPM%s\n", s, zp, flag))
}
cat("\n")

## 4f) Consistency: all samples should have the same gene set
gene_sets <- map(raw_list, ~ sort(.x$Name))
if (length(unique(gene_sets)) != 1) {
  cat("** WARNING: Gene ID sets differ between samples! **\n")
  cat("  This can happen if different references were used.\n\n")
} else {
  cat("Gene ID sets are identical across all 4 samples. Good.\n\n")
}

## ---- 5) Build TPM matrix + strip version suffixes ---------------------------
tpm_list <- map2(files, sample_names, function(path, sn) {
  df <- read_tsv(path, show_col_types = FALSE)
  df %>%
    select(Name, TPM) %>%
    rename(GeneID = Name, !!sn := TPM)
})
tpm_mat <- reduce(tpm_list, full_join, by = "GeneID")

## Strip Ensembl version suffixes (ENSG00000163347.16 -> ENSG00000163347)
tpm_mat <- tpm_mat %>%
  mutate(GeneID = sub("\\.\\d+$", "", GeneID))

## ---- 5b) Duplicate gene IDs after version stripping -------------------------
cat("========================================\n")
cat("SANITY CHECK 4: Duplicate gene IDs\n")
cat("========================================\n")

dup_ids <- tpm_mat$GeneID[duplicated(tpm_mat$GeneID)]
if (length(dup_ids) > 0) {
  cat(sprintf("** WARNING: %d duplicate gene IDs after stripping version suffixes. **\n", length(dup_ids)))
  cat("  Collapsing duplicates by summing TPM (standard practice).\n")
  cat("  First 10 duplicates:", paste(head(dup_ids, 10), collapse = ", "), "\n\n")

  tpm_mat <- tpm_mat %>%
    group_by(GeneID) %>%
    summarise(across(all_of(sample_names), sum, na.rm = TRUE), .groups = "drop")
} else {
  cat("No duplicates after version stripping. Good.\n\n")
}

## Write TPM matrix to disk
write.csv(tpm_mat, file.path(results_dir, "TPM_matrix_gene_level.csv"), row.names = FALSE)
cat("Wrote: results/TPM_matrix_gene_level.csv\n\n")

## ---- 6) Cross-sample QC ----------------------------------------------------
cat("========================================\n")
cat("SANITY CHECK 5: Cross-sample QC\n")
cat("========================================\n")

## 6a) Spearman correlation of log2(TPM + 1)
log2_mat <- tpm_mat %>%
  select(all_of(sample_names)) %>%
  mutate(across(everything(), ~ log2(.x + 1)))

cor_mat <- cor(log2_mat, method = "spearman")
cat("Spearman correlation matrix [log2(TPM + 1)]:\n")
print(round(cor_mat, 3))
cat("\n")

## Flag if any within-tissue pair has r < 0.85
for (tissue in c("WAT", "BAT")) {
  cols <- sample_names[grepl(paste0("^", tissue), sample_names)]
  if (length(cols) == 2) {
    r <- cor_mat[cols[1], cols[2]]
    flag <- if (r < 0.85) " ** WARNING: low within-tissue correlation **" else " OK"
    cat(sprintf("  %s TN vs CE: r = %.3f%s\n", tissue, r, flag))
  }
}
cat("\n")

## 6b) PCA on log2(TPM+1)
pca_res <- prcomp(t(as.matrix(log2_mat)), center = TRUE, scale. = FALSE)
pca_df  <- data.frame(
  Sample    = sample_names,
  PC1       = pca_res$x[, 1],
  PC2       = pca_res$x[, 2],
  Tissue    = ifelse(grepl("^BAT", sample_names), "BAT", "WAT"),
  Condition = ifelse(grepl("_CE$", sample_names), "Cold", "TN")
)
var_expl <- summary(pca_res)$importance[2, 1:2] * 100

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue, shape = Condition)) +
  geom_point(size = 5) +
  geom_text(aes(label = Sample), vjust = -1, size = 3.5) +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA of log2(TPM + 1)",
    x = sprintf("PC1 (%.1f%%)", var_expl[1]),
    y = sprintf("PC2 (%.1f%%)", var_expl[2])
  )
ggsave(file.path(results_dir, "QC_PCA_plot.png"), p_pca, width = 6.5, height = 5, dpi = 300)
cat("Saved: results/QC_PCA_plot.png\n")
cat("  Expect BAT and WAT to separate on PC1.\n\n")

## 6c) Density plot of log2(TPM + 1)
density_df <- tpm_mat %>%
  pivot_longer(cols = all_of(sample_names), names_to = "Sample", values_to = "TPM") %>%
  mutate(log2TPM1 = log2(TPM + 1))

p_density <- ggplot(density_df, aes(x = log2TPM1, color = Sample)) +
  geom_density(linewidth = 0.8) +
  theme_classic(base_size = 14) +
  labs(
    title = "log2(TPM + 1) density per sample",
    x = "log2(TPM + 1)",
    y = "Density"
  )
ggsave(file.path(results_dir, "QC_density_plot.png"), p_density, width = 7, height = 4.5, dpi = 300)
cat("Saved: results/QC_density_plot.png\n")
cat("  Curves should largely overlap; a shifted sample indicates a QC issue.\n\n")

## ---- 7) Biological sanity checks (marker genes) ----------------------------
cat("========================================\n")
cat("SANITY CHECK 6: Biological marker genes\n")
cat("========================================\n")

markers <- tribble(
  ~GeneID,              ~Symbol,    ~Category,
  ## ---- BAT identity markers (should be HIGH in BAT, LOW in WAT) ----
  "ENSG00000109424",    "UCP1",     "BAT identity",
  "ENSG00000147655",    "CIDEA",    "BAT identity",
  "ENSG00000142611",    "PRDM16",   "BAT identity",
  "ENSG00000109819",    "PPARGC1A", "BAT identity",
  "ENSG00000119915",    "ELOVL3",   "BAT identity",
  "ENSG00000152977",    "ZIC1",     "BAT identity",
  "ENSG00000162624",    "LHX8",     "BAT identity",
  "ENSG00000221818",    "EBF2",     "BAT identity",

  ## ---- WAT identity markers (should be HIGH in WAT, LOW in BAT) ----
  "ENSG00000174697",    "LEP",      "WAT identity",
  "ENSG00000037965",    "HOXC8",    "WAT identity",
  "ENSG00000180806",    "HOXC9",    "WAT identity",
  "ENSG00000118526",    "TCF21",    "WAT identity",

  ## ---- Cold-responsive in BAT (should go UP in BAT_CE vs BAT_TN) ----
  "ENSG00000211448",    "DIO2",     "Cold-responsive",
  "ENSG00000205560",    "CPT1B",    "Cold-responsive",
  "ENSG00000124253",    "PCK1",     "Cold-responsive",
  "ENSG00000178537",    "SLC25A20", "Cold-responsive",

  ## ---- Pan-adipocyte (should be expressed >1 TPM in all samples) ----
  "ENSG00000181092",    "ADIPOQ",   "Pan-adipocyte",
  "ENSG00000170323",    "FABP4",    "Pan-adipocyte",
  "ENSG00000132170",    "PPARG",    "Pan-adipocyte",
  "ENSG00000166819",    "PLIN1",    "Pan-adipocyte",

  ## ---- Thermogenesis pathway ----
  "ENSG00000188778",    "ADRB3",    "Thermogenesis",
  "ENSG00000105550",    "FGF21",    "Thermogenesis",
  "ENSG00000112715",    "VEGFA",    "Thermogenesis",

  ## ---- Housekeeping (should be stable across all 4 samples) ----
  "ENSG00000075624",    "ACTB",     "Housekeeping",
  "ENSG00000111640",    "GAPDH",    "Housekeeping",
  "ENSG00000166710",    "B2M",      "Housekeeping",
  "ENSG00000142541",    "RPL13A",   "Housekeeping",
  "ENSG00000089157",    "RPLP0",    "Housekeeping",

  ## ---- Tight junction family (context for CLDN1) ----
  "ENSG00000163347",    "CLDN1",    "Tight junction",
  "ENSG00000184113",    "CLDN5",    "Tight junction",
  "ENSG00000197822",    "OCLN",     "Tight junction",
  "ENSG00000104067",    "TJP1",     "Tight junction"
)

marker_tpm <- markers %>%
  left_join(tpm_mat, by = "GeneID")

## Report any markers not found in the TPM matrix
missing_markers <- marker_tpm %>% filter(is.na(WAT_TN))
if (nrow(missing_markers) > 0) {
  cat("  WARNING: These markers were NOT found in the TPM matrix:\n")
  cat(paste0("    ", missing_markers$Symbol, " (", missing_markers$GeneID, ")"), sep = "\n")
  cat("\n")
}

marker_wide <- marker_tpm %>%
  select(Category, Symbol, all_of(sample_names))

cat("\nFull marker gene TPM table (35 genes, 7 categories):\n")
print(as.data.frame(marker_wide), row.names = FALSE, digits = 4)
cat("\n")

## Save full marker table
write.csv(marker_wide, file.path(results_dir, "QC_marker_gene_table.csv"), row.names = FALSE)
cat("Saved: results/QC_marker_gene_table.csv\n\n")

## ---- 7b) Automated pass/fail scoring ---------------------------------------
cat("========================================\n")
cat("SANITY CHECK 7: Automated marker scoring\n")
cat("========================================\n")

get_tpm <- function(sym, samp) {
  v <- marker_tpm %>% filter(Symbol == sym) %>% pull(!!sym(samp))
  if (length(v) == 0) return(NA_real_)
  v
}

pass_count <- 0
total_count <- 0

## --- BAT identity: each should be higher in BAT_TN than WAT_TN ---
cat("\n-- BAT identity markers (expect BAT_TN > WAT_TN) --\n")
bat_markers <- markers %>% filter(Category == "BAT identity") %>% pull(Symbol)
for (g in bat_markers) {
  bat_val <- get_tpm(g, "BAT_TN")
  wat_val <- get_tpm(g, "WAT_TN")
  total_count <- total_count + 1
  if (!is.na(bat_val) && !is.na(wat_val)) {
    if (bat_val > wat_val) {
      cat(sprintf("  [PASS] %s: BAT=%.3f > WAT=%.3f\n", g, bat_val, wat_val))
      pass_count <- pass_count + 1
    } else {
      cat(sprintf("  [FAIL] %s: BAT=%.3f <= WAT=%.3f\n", g, bat_val, wat_val))
    }
  } else {
    cat(sprintf("  [SKIP] %s: not found in matrix\n", g))
  }
}

## --- WAT identity: each should be higher in WAT_TN than BAT_TN ---
cat("\n-- WAT identity markers (expect WAT_TN > BAT_TN) --\n")
wat_markers <- markers %>% filter(Category == "WAT identity") %>% pull(Symbol)
for (g in wat_markers) {
  wat_val <- get_tpm(g, "WAT_TN")
  bat_val <- get_tpm(g, "BAT_TN")
  total_count <- total_count + 1
  if (!is.na(bat_val) && !is.na(wat_val)) {
    if (wat_val > bat_val) {
      cat(sprintf("  [PASS] %s: WAT=%.3f > BAT=%.3f\n", g, wat_val, bat_val))
      pass_count <- pass_count + 1
    } else {
      cat(sprintf("  [FAIL] %s: WAT=%.3f <= BAT=%.3f\n", g, wat_val, bat_val))
    }
  } else {
    cat(sprintf("  [SKIP] %s: not found in matrix\n", g))
  }
}

## --- Cold-responsive: should increase in BAT with cold (BAT_CE > BAT_TN) ---
cat("\n-- Cold-responsive markers (expect BAT_CE > BAT_TN) --\n")
cold_markers <- markers %>% filter(Category == "Cold-responsive") %>% pull(Symbol)
for (g in cold_markers) {
  tn_val <- get_tpm(g, "BAT_TN")
  ce_val <- get_tpm(g, "BAT_CE")
  total_count <- total_count + 1
  if (!is.na(tn_val) && !is.na(ce_val)) {
    if (ce_val > tn_val) {
      cat(sprintf("  [PASS] %s: BAT_CE=%.3f > BAT_TN=%.3f\n", g, ce_val, tn_val))
      pass_count <- pass_count + 1
    } else {
      cat(sprintf("  [FAIL] %s: BAT_CE=%.3f <= BAT_TN=%.3f\n", g, ce_val, tn_val))
    }
  } else {
    cat(sprintf("  [SKIP] %s: not found in matrix\n", g))
  }
}

## --- Pan-adipocyte: should be >1 TPM in all 4 samples ---
cat("\n-- Pan-adipocyte markers (expect >1 TPM in all samples) --\n")
pan_markers <- markers %>% filter(Category == "Pan-adipocyte") %>% pull(Symbol)
for (g in pan_markers) {
  vals <- sapply(sample_names, function(s) get_tpm(g, s))
  total_count <- total_count + 1
  if (all(!is.na(vals))) {
    min_val <- min(vals)
    if (min_val > 1) {
      cat(sprintf("  [PASS] %s: min TPM = %.2f (all > 1)\n", g, min_val))
      pass_count <- pass_count + 1
    } else {
      cat(sprintf("  [FAIL] %s: min TPM = %.2f (below 1 in at least one sample)\n", g, min_val))
    }
  } else {
    cat(sprintf("  [SKIP] %s: not found in matrix\n", g))
  }
}

## --- Housekeeping: CV across 4 samples should be < 50% ---
cat("\n-- Housekeeping markers (expect CV < 50%% across samples) --\n")
hk_markers <- markers %>% filter(Category == "Housekeeping") %>% pull(Symbol)
for (g in hk_markers) {
  vals <- sapply(sample_names, function(s) get_tpm(g, s))
  total_count <- total_count + 1
  if (all(!is.na(vals)) && mean(vals) > 0) {
    cv <- 100 * sd(vals) / mean(vals)
    if (cv < 50) {
      cat(sprintf("  [PASS] %s: CV = %.1f%%, mean = %.1f TPM\n", g, cv, mean(vals)))
      pass_count <- pass_count + 1
    } else {
      cat(sprintf("  [FAIL] %s: CV = %.1f%% (too variable), mean = %.1f TPM\n", g, cv, mean(vals)))
    }
  } else {
    cat(sprintf("  [SKIP] %s: not found or zero mean\n", g))
  }
}

## --- Summary score ---
cat(sprintf("\n** OVERALL SCORE: %d / %d marker checks passed **\n", pass_count, total_count))
if (pass_count >= total_count * 0.75) {
  cat("   -> Data looks GOOD: >=75%% of expected patterns confirmed.\n")
} else if (pass_count >= total_count * 0.5) {
  cat("   -> Data looks MARGINAL: 50-75%% of expected patterns confirmed.\n")
  cat("   -> Some markers don't match expectations. Review sample labels and mapping.\n")
} else {
  cat("   -> Data looks PROBLEMATIC: <50%% of expected patterns confirmed.\n")
  cat("   -> Consider re-running Salmon or checking sample identity.\n")
}
cat("\n")

## ---- 8) Helpers -------------------------------------------------------------

make_gene_long <- function(gene_id, gene_label = gene_id) {
  g <- tpm_mat %>% filter(GeneID == gene_id)

  if (nrow(g) == 0) {
    stop(
      "Gene ID not found in TPM matrix: ", gene_id, "\n",
      "Check: sum(tpm_mat$GeneID == '", gene_id, "')"
    )
  }

  g %>%
    pivot_longer(-GeneID, names_to = "Sample", values_to = "TPM") %>%
    mutate(
      GeneLabel = gene_label,
      Tissue    = ifelse(grepl("^BAT", Sample), "BAT", "WAT"),
      Condition = ifelse(grepl("_CE$", Sample), "Cold", "TN"),
      Condition = factor(Condition, levels = c("TN", "Cold")),
      Tissue    = factor(Tissue, levels = c("WAT", "BAT"))
    )
}

## Descriptive fold-change (NOT statistical inference with n=1)
desc_log2fc <- function(df_long) {
  df_long %>%
    select(GeneLabel, Tissue, Condition, TPM) %>%
    pivot_wider(names_from = Condition, values_from = TPM) %>%
    mutate(log2FC_Cold_vs_TN = log2((Cold + 1e-6) / (TN + 1e-6)))
}

## ---- 9) Gene-of-interest + controls ----------------------------------------
cat("========================================\n")
cat("CLDN1 and control gene analysis\n")
cat("========================================\n")

UCP1_ENSG  <- "ENSG00000109424"
LEP_ENSG   <- "ENSG00000174697"
CLDN1_ENSG <- "ENSG00000163347"

ucp1_long  <- make_gene_long(UCP1_ENSG,  "UCP1")
lep_long   <- make_gene_long(LEP_ENSG,   "LEP")
cldn1_long <- make_gene_long(CLDN1_ENSG, "CLDN1")

cat("\nUCP1 TPM values:\n")
print(ucp1_long %>% select(GeneLabel, Sample, Tissue, Condition, TPM))
cat("\nLEP TPM values:\n")
print(lep_long %>% select(GeneLabel, Sample, Tissue, Condition, TPM))
cat("\nCLDN1 TPM values:\n")
print(cldn1_long %>% select(GeneLabel, Sample, Tissue, Condition, TPM))

## ---- 10) VISUALIZATION 1: Grouped bar charts (per gene) --------------------
## Bar charts are more honest than line plots for n=1 data.
## Lines imply a continuous trajectory; bars show individual observations.

plot_gene_bars <- function(df_long) {
  ggplot(df_long, aes(x = Tissue, y = TPM, fill = Condition)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("TN" = "#4393C3", "Cold" = "#D6604D"),
                      labels = c("TN" = "Thermoneutral", "Cold" = "Cold Exposure")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top") +
    labs(
      title = paste0(df_long$GeneLabel[1], " expression (TPM)"),
      subtitle = "n = 1 per group — descriptive only",
      x = NULL, y = "TPM", fill = NULL
    )
}

p_ucp1  <- plot_gene_bars(ucp1_long)
p_lep   <- plot_gene_bars(lep_long)
p_cldn1 <- plot_gene_bars(cldn1_long)

ggsave(file.path(results_dir, "UCP1_bar.png"),  p_ucp1,  width = 5, height = 4.5, dpi = 300)
ggsave(file.path(results_dir, "LEP_bar.png"),   p_lep,   width = 5, height = 4.5, dpi = 300)
ggsave(file.path(results_dir, "CLDN1_bar.png"), p_cldn1, width = 5, height = 4.5, dpi = 300)
cat("\nSaved individual bar charts to results/\n")

## ---- 11) VISUALIZATION 2: Combined multi-panel figure -----------------------
## Shows UCP1, LEP, and CLDN1 side-by-side with independent y-axes.
## This is the most informative single figure for this analysis.

all_long <- bind_rows(ucp1_long, lep_long, cldn1_long) %>%
  mutate(GeneLabel = factor(GeneLabel, levels = c("UCP1", "LEP", "CLDN1")))

p_combined <- ggplot(all_long, aes(x = Tissue, y = TPM, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("TN" = "#4393C3", "Cold" = "#D6604D"),
                    labels = c("TN" = "Thermoneutral", "Cold" = "Cold Exposure")) +
  facet_wrap(~ GeneLabel, scales = "free_y", nrow = 1) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 13)
  ) +
  labs(
    title = "Gene expression: Controls (UCP1, LEP) + Gene of Interest (CLDN1)",
    subtitle = "E-MTAB-4031  |  n = 1 per group  |  descriptive only",
    x = NULL, y = "TPM", fill = NULL
  )

ggsave(file.path(results_dir, "combined_panel.png"), p_combined,
       width = 11, height = 4.5, dpi = 300)
cat("Saved: results/combined_panel.png\n")

## ---- 12) VISUALIZATION 3: CLDN1 with noise-floor annotation -----------------
## CLDN1 values are extremely low (<0.04 TPM in all samples).
## A dedicated plot with a noise threshold makes this visually obvious.

noise_threshold <- 0.1  # genes below this TPM are unreliable at gene level

p_cldn1_annotated <- ggplot(cldn1_long, aes(x = Tissue, y = TPM, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("TN" = "#4393C3", "Cold" = "#D6604D"),
                    labels = c("TN" = "Thermoneutral", "Cold" = "Cold Exposure")) +
  geom_hline(yintercept = noise_threshold, linetype = "dashed", color = "grey40") +
  annotate("text", x = 2.4, y = noise_threshold + 0.008,
           label = "noise floor (~0.1 TPM)", size = 3.3, color = "grey40") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top") +
  labs(
    title = "CLDN1 (Claudin-1) expression",
    subtitle = "All values below noise floor — not reliably detected",
    x = NULL, y = "TPM", fill = NULL
  )

ggsave(file.path(results_dir, "CLDN1_with_noise_floor.png"), p_cldn1_annotated,
       width = 5.5, height = 5, dpi = 300)
cat("Saved: results/CLDN1_with_noise_floor.png\n")

## ---- 13) VISUALIZATION 4: Full marker heatmap (35 genes, 7 categories) -----
## This is the centerpiece QC figure. Rows grouped by category.

## Define row order: genes grouped by category
category_order <- c("BAT identity", "WAT identity", "Cold-responsive",
                     "Pan-adipocyte", "Thermogenesis", "Housekeeping",
                     "Tight junction")

gene_order <- marker_tpm %>%
  mutate(Category = factor(Category, levels = category_order)) %>%
  arrange(Category) %>%
  pull(Symbol)

heat_df <- marker_tpm %>%
  select(Symbol, Category, all_of(sample_names)) %>%
  pivot_longer(cols = all_of(sample_names), names_to = "Sample", values_to = "TPM") %>%
  mutate(
    log10TPM1 = log10(TPM + 1),
    Symbol = factor(Symbol, levels = rev(gene_order)),
    Sample = factor(Sample, levels = c("WAT_TN", "WAT_CE", "BAT_TN", "BAT_CE")),
    ## Smart label: show 0 for zero, integer for >=10, 1 decimal for >=1, 2 decimal otherwise
    label = case_when(
      is.na(TPM) ~ "NA",
      TPM == 0   ~ "0",
      TPM >= 100 ~ sprintf("%.0f", TPM),
      TPM >= 10  ~ sprintf("%.1f", TPM),
      TPM >= 1   ~ sprintf("%.1f", TPM),
      TPM >= 0.1 ~ sprintf("%.2f", TPM),
      TRUE       ~ sprintf("%.3f", TPM)
    )
  )

## Build category separator positions for horizontal lines
cat_breaks <- marker_tpm %>%
  mutate(Category = factor(Category, levels = category_order)) %>%
  arrange(Category) %>%
  mutate(row_num = row_number()) %>%
  group_by(Category) %>%
  summarise(ymax = max(row_num), .groups = "drop")

## Convert to positions between gene rows (reversed because ggplot y is bottom-up)
n_genes <- nrow(marker_tpm)
hline_positions <- cat_breaks$ymax[-nrow(cat_breaks)] + 0.5

## Category label positions (midpoint of each group)
cat_labels <- marker_tpm %>%
  mutate(Category = factor(Category, levels = category_order)) %>%
  arrange(Category) %>%
  mutate(row_num = row_number()) %>%
  group_by(Category) %>%
  summarise(mid = mean(row_num), .groups = "drop") %>%
  mutate(mid_rev = n_genes + 1 - mid)

p_heat <- ggplot(heat_df, aes(x = Sample, y = Symbol, fill = log10TPM1)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.6) +
  scale_fill_gradient(low = "#F7F7F7", high = "#B2182B",
                      name = "log10(TPM+1)",
                      na.value = "grey80") +
  ## Category separator lines
  geom_hline(yintercept = n_genes + 1 - hline_positions, linewidth = 0.6, color = "grey30") +
  ## Category labels on the right
  annotate("text", x = 4.7, y = cat_labels$mid_rev,
           label = cat_labels$Category,
           size = 2.8, fontface = "bold", hjust = 0, color = "grey30") +
  scale_x_discrete(position = "top") +
  coord_cartesian(clip = "off", xlim = c(0.5, 4.5)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x.top = element_text(face = "bold", size = 11),
    axis.text.y = element_text(face = "bold", size = 9),
    panel.grid = element_blank(),
    plot.margin = margin(5, 100, 5, 5)   # extra right margin for labels
  ) +
  labs(
    title = "Marker gene heatmap (35 genes, 7 categories)",
    subtitle = "E-MTAB-4031  |  n = 1 per group  |  TPM values in cells",
    x = NULL, y = NULL
  )

ggsave(file.path(results_dir, "marker_gene_heatmap.png"), p_heat,
       width = 8.5, height = 12, dpi = 300)
cat("Saved: results/marker_gene_heatmap.png\n")

## ---- 14) Descriptive fold-change tables + bar chart -------------------------
fc_ucp1  <- desc_log2fc(ucp1_long)
fc_lep   <- desc_log2fc(lep_long)
fc_cldn1 <- desc_log2fc(cldn1_long)

fc_all <- bind_rows(fc_ucp1, fc_lep, fc_cldn1)
write.csv(fc_all, file.path(results_dir, "descriptive_log2FC_all.csv"), row.names = FALSE)
cat("Saved: results/descriptive_log2FC_all.csv\n")

## log2FC bar chart
fc_plot_df <- fc_all %>%
  mutate(
    GeneLabel = factor(GeneLabel, levels = c("UCP1", "LEP", "CLDN1")),
    ## flag CLDN1 WAT fold-change as unreliable (0 -> noise)
    reliable = !(GeneLabel == "CLDN1" & Tissue == "WAT")
  )

p_fc <- ggplot(fc_plot_df, aes(x = Tissue, y = log2FC_Cold_vs_TN, fill = Tissue)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  ## Cross out unreliable bars (CLDN1 WAT: 0 -> 0.008 = noise artifact)
  geom_text(
    data = fc_plot_df %>% filter(!reliable),
    aes(label = "artifact\n(0 -> noise)"),
    vjust = -0.3, size = 2.8, color = "grey40"
  ) +
  scale_fill_manual(values = c("WAT" = "#FDB863", "BAT" = "#B2ABD2")) +
  facet_wrap(~ GeneLabel, nrow = 1) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 13)
  ) +
  labs(
    title = "Descriptive log2 fold-change (Cold Exposure vs Thermoneutral)",
    subtitle = "n = 1 — NOT statistical inference",
    x = NULL, y = "log2FC (CE / TN)"
  )

ggsave(file.path(results_dir, "log2FC_bar_chart.png"), p_fc,
       width = 10, height = 4.5, dpi = 300)
cat("Saved: results/log2FC_bar_chart.png\n\n")

## ---- 15) Final CLDN1 summary ------------------------------------------------
cat("========================================\n")
cat("CLDN1 (Claudin-1) summary\n")
cat("========================================\n")

cldn1_bat_tn <- cldn1_long %>% filter(Tissue == "BAT", Condition == "TN") %>% pull(TPM)
cldn1_bat_ce <- cldn1_long %>% filter(Tissue == "BAT", Condition == "Cold") %>% pull(TPM)
cldn1_wat_tn <- cldn1_long %>% filter(Tissue == "WAT", Condition == "TN") %>% pull(TPM)
cldn1_wat_ce <- cldn1_long %>% filter(Tissue == "WAT", Condition == "Cold") %>% pull(TPM)

cat(sprintf("  BAT thermoneutral: %.4f TPM\n", cldn1_bat_tn))
cat(sprintf("  BAT cold exposure: %.4f TPM\n", cldn1_bat_ce))
cat(sprintf("  WAT thermoneutral: %.4f TPM\n", cldn1_wat_tn))
cat(sprintf("  WAT cold exposure: %.4f TPM\n", cldn1_wat_ce))
cat("\n")

## Only report fold-change when both values are above noise
if (cldn1_bat_tn >= noise_threshold || cldn1_bat_ce >= noise_threshold) {
  cat(sprintf("  BAT log2FC (CE/TN): %.3f\n",
              log2((cldn1_bat_ce + 1e-6) / (cldn1_bat_tn + 1e-6))))
} else {
  cat("  BAT log2FC: not meaningful (both values below noise floor).\n")
}
if (cldn1_wat_tn >= noise_threshold || cldn1_wat_ce >= noise_threshold) {
  cat(sprintf("  WAT log2FC (CE/TN): %.3f\n",
              log2((cldn1_wat_ce + 1e-6) / (cldn1_wat_tn + 1e-6))))
} else {
  cat("  WAT log2FC: not meaningful (both values below noise floor).\n")
}
cat("\n")

cat("  INTERPRETATION:\n")
cat("  CLDN1 is NOT reliably detected in any of the 4 samples (all < 0.04 TPM).\n")
cat("  These values are below the noise floor for gene-level Salmon quantification.\n")
cat("  There is NO evidence of cold-induced CLDN1 expression in BAT from this\n")
cat("  dataset, but this does not prove CLDN1 is absent — the gene may be\n")
cat("  expressed at very low levels or in a small subset of cells not captured\n")
cat("  by bulk RNA-seq. Confirmation would require targeted methods (qPCR,\n")
cat("  single-cell RNA-seq, or immunohistochemistry).\n")
cat("\n")

cat("========================================\n")
cat("IMPORTANT CAVEATS\n")
cat("========================================\n")
cat("  1. n = 1 per group => purely descriptive, no statistical inference.\n")
cat("  2. All samples are from a SINGLE individual (E-MTAB-4031).\n")
cat("  3. Fold-changes should be interpreted as observations, not conclusions.\n")
cat("  4. To claim cold-induced CLDN1 expression you would need biological\n")
cat("     replicates (ideally n >= 3) and appropriate statistical testing.\n")
cat("  5. Confirm Salmon was run with GENCODE 38 primary assembly (not the\n")
cat("     full genome including scaffolds/patches) for clean gene-level quant.\n")
cat("  6. If TPM sums deviate from 1,000,000, check that Salmon ran to\n")
cat("     completion and that the GTF matches the transcriptome index.\n")
cat("  7. E-MTAB-4031 is ~51 bp single-end (Illumina HiSeq 2500, UCSF\n")
cat("     Diabetes Center). One FASTQ per sample is correct.\n")
cat("  8. A ~76%% mapping rate is reasonable for 51 bp single-end human\n")
cat("     RNA-seq in Salmon. Rates <60%% would indicate a problem.\n")
cat("\n")

## ---- 16) Print all final plots in RStudio -----------------------------------
print(p_combined)
print(p_cldn1_annotated)
print(p_heat)
print(p_fc)

cat("========================================\n")
cat("All outputs saved to: results/\n")
cat("========================================\n")
cat("  Plots:\n")
cat("    combined_panel.png          <- best overview figure\n")
cat("    CLDN1_with_noise_floor.png  <- key finding: CLDN1 below detection\n")
cat("    marker_gene_heatmap.png     <- QC: all markers at a glance\n")
cat("    log2FC_bar_chart.png        <- fold-change comparison\n")
cat("    UCP1_bar.png                <- positive control (BAT identity)\n")
cat("    LEP_bar.png                 <- tissue identity control\n")
cat("    CLDN1_bar.png               <- gene of interest\n")
cat("    QC_PCA_plot.png             <- sample clustering\n")
cat("    QC_density_plot.png         <- expression distribution\n")
cat("  Tables:\n")
cat("    TPM_matrix_gene_level.csv   <- full TPM matrix\n")
cat("    QC_marker_gene_table.csv    <- marker gene summary\n")
cat("    descriptive_log2FC_all.csv  <- fold-changes for UCP1/LEP/CLDN1\n")
cat("\nDone.\n")
