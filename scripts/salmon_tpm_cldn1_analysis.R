#!/usr/bin/env Rscript

## =========================================================
## Salmon (gene-level) -> TPM matrix -> Sanity Checks + CLDN1
## Human BAT/WAT  Thermoneutral vs Cold Exposure (n = 1/group)
## Dataset: E-MTAB-4031 (Jespersen / Scheele, Rigshospitalet)
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

## ---- 0) Working directory ---------------------------------------------------
setwd("C:/Users/lking/OneDrive - Louisiana State University/PBRC/Bioinformatics/Maria_Human_BAT/FASTQ")

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

## ---- 2b) Paired-end / library-type reminder ---------------------------------
cat("========================================\n")
cat("SANITY CHECK 1b: Paired-end handling\n")
cat("========================================\n")
cat("E-MTAB-4031 was sequenced with 100 bp paired-end reads (Illumina).\n")
cat("You should have provided BOTH R1 and R2 FASTQ files to Salmon.\n")
cat("\n")
cat("  In Galaxy Salmon quant:\n")
cat("    - Select 'Paired-end' under library type\n")
cat("    - Provide both mate 1 (R1) and mate 2 (R2) files\n")
cat("    - If you only uploaded 1 FASTQ per sample, Salmon treated it as\n")
cat("      single-end, which wastes half your data and lowers mapping rates.\n")
cat("\n")
cat("  How to check: Go back to your Galaxy history. For each sample you should\n")
cat("  have downloaded TWO .fastq.gz files from ENA (e.g. *_1.fastq.gz and\n")
cat("  *_2.fastq.gz). If you only see one file per sample, re-download from:\n")
cat("    https://www.ebi.ac.uk/ena/browser/view/E-MTAB-4031\n")
cat("\n")
cat("  Expected Salmon mapping rate for paired-end human RNA-seq: 75-90%%.\n")
cat("  If you ran paired-end correctly and got ~76%%, that is acceptable.\n")
cat("  If you ran single-end only (R1), ~76%% is expected but suboptimal.\n")
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
write.csv(tpm_mat, "TPM_matrix_gene_level.csv", row.names = FALSE)
cat("Wrote: TPM_matrix_gene_level.csv\n\n")

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
ggsave("QC_PCA_plot.png", p_pca, width = 6.5, height = 5, dpi = 300)
cat("Saved: QC_PCA_plot.png\n")
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
ggsave("QC_density_plot.png", p_density, width = 7, height = 4.5, dpi = 300)
cat("Saved: QC_density_plot.png\n")
cat("  Curves should largely overlap; a shifted sample indicates a QC issue.\n\n")

## ---- 7) Biological sanity checks (marker genes) ----------------------------
cat("========================================\n")
cat("SANITY CHECK 6: Biological marker genes\n")
cat("========================================\n")

markers <- tribble(
  ~GeneID,              ~Symbol,  ~Expectation,
  "ENSG00000109424",    "UCP1",   "High in BAT, low/absent in WAT; may increase with cold in BAT",
  "ENSG00000174697",    "LEP",    "High in WAT, low in BAT",
  "ENSG00000163347",    "CLDN1",  "Gene of interest (claudin-1)",
  "ENSG00000147655",    "CIDEA",  "BAT marker; higher in BAT than WAT",
  "ENSG00000211448",    "DIO2",   "Cold-responsive in BAT (type 2 deiodinase)",
  "ENSG00000181092",    "ADIPOQ", "Adiponectin; expressed in both, often higher in WAT",
  "ENSG00000006025",    "OSBPL7", "Negative control – housekeeping-like, should be stable"
)

marker_tpm <- markers %>%
  left_join(tpm_mat, by = "GeneID")

cat("\nMarker gene TPM values:\n")
marker_long <- marker_tpm %>%
  pivot_longer(
    cols = all_of(sample_names),
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  select(Symbol, Sample, TPM, Expectation)

marker_wide <- marker_tpm %>%
  select(Symbol, all_of(sample_names), Expectation)

print(as.data.frame(marker_wide), row.names = FALSE)
cat("\n")

## Biological plausibility flags
get_tpm <- function(sym, samp) {
  v <- marker_tpm %>% filter(Symbol == sym) %>% pull(!!sym(samp))
  if (length(v) == 0) return(NA_real_)
  v
}

ucp1_bat_tn <- get_tpm("UCP1", "BAT_TN")
ucp1_wat_tn <- get_tpm("UCP1", "WAT_TN")
lep_wat_tn  <- get_tpm("LEP",  "WAT_TN")
lep_bat_tn  <- get_tpm("LEP",  "BAT_TN")

cat("Biological plausibility checks:\n")

if (!is.na(ucp1_bat_tn) && !is.na(ucp1_wat_tn)) {
  if (ucp1_bat_tn > ucp1_wat_tn) {
    cat("  [PASS] UCP1 is higher in BAT_TN than WAT_TN.\n")
  } else {
    cat("  [WARN] UCP1 is NOT higher in BAT than WAT — samples may be swapped!\n")
  }
} else {
  cat("  [WARN] UCP1 not found in TPM matrix — check gene ID ENSG00000109424.\n")
}

if (!is.na(lep_wat_tn) && !is.na(lep_bat_tn)) {
  if (lep_wat_tn > lep_bat_tn) {
    cat("  [PASS] LEP is higher in WAT_TN than BAT_TN.\n")
  } else {
    cat("  [WARN] LEP is NOT higher in WAT — unexpected, check sample labels.\n")
  }
} else {
  cat("  [WARN] LEP not found in TPM matrix — check gene ID ENSG00000174697.\n")
}

## CIDEA check
cidea_bat <- get_tpm("CIDEA", "BAT_TN")
cidea_wat <- get_tpm("CIDEA", "WAT_TN")
if (!is.na(cidea_bat) && !is.na(cidea_wat)) {
  if (cidea_bat > cidea_wat) {
    cat("  [PASS] CIDEA is higher in BAT_TN than WAT_TN.\n")
  } else {
    cat("  [WARN] CIDEA is NOT higher in BAT — review sample identity.\n")
  }
}

## DIO2 cold-responsiveness in BAT
dio2_bat_tn <- get_tpm("DIO2", "BAT_TN")
dio2_bat_ce <- get_tpm("DIO2", "BAT_CE")
if (!is.na(dio2_bat_tn) && !is.na(dio2_bat_ce)) {
  cat(sprintf("  [INFO] DIO2 in BAT: TN = %.2f, CE = %.2f (cold-responsive marker).\n",
              dio2_bat_tn, dio2_bat_ce))
}
cat("\n")

## Save marker table
write.csv(marker_wide, "QC_marker_gene_table.csv", row.names = FALSE)
cat("Saved: QC_marker_gene_table.csv\n\n")

## ---- 8) Helpers for gene plots and fold-change ------------------------------

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

plot_gene_lines <- function(df_long) {
  ggplot(df_long, aes(x = Condition, y = TPM, group = Tissue, color = Tissue)) +
    geom_point(size = 4) +
    geom_line(linewidth = 1.2) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(df_long$GeneLabel[1], " expression (TPM)  [n = 1/group]"),
      subtitle = "Descriptive only — no statistical inference",
      x = NULL,
      y = "TPM"
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

## ---- 10) Plots + fold-change tables ----------------------------------------
p_ucp1  <- plot_gene_lines(ucp1_long)
p_lep   <- plot_gene_lines(lep_long)
p_cldn1 <- plot_gene_lines(cldn1_long)

ggsave("UCP1_TPM_control_plot.png",   p_ucp1,  width = 5.5, height = 4.5, dpi = 300)
ggsave("LEP_TPM_WAT_marker_plot.png", p_lep,   width = 5.5, height = 4.5, dpi = 300)
ggsave("CLDN1_TPM_plot.png",          p_cldn1, width = 5.5, height = 4.5, dpi = 300)

cat("\nSaved gene plots:\n")
cat("  - UCP1_TPM_control_plot.png\n")
cat("  - LEP_TPM_WAT_marker_plot.png\n")
cat("  - CLDN1_TPM_plot.png\n\n")

## Descriptive fold-changes
fc_ucp1  <- desc_log2fc(ucp1_long)
fc_lep   <- desc_log2fc(lep_long)
fc_cldn1 <- desc_log2fc(cldn1_long)

write.csv(fc_ucp1,  "UCP1_descriptive_log2FC.csv",  row.names = FALSE)
write.csv(fc_lep,   "LEP_descriptive_log2FC.csv",   row.names = FALSE)
write.csv(fc_cldn1, "CLDN1_descriptive_log2FC.csv", row.names = FALSE)

cat("Wrote descriptive fold-change tables:\n")
cat("  - UCP1_descriptive_log2FC.csv\n")
cat("  - LEP_descriptive_log2FC.csv\n")
cat("  - CLDN1_descriptive_log2FC.csv\n\n")

## ---- 11) Final CLDN1 summary ------------------------------------------------
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
cat(sprintf("  BAT log2FC (CE/TN): %.3f\n", log2((cldn1_bat_ce + 1e-6) / (cldn1_bat_tn + 1e-6))))
cat(sprintf("  WAT log2FC (CE/TN): %.3f\n", log2((cldn1_wat_ce + 1e-6) / (cldn1_wat_tn + 1e-6))))
cat("\n")

if (cldn1_bat_ce > cldn1_bat_tn) {
  cat("  -> CLDN1 TPM is HIGHER in BAT after cold exposure.\n")
} else if (cldn1_bat_ce < cldn1_bat_tn) {
  cat("  -> CLDN1 TPM is LOWER in BAT after cold exposure.\n")
} else {
  cat("  -> CLDN1 TPM is unchanged in BAT after cold exposure.\n")
}

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
cat("  7. E-MTAB-4031 is 100 bp paired-end (Illumina). Verify you gave\n")
cat("     Salmon BOTH R1 and R2 FASTQs per sample, not just one file.\n")
cat("  8. The original study (Sun et al. 2018, Nat Commun) used TopHat/hg19.\n")
cat("     Re-quantifying with Salmon/GENCODE38 is valid but expect minor\n")
cat("     differences from published values.\n")
cat("  9. A mapping rate of ~76%% is acceptable for Salmon with human data.\n")
cat("     Rates <60%% would indicate a problem (wrong reference, truncated\n")
cat("     files, or adapter contamination). Rates >85%% are ideal.\n")
cat("\n")

print(p_ucp1)
print(p_lep)
print(p_cldn1)

cat("Done.\n")
