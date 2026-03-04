#!/usr/bin/env Rscript

## =========================================================
## Human BAT/WAT (E-MTAB-4031) — Focused Salmon Analysis
##
## Key changes from previous version:
##   - Focused on CLDN1 specifically (not all claudins)
##   - Added "reality check" figure: CLDN1 alongside reference
##     genes so you can see the scale difference
##   - Kept FA metabolism heatmap for paper replication
##   - Added raw-expression heatmap option (no z-score) so you
##     can see absolute levels, not just relative patterns
##   - Removed the broad claudin family scan
##
## The original paper (Jespersen et al.) used:
##   TopHat2 + Cufflinks/Cuffdiff, paired-end, hg19, older GENCODE
## You used:
##   Salmon, single-end, latest GENCODE
## Expect similar trends for well-expressed genes, but low-
## expression genes (like CLDN1) may differ at the noise floor.
##
## Usage:
##   Rscript scripts/salmon_heatmap_analysis.R
## =========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
})

## ---- Configuration ---------------------------------------------------------
## Point this to wherever your 4 Salmon gene quant .tabular files live
DATA_DIR <- "FASTQ"

RESULTS_DIR <- "results"
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

sample_info <- tibble(
  file      = c("WAT_TN gene quantification.tabular",
                "WAT_CE gene quantification.tabular",
                "BAT_TN gene quantification.tabular",
                "BAT_CE gene quantification.tabular"),
  sample    = c("WAT_TN", "WAT_CE", "BAT_TN", "BAT_CE"),
  tissue    = c("WAT", "WAT", "BAT", "BAT"),
  condition = c("TN", "CE", "TN", "CE")
)

sample_order <- c("WAT_TN", "WAT_CE", "BAT_TN", "BAT_CE")

cat("============================================================\n")
cat("Human BAT/WAT Transcriptional Analysis (E-MTAB-4031)\n")
cat("============================================================\n\n")

## ---- 1. Read Salmon quant files --------------------------------------------
cat("--- Reading Salmon gene quantification files ---\n")

read_salmon <- function(filepath, sample_name) {
  if (!file.exists(filepath))
    stop("'", filepath, "' does not exist in current working directory ('",
         getwd(), "').")
  df <- read_tsv(filepath, show_col_types = FALSE)
  required <- c("Name", "Length", "EffectiveLength", "TPM", "NumReads")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0)
    stop("File '", filepath, "' missing columns: ", paste(missing, collapse = ", "))
  df %>%
    mutate(
      GeneID = sub("\\.\\d+$", "", Name),
      sample = sample_name
    )
}

quant_list <- map2(
  file.path(DATA_DIR, sample_info$file),
  sample_info$sample,
  read_salmon
)
names(quant_list) <- sample_info$sample

cat("All 4 files read successfully.\n\n")

## ---- 2. Per-sample QC ------------------------------------------------------
cat("--- Per-sample QC ---\n")

lib_sizes <- map_dfr(quant_list, function(df) {
  tibble(
    sample          = df$sample[1],
    total_reads     = sum(df$NumReads, na.rm = TRUE),
    genes_detected  = sum(df$TPM > 0, na.rm = TRUE),
    genes_above_1TPM = sum(df$TPM > 1, na.rm = TRUE),
    tpm_sum         = round(sum(df$TPM, na.rm = TRUE), 1)
  )
})
print(lib_sizes, width = Inf)

cat("\n")
cat("NOTE: BAT_CE has the most reads but fewest genes >1 TPM.\n")
cat("This means cold-exposed BAT concentrates reads into a small\n")
cat("set of thermogenic genes (UCP1 etc). This is expected biology,\n")
cat("but it distorts z-score heatmaps — most genes will appear\n")
cat("'downregulated' in BAT_CE relative to other samples.\n\n")

## ---- 3. Build TPM + FPKM matrices -----------------------------------------
tpm_wide <- map(quant_list, ~ .x %>% select(GeneID, TPM)) %>%
  imap(~ rename(.x, !!.y := TPM)) %>%
  reduce(full_join, by = "GeneID") %>%
  mutate(across(all_of(sample_order), ~ replace_na(.x, 0)))

fpkm_wide <- {
  parts <- map(sample_info$sample, function(sn) {
    df <- quant_list[[sn]]
    libsize <- sum(df$NumReads, na.rm = TRUE)
    df %>%
      mutate(
        FPKM = (NumReads * 1e9) / (pmax(EffectiveLength, 1) * pmax(libsize, 1))
      ) %>%
      select(GeneID, FPKM) %>%
      rename(!!sn := FPKM)
  })
  reduce(parts, full_join, by = "GeneID") %>%
    mutate(across(all_of(sample_order), ~ replace_na(.x, 0)))
}

write_csv(tpm_wide,  file.path(RESULTS_DIR, "TPM_matrix.csv"))
write_csv(fpkm_wide, file.path(RESULTS_DIR, "FPKM_matrix.csv"))
cat("Saved: TPM_matrix.csv, FPKM_matrix.csv\n\n")

## ---- 4. QC: Spearman correlation + PCA -------------------------------------
log2tpm <- tpm_wide %>%
  select(all_of(sample_order)) %>%
  mutate(across(everything(), ~ log2(.x + 1)))

cor_mat <- cor(log2tpm, method = "spearman")
cat("--- Spearman correlation [log2(TPM+1)] ---\n")
print(round(cor_mat, 3))
cat("\n")

pca <- prcomp(t(as.matrix(log2tpm)), center = TRUE, scale. = FALSE)
pca_df <- tibble(
  sample    = sample_order,
  PC1       = pca$x[, 1],
  PC2       = pca$x[, 2],
  tissue    = sample_info$tissue,
  condition = sample_info$condition
)
ve <- summary(pca)$importance[2, 1:2] * 100

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = tissue, shape = condition)) +
  geom_point(size = 5) +
  geom_text(aes(label = sample), vjust = -1.2, size = 3.5) +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA — log2(TPM+1)",
    x = sprintf("PC1 (%.1f%%)", ve[1]),
    y = sprintf("PC2 (%.1f%%)", ve[2])
  )
ggsave(file.path(RESULTS_DIR, "PCA_plot.png"), p_pca,
       width = 7, height = 5, dpi = 300)
cat("Saved: PCA_plot.png\n\n")

## ---- 5. Helpers -------------------------------------------------------------
zscore_row <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

## ---- 6. Paper-style heatmap builder -----------------------------------------
## Blue-white-MAGENTA, z-scored log2(FPKM+1), clipped [-1,1]

make_paper_heatmap <- function(gene_table, fpkm_matrix, title, filename,
                               width = 6, height = NULL,
                               z_limits = c(-1, 1),
                               preserve_order = FALSE) {

  panel <- gene_table %>%
    inner_join(fpkm_matrix, by = "GeneID") %>%
    filter(if_any(all_of(sample_order), ~ .x > 0))

  if (nrow(panel) == 0) {
    cat("  WARNING: No genes found for '", title, "'\n")
    return(invisible(NULL))
  }

  zmat <- panel %>%
    select(all_of(sample_order)) %>%
    mutate(across(everything(), ~ log2(.x + 1))) %>%
    as.matrix()
  zmat <- t(apply(zmat, 1, zscore_row))
  colnames(zmat) <- sample_order

  panel_long <- tibble(
    Symbol   = rep(panel$Symbol, each = length(sample_order)),
    Category = rep(panel$Category, each = length(sample_order)),
    Sample   = rep(sample_order, nrow(panel)),
    z        = as.vector(t(zmat))
  ) %>%
    mutate(z_clip = pmax(z_limits[1], pmin(z_limits[2], z)))

  if (preserve_order) {
    gene_levels <- panel$Symbol
  } else {
    cat_levels <- unique(gene_table$Category)
    gene_levels <- panel %>%
      mutate(Category = factor(Category, levels = cat_levels)) %>%
      arrange(Category, Symbol) %>%
      pull(Symbol)
  }

  panel_long <- panel_long %>%
    mutate(
      Symbol = factor(Symbol, levels = rev(gene_levels)),
      Sample = factor(Sample, levels = sample_order)
    )

  if (is.null(height)) height <- max(4, length(gene_levels) * 0.38 + 2)

  p <- ggplot(panel_long, aes(x = Sample, y = Symbol, fill = z_clip)) +
    geom_tile(color = "black", linewidth = 0.4) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#C51B7D",
      midpoint = 0, limits = z_limits,
      oob = scales::squish,
      name = "z-score"
    ) +
    scale_x_discrete(position = "top") +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid    = element_blank(),
      axis.title    = element_blank(),
      axis.text.x.top = element_text(angle = 45, hjust = 0,
                                     face = "bold", size = 11),
      axis.text.y   = element_text(face = "italic", size = 10),
      legend.position = "right"
    ) +
    labs(title = title,
         subtitle = sprintf("z-scored log2(FPKM+1), clipped [%g, %g]",
                            z_limits[1], z_limits[2]))

  ggsave(file.path(RESULTS_DIR, filename), p,
         width = width, height = height, dpi = 300)
  cat("  Saved:", filename, "\n")
  invisible(p)
}

## ---- 7. Gene panels ---------------------------------------------------------

## 7a. Paper FA / lipid metabolism panel — EXACT order from the figure
fa_genes_paper <- tribble(
  ~Symbol,      ~GeneID,            ~Category,
  "HMGCS2",     "ENSG00000134240",  "FA / lipid metabolism",
  "ACADVL",     "ENSG00000072778",  "FA / lipid metabolism",
  "ECHS1",      "ENSG00000127884",  "FA / lipid metabolism",
  "DGAT1",      "ENSG00000185000",  "FA / lipid metabolism",
  "DGAT2",      "ENSG00000062282",  "FA / lipid metabolism",
  "SLC25A20",   "ENSG00000178537",  "FA / lipid metabolism",
  "HADHB",      "ENSG00000138029",  "FA / lipid metabolism",
  "ECI1",       "ENSG00000167969",  "FA / lipid metabolism",
  "CPT1B",      "ENSG00000205560",  "FA / lipid metabolism",
  "AGPAT3",     "ENSG00000160216",  "FA / lipid metabolism",
  "PPARA",      "ENSG00000186951",  "FA / lipid metabolism",
  "ACLY",       "ENSG00000131473",  "FA / lipid metabolism",
  "DECR1",      "ENSG00000104325",  "FA / lipid metabolism",
  "GK",         "ENSG00000198814",  "FA / lipid metabolism",
  "ACADM",      "ENSG00000117054",  "FA / lipid metabolism",
  "BDH1",       "ENSG00000161267",  "FA / lipid metabolism",
  "DLD",        "ENSG00000091140",  "FA / lipid metabolism",
  "ACAA2",      "ENSG00000167315",  "FA / lipid metabolism"
)

fa_genes_extra <- tribble(
  ~Symbol,  ~GeneID,            ~Category,
  "CD36",   "ENSG00000135218",  "FA uptake (paper text)",
  "LPL",    "ENSG00000175445",  "FA uptake (paper text)"
)

## 7b. Core identity markers
identity_genes <- tribble(
  ~Symbol,      ~GeneID,            ~Category,
  "UCP1",       "ENSG00000109424",  "BAT / thermogenesis",
  "CIDEA",      "ENSG00000147655",  "BAT / thermogenesis",
  "PRDM16",     "ENSG00000142611",  "BAT / thermogenesis",
  "PPARGC1A",   "ENSG00000109819",  "BAT / thermogenesis",
  "ELOVL3",     "ENSG00000119915",  "BAT / thermogenesis",
  "DIO2",       "ENSG00000211448",  "BAT / thermogenesis",
  "LEP",        "ENSG00000174697",  "WAT marker",
  "HOXC8",      "ENSG00000037965",  "WAT marker",
  "HOXC9",      "ENSG00000180806",  "WAT marker",
  "ADIPOQ",     "ENSG00000181092",  "Pan-adipocyte",
  "FABP4",      "ENSG00000170323",  "Pan-adipocyte",
  "PPARG",      "ENSG00000132170",  "Pan-adipocyte"
)

## ---- 8. Generate heatmaps --------------------------------------------------
cat("\n--- Generating heatmaps ---\n")

## Paper-exact FA metabolism
p_fa <- make_paper_heatmap(
  fa_genes_paper, fpkm_wide,
  "FA / Lipid Metabolism (paper figure replication)",
  "heatmap_FA_metabolism.png", width = 6,
  preserve_order = TRUE
)

## Extended with CD36 + LPL
p_fa_extra <- make_paper_heatmap(
  bind_rows(fa_genes_paper, fa_genes_extra), fpkm_wide,
  "FA / Lipid Metabolism + CD36/LPL",
  "heatmap_FA_metabolism_extended.png", width = 6,
  preserve_order = TRUE
)

## Identity markers
p_id <- make_paper_heatmap(
  identity_genes, fpkm_wide,
  "BAT/WAT Identity Markers",
  "heatmap_identity_markers.png", width = 6
)

## ---- 9. CLDN1 focused analysis ---------------------------------------------
cat("\n--- CLDN1 focused analysis ---\n")

## 9a. Raw evidence directly from each Salmon file
cldn1_evidence <- map2_dfr(
  file.path(DATA_DIR, sample_info$file),
  sample_info$sample,
  function(path, sn) {
    read_tsv(path, show_col_types = FALSE) %>%
      mutate(Name = sub("\\.\\d+$", "", Name)) %>%
      filter(Name == "ENSG00000163347") %>%
      transmute(
        sample          = sn,
        GeneID          = Name,
        TPM             = as.numeric(TPM),
        NumReads        = as.numeric(NumReads),
        Length          = as.numeric(Length),
        EffectiveLength = as.numeric(EffectiveLength)
      )
  }
)

cat("\nCLDN1 raw Salmon quantification:\n")
print(cldn1_evidence, width = Inf)
cat("\n")

max_cldn1_reads <- max(cldn1_evidence$NumReads, na.rm = TRUE)
max_cldn1_tpm   <- max(cldn1_evidence$TPM, na.rm = TRUE)

if (max_cldn1_reads <= 10) {
  cat("*** CLDN1 has 0-", max_cldn1_reads, " reads across all samples. ***\n",
      "*** This is below reliable quantification (~0.03 TPM max). ***\n",
      "*** Differences between samples are stochastic noise,     ***\n",
      "*** not biology. CLDN1 is not expressed in this dataset.   ***\n\n",
      sep = "")
}

write_csv(cldn1_evidence, file.path(RESULTS_DIR, "CLDN1_evidence.csv"))

## 9b. CLDN1 in context: bar plot alongside well-expressed reference genes
##     This is the KEY figure — shows CLDN1's scale relative to genes
##     you know are real. Much more informative than a heatmap.
context_genes <- tribble(
  ~Symbol, ~GeneID,            ~Category,
  "UCP1",  "ENSG00000109424",  "BAT marker (positive ctrl)",
  "LEP",   "ENSG00000174697",  "WAT marker (positive ctrl)",
  "ADIPOQ","ENSG00000181092",  "Pan-adipocyte (positive ctrl)",
  "FABP4", "ENSG00000170323",  "Pan-adipocyte (positive ctrl)",
  "CLDN1", "ENSG00000163347",  "Gene of interest"
)

context_df <- context_genes %>%
  left_join(tpm_wide, by = "GeneID") %>%
  pivot_longer(cols = all_of(sample_order),
               names_to = "sample", values_to = "TPM") %>%
  mutate(
    tissue = ifelse(grepl("^BAT", sample), "BAT", "WAT"),
    sample = factor(sample, levels = sample_order),
    Symbol = factor(Symbol, levels = context_genes$Symbol),
    is_cldn1 = Symbol == "CLDN1"
  )

p_context <- ggplot(context_df, aes(x = sample, y = TPM, fill = tissue)) +
  geom_col(width = 0.7) +
  facet_wrap(~ Symbol, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("BAT" = "#D6604D", "WAT" = "#4393C3")) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
    strip.text    = element_text(face = "bold.italic", size = 10),
    legend.position = "top"
  ) +
  labs(
    title    = "CLDN1 in context: TPM alongside reference genes",
    subtitle = "Note the y-axis scale — CLDN1 is orders of magnitude below detectable genes",
    x = NULL, y = "TPM", fill = NULL
  )

ggsave(file.path(RESULTS_DIR, "CLDN1_in_context.png"), p_context,
       width = 12, height = 5, dpi = 300)
cat("Saved: CLDN1_in_context.png\n")

## 9c. CLDN1-only close-up with read counts annotated
##     Zoomed in so you can see the tiny values, but with read counts
##     on the bars so it's clear how little evidence there is.
cldn1_bar_df <- cldn1_evidence %>%
  mutate(
    tissue = ifelse(grepl("^BAT", sample), "BAT", "WAT"),
    sample = factor(sample, levels = sample_order),
    read_label = paste0(NumReads, " read", ifelse(NumReads != 1, "s", ""))
  )

p_cldn1_zoom <- ggplot(cldn1_bar_df, aes(x = sample, y = TPM, fill = tissue)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = read_label), vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("BAT" = "#D6604D", "WAT" = "#4393C3")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top") +
  labs(
    title    = "CLDN1 (Claudin-1) expression — zoomed in",
    subtitle = "0–3 reads per sample; below reliable detection threshold",
    x = NULL, y = "TPM", fill = NULL
  )

ggsave(file.path(RESULTS_DIR, "CLDN1_zoomed.png"), p_cldn1_zoom,
       width = 6, height = 5, dpi = 300)
cat("Saved: CLDN1_zoomed.png\n")

## ---- 10. Positive controls bar plot ----------------------------------------
ctrl_genes <- tribble(
  ~Symbol, ~GeneID,
  "UCP1",  "ENSG00000109424",
  "LEP",   "ENSG00000174697"
)

ctrl_df <- ctrl_genes %>%
  left_join(tpm_wide, by = "GeneID") %>%
  pivot_longer(cols = all_of(sample_order),
               names_to = "sample", values_to = "TPM") %>%
  mutate(
    tissue = ifelse(grepl("^BAT", sample), "BAT", "WAT"),
    sample = factor(sample, levels = sample_order),
    Symbol = factor(Symbol, levels = c("UCP1", "LEP"))
  )

p_ctrl <- ggplot(ctrl_df, aes(x = sample, y = TPM, fill = tissue)) +
  geom_col(width = 0.7) +
  facet_wrap(~ Symbol, scales = "free_y") +
  scale_fill_manual(values = c("BAT" = "#D6604D", "WAT" = "#4393C3")) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text  = element_text(face = "bold.italic", size = 13),
    legend.position = "top"
  ) +
  labs(
    title    = "Positive controls — tissue identity",
    subtitle = "UCP1 high in BAT-CE = data is working; LEP high in WAT = correct tissue ID",
    x = NULL, y = "TPM", fill = NULL
  )

ggsave(file.path(RESULTS_DIR, "positive_controls.png"), p_ctrl,
       width = 8, height = 5, dpi = 300)
cat("Saved: positive_controls.png\n")

## ---- 11. Descriptive fold-change table for FA panel -------------------------
fa_fc <- bind_rows(fa_genes_paper, fa_genes_extra) %>%
  inner_join(fpkm_wide, by = "GeneID") %>%
  filter(if_any(all_of(sample_order), ~ .x > 0)) %>%
  mutate(
    BAT_log2FC_CE_vs_TN = log2((BAT_CE + 0.01) / (BAT_TN + 0.01)),
    WAT_log2FC_CE_vs_TN = log2((WAT_CE + 0.01) / (WAT_TN + 0.01))
  ) %>%
  select(Symbol, Category, all_of(sample_order),
         BAT_log2FC_CE_vs_TN, WAT_log2FC_CE_vs_TN)

write_csv(fa_fc, file.path(RESULTS_DIR, "FA_panel_FPKM_and_log2FC.csv"))
cat("Saved: FA_panel_FPKM_and_log2FC.csv\n")

## ---- 12. Summary ------------------------------------------------------------
cat("\n============================================================\n")
cat("DONE — key outputs in", RESULTS_DIR, "/\n")
cat("============================================================\n")
cat("\n")
cat("  Matrices:\n")
cat("    TPM_matrix.csv\n")
cat("    FPKM_matrix.csv\n")
cat("\n")
cat("  Heatmaps (blue-white-magenta, z-scored log2(FPKM+1)):\n")
cat("    heatmap_FA_metabolism.png          <- paper figure replication\n")
cat("    heatmap_FA_metabolism_extended.png <- + CD36/LPL from paper text\n")
cat("    heatmap_identity_markers.png       <- BAT/WAT identity markers\n")
cat("\n")
cat("  CLDN1 analysis (YOUR KEY FIGURES):\n")
cat("    CLDN1_in_context.png   <- CLDN1 alongside reference genes (scale!)\n")
cat("    CLDN1_zoomed.png       <- close-up with read counts on bars\n")
cat("    CLDN1_evidence.csv     <- raw Salmon evidence for CLDN1\n")
cat("\n")
cat("  Other:\n")
cat("    PCA_plot.png\n")
cat("    positive_controls.png\n")
cat("    FA_panel_FPKM_and_log2FC.csv\n")
cat("============================================================\n")
cat("\n")
cat("INTERPRETATION GUIDE:\n")
cat("  1. Check positive_controls.png first — if UCP1 is high in\n")
cat("     BAT_CE and LEP is high in WAT, your data is valid.\n")
cat("  2. CLDN1 has 0-3 reads across all samples. This is stochastic\n")
cat("     noise, not a biological signal. CLDN1 is not meaningfully\n")
cat("     expressed in bulk adipose tissue in this dataset.\n")
cat("  3. The FA heatmaps may not perfectly match the paper because\n")
cat("     you used a different pipeline (Salmon SE + latest GENCODE\n")
cat("     vs TopHat2 PE + older GENCODE). Trends should be similar\n")
cat("     for well-expressed genes; exact fold-changes will differ.\n")
cat("  4. With n=1 per condition (no biological replicates), these\n")
cat("     are descriptive patterns only — no statistical significance\n")
cat("     can be claimed.\n")
cat("============================================================\n")
