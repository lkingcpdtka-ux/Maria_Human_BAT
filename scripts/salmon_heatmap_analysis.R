#!/usr/bin/env Rscript

## =========================================================
## Human BAT/WAT (E-MTAB-4031) — Improved Salmon Analysis
##
## Fixes vs original script:
##   1. FPKM uses EffectiveLength (not Length) — matches how
##      Cuffdiff computes FPKM internally
##   2. Blue-white-magenta palette matching MeV's actual rendering
##      (paper caption says "red" but the figure is pink/magenta)
##   3. Added CD36 + LPL to FA panel (paper highlights these)
##   4. Expanded claudin family (CLDN1–15 + TJ scaffolds)
##   5. Trimmed marker panel to essentials
##   6. Cleaner structure, fewer redundant outputs
##
## Usage:
##   Rscript scripts/salmon_heatmap_analysis.R
##
## Before running:
##   - Set DATA_DIR below to your Salmon gene quant folder
##   - Install: readr, dplyr, purrr, tidyr, ggplot2
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
DATA_DIR <- "data/salmon_quant"

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

## ---- 3. Build TPM + corrected FPKM matrices --------------------------------
## TPM: directly from Salmon (already normalized)
## FPKM: (NumReads * 1e9) / (EffectiveLength * LibrarySize)
##        ^^^^^^^^^^^^^^^ KEY FIX: use EffectiveLength, not Length

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
        ## EffectiveLength accounts for fragment-length + positional bias
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

## ---- 5. Helper: z-score per row --------------------------------------------
zscore_row <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

## ---- 6. Paper-style heatmap builder ----------------------------------------
## Blue-white-MAGENTA, z-scored log2(FPKM+1), clipped [-1,1]
## Paper caption says "blue-white-red" but MeV renders as magenta/pink.
## Paper legend clearly shows z-score range -1.0 to 1.0.

make_paper_heatmap <- function(gene_table, fpkm_matrix, title, filename,
                               width = 6, height = NULL,
                               z_limits = c(-1, 1),
                               preserve_order = FALSE) {

  panel <- gene_table %>%
    inner_join(fpkm_matrix, by = "GeneID") %>%
    ## drop genes that are zero in ALL samples (unmapped)
    filter(if_any(all_of(sample_order), ~ .x > 0))

  if (nrow(panel) == 0) {
    cat("  WARNING: No genes found for '", title, "'\n")
    return(invisible(NULL))
  }

  ## z-score of log2(FPKM+1), per gene across the 4 samples
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

  ## Gene order: preserve input order if requested, otherwise sort by category
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
      ## MeV's "blue-white-red" actually renders as magenta/pink
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

## ---- 7. Gene panels --------------------------------------------------------

## 7a. Paper FA / lipid metabolism panel — EXACT order from the figure
##     Top-to-bottom as shown in the paper's heatmap
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

## 7a-extra. CD36 + LPL — mentioned in the paper text as cold-induced
##           but not in the heatmap figure; kept separate
fa_genes_extra <- tribble(
  ~Symbol,  ~GeneID,            ~Category,
  "CD36",   "ENSG00000135218",  "FA uptake (paper text)",
  "LPL",    "ENSG00000175445",  "FA uptake (paper text)"
)

## 7b. Core identity markers (trimmed — just the essentials)
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

## 7c. Claudin / tight junction family — broad scan
##     This lets you see which claudins actually have signal
claudin_genes <- tribble(
  ~Symbol,    ~GeneID,            ~Category,
  "CLDN1",    "ENSG00000163347",  "Claudin",
  "CLDN3",    "ENSG00000165215",  "Claudin",
  "CLDN4",    "ENSG00000189143",  "Claudin",
  "CLDN5",    "ENSG00000184113",  "Claudin",
  "CLDN6",    "ENSG00000184697",  "Claudin",
  "CLDN7",    "ENSG00000181885",  "Claudin",
  "CLDN8",    "ENSG00000156284",  "Claudin",
  "CLDN9",    "ENSG00000213937",  "Claudin",
  "CLDN10",   "ENSG00000134873",  "Claudin",
  "CLDN11",   "ENSG00000013297",  "Claudin",
  "CLDN12",   "ENSG00000107959",  "Claudin",
  "CLDN15",   "ENSG00000106404",  "Claudin",
  "CLDN16",   "ENSG00000113946",  "Claudin",
  "CLDN18",   "ENSG00000066405",  "Claudin",
  "CLDN19",   "ENSG00000164007",  "Claudin",
  "CLDN23",   "ENSG00000253958",  "Claudin",
  "OCLN",     "ENSG00000197822",  "TJ scaffold",
  "TJP1",     "ENSG00000104067",  "TJ scaffold",
  "TJP2",     "ENSG00000119139",  "TJ scaffold",
  "F11R",     "ENSG00000158769",  "TJ scaffold"
)

## ---- 8. Generate heatmaps --------------------------------------------------
cat("\n--- Generating heatmaps ---\n")

## Paper-exact: 18 genes in the paper's row order, [-1,1] clipping
p_fa <- make_paper_heatmap(
  fa_genes_paper, fpkm_wide,
  "FA / Lipid Metabolism (paper figure replication)",
  "heatmap_FA_metabolism.png", width = 6,
  preserve_order = TRUE
)

## CD36 + LPL extra (not in the figure but discussed in text)
p_fa_extra <- make_paper_heatmap(
  bind_rows(fa_genes_paper, fa_genes_extra), fpkm_wide,
  "FA / Lipid Metabolism + CD36/LPL",
  "heatmap_FA_metabolism_extended.png", width = 6,
  preserve_order = TRUE
)

p_id <- make_paper_heatmap(
  identity_genes, fpkm_wide,
  "BAT/WAT Identity Markers",
  "heatmap_identity_markers.png", width = 6
)

p_cldn <- make_paper_heatmap(
  claudin_genes, fpkm_wide,
  "Claudin / Tight Junction Family",
  "heatmap_claudin_family.png", width = 6
)

## Combined: identity + claudins on one figure
combined_genes <- bind_rows(identity_genes, claudin_genes)
p_combo <- make_paper_heatmap(
  combined_genes, fpkm_wide,
  "Identity Markers + Claudin Family",
  "heatmap_identity_and_claudins.png", width = 6.5
)

## ---- 9. Claudin family expression table ------------------------------------
cat("\n--- Claudin family expression summary (TPM) ---\n")

claudin_tpm <- claudin_genes %>%
  left_join(tpm_wide, by = "GeneID") %>%
  select(Symbol, Category, all_of(sample_order)) %>%
  mutate(
    max_TPM  = pmax(WAT_TN, WAT_CE, BAT_TN, BAT_CE, na.rm = TRUE),
    detected = max_TPM > 0.1
  ) %>%
  arrange(desc(max_TPM))

print(claudin_tpm %>% select(-detected), width = Inf)
cat("\n")

write_csv(claudin_tpm, file.path(RESULTS_DIR, "claudin_family_TPM.csv"))
cat("Saved: claudin_family_TPM.csv\n")

## ---- 10. CLDN1 raw evidence from each file ---------------------------------
cat("\n--- CLDN1 raw evidence (directly from Salmon files) ---\n")

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

print(cldn1_evidence, width = Inf)
cat("\n")
write_csv(cldn1_evidence, file.path(RESULTS_DIR, "CLDN1_evidence.csv"))

## ---- 11. Bar plots for detected claudins -----------------------------------
detected_symbols <- claudin_tpm %>% filter(detected) %>% pull(Symbol)

if (length(detected_symbols) > 0) {
  cldn_bar_df <- claudin_tpm %>%
    filter(Symbol %in% detected_symbols) %>%
    select(Symbol, all_of(sample_order)) %>%
    pivot_longer(cols = all_of(sample_order),
                 names_to = "sample", values_to = "TPM") %>%
    mutate(
      tissue    = ifelse(grepl("^BAT", sample), "BAT", "WAT"),
      condition = ifelse(grepl("_CE$", sample), "Cold", "TN"),
      sample    = factor(sample, levels = sample_order)
    )

  p_cldn_bar <- ggplot(cldn_bar_df,
                       aes(x = sample, y = TPM, fill = tissue)) +
    geom_col(width = 0.7) +
    facet_wrap(~ Symbol, scales = "free_y") +
    scale_fill_manual(values = c("BAT" = "#D6604D", "WAT" = "#4393C3")) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
      strip.text    = element_text(face = "bold.italic", size = 11),
      legend.position = "top"
    ) +
    labs(title = "Detected claudin / TJ genes (TPM > 0.1)",
         x = NULL, y = "TPM", fill = NULL)

  ggsave(file.path(RESULTS_DIR, "claudin_detected_barplots.png"),
         p_cldn_bar, width = 10, height = 7, dpi = 300)
  cat("Saved: claudin_detected_barplots.png\n")
}

## ---- 12. Positive controls bar plot ----------------------------------------
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
    subtitle = "UCP1 should be high in BAT-CE; LEP should be high in WAT",
    x = NULL, y = "TPM", fill = NULL
  )

ggsave(file.path(RESULTS_DIR, "positive_controls.png"), p_ctrl,
       width = 8, height = 5, dpi = 300)
cat("Saved: positive_controls.png\n")

## ---- 13. Descriptive fold-change table for FA panel -------------------------
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

## ---- 14. Summary ------------------------------------------------------------
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
cat("    heatmap_identity_markers.png       <- BAT/WAT/adipocyte markers\n")
cat("    heatmap_claudin_family.png         <- all claudins scanned\n")
cat("    heatmap_identity_and_claudins.png  <- combined\n")
cat("\n")
cat("  Claudin analysis:\n")
cat("    claudin_family_TPM.csv             <- TPM for every claudin\n")
cat("    claudin_detected_barplots.png      <- bar charts for detected ones\n")
cat("    CLDN1_evidence.csv                 <- raw Salmon evidence for CLDN1\n")
cat("\n")
cat("  Other:\n")
cat("    PCA_plot.png\n")
cat("    positive_controls.png\n")
cat("    FA_panel_FPKM_and_log2FC.csv\n")
cat("============================================================\n")
