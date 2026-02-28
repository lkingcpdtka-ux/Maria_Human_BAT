#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

read_tabular <- function(path) {
  readr::read_delim(
    file = path,
    delim = NULL,
    comment = "#",
    show_col_types = FALSE,
    progress = FALSE,
    col_types = cols(.default = col_character())
  )
}

ensure_output_dirs <- function(out_dir) {
  dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "csv"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
}

find_expression_sample_columns <- function(df) {
  frac_numeric <- sapply(df, function(col) {
    suppressWarnings(mean(!is.na(as.numeric(col))))
  })
  names(frac_numeric[frac_numeric > 0.8])
}

find_probe_id_column <- function(df, sample_cols) {
  non_sample <- setdiff(colnames(df), sample_cols)
  if (length(non_sample) > 0) non_sample[[1]] else colnames(df)[[1]]
}

find_gene_symbol_column <- function(df) {
  cols <- colnames(df)
  idx <- str_which(str_to_lower(cols), "gene symbol|symbol|hgnc")
  if (length(idx) > 0) cols[[idx[[1]]]] else NA_character_
}

find_best_sdrf_sample_column <- function(sdrf, expr_samples) {
  overlaps <- sapply(sdrf, function(col) sum(col %in% expr_samples))
  if (max(overlaps) == 0) return(NA_character_)
  names(overlaps)[which.max(overlaps)]
}

infer_groups <- function(sdrf, sample_col) {
  collapsed <- apply(replace_na(sdrf, ""), 1, paste, collapse = " | ") |> str_to_lower()

  tissue_group <- ifelse(
    str_detect(collapsed, "\\bbat\\b|brown adipose|supraclavicular"), "BAT",
    ifelse(str_detect(collapsed, "\\bwat\\b|white adipose"), "WAT", "UNKNOWN")
  )

  temperature_group <- ifelse(
    str_detect(collapsed, "cold|cooling|16c|17c|18c|19c"), "COLD",
    ifelse(
      str_detect(collapsed, "control|thermoneutral|warm|24c|baseline|room temp"),
      "CONTROL",
      "UNKNOWN"
    )
  )

  tibble(
    sample_id = as.character(sdrf[[sample_col]]),
    tissue_group = tissue_group,
    temperature_group = temperature_group,
    inference_source = "auto_from_sdrf_text"
  )
}

option_list <- list(
  make_option("--expression", type = "character", help = "Path to processed expression matrix (txt/tsv)."),
  make_option("--sdrf", type = "character", help = "Path to SDRF sample annotation file."),
  make_option("--sample-map", type = "character", default = NULL,
              help = "Optional CSV/TSV with columns: sample_id,tissue_group,temperature_group."),
  make_option("--platform-annotation", type = "character", default = NULL,
              help = "Optional annotation file with columns probe_id and gene_symbol."),
  make_option("--out-dir", type = "character", default = "results/cldn1_bat_cold",
              help = "Output folder [default %default].")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$expression) || is.null(opt$sdrf)) {
  stop("--expression and --sdrf are required.")
}

out_dir <- opt$out_dir
ensure_output_dirs(out_dir)

log_lines <- c("=== CLDN1 BAT cold analysis: sanity checks ===")

expr <- read_tabular(opt$expression)
sdrf <- read_tabular(opt$sdrf)

log_lines <- c(log_lines,
               sprintf("Expression matrix shape: (%d, %d)", nrow(expr), ncol(expr)),
               sprintf("SDRF shape: (%d, %d)", nrow(sdrf), ncol(sdrf)),
               sprintf("Expression columns (first 12): %s", paste(head(colnames(expr), 12), collapse = ", ")),
               sprintf("SDRF columns: %s", paste(colnames(sdrf), collapse = ", ")))

sample_cols <- find_expression_sample_columns(expr)
if (length(sample_cols) < 2) {
  stop("Could not identify sample columns in expression matrix. Check format and delimiters.")
}

probe_col <- find_probe_id_column(expr, sample_cols)
symbol_col <- find_gene_symbol_column(expr)

log_lines <- c(log_lines,
               sprintf("Detected expression sample columns: %d", length(sample_cols)),
               sprintf("Detected probe ID column: %s", probe_col),
               sprintf("Detected gene symbol column: %s", symbol_col))

sdrf_sample_col <- find_best_sdrf_sample_column(sdrf, sample_cols)
if (is.na(sdrf_sample_col)) {
  stop("No SDRF column overlaps expression sample names. Provide --sample-map.")
}
log_lines <- c(log_lines, sprintf("Detected SDRF sample ID column: %s", sdrf_sample_col))

if (!is.null(opt$`sample-map`)) {
  sample_map <- read_tabular(opt$`sample-map`)
  required <- c("sample_id", "tissue_group", "temperature_group")
  if (!all(required %in% colnames(sample_map))) {
    stop("Sample map is missing required columns: sample_id,tissue_group,temperature_group")
  }
  sample_map <- sample_map |> mutate(inference_source = "user_provided")
  log_lines <- c(log_lines, "Using user-provided sample map.")
} else {
  sample_map <- infer_groups(sdrf, sdrf_sample_col)
  log_lines <- c(log_lines, "Using auto-inferred sample groups from SDRF text.")
}

write_csv(sample_map, file.path(out_dir, "csv", "inferred_or_input_sample_map.csv"))

map_clean <- sample_map |>
  filter(sample_id %in% sample_cols) |>
  distinct(sample_id, .keep_all = TRUE)

overlap <- nrow(map_clean)
log_lines <- c(log_lines,
               sprintf("Sample map vs expression overlap: %d samples", overlap),
               sprintf("Counts by tissue_group: %s",
                       paste(capture.output(print(table(map_clean$tissue_group, useNA = "ifany"))), collapse = " ")),
               sprintf("Counts by temperature_group: %s",
                       paste(capture.output(print(table(map_clean$temperature_group, useNA = "ifany"))), collapse = " ")))

expr_long <- expr |>
  select(all_of(c(probe_col, sample_cols))) |>
  pivot_longer(cols = all_of(sample_cols), names_to = "sample_id", values_to = "expression") |>
  mutate(expression = suppressWarnings(as.numeric(expression)))

if (!is.na(symbol_col)) {
  symbols <- expr |>
    select(all_of(c(probe_col, symbol_col))) |>
    distinct() |>
    rename(gene_symbol = !!sym(symbol_col))
} else if (!is.null(opt$`platform-annotation`)) {
  annot <- read_tabular(opt$`platform-annotation`)
  if (!all(c("probe_id", "gene_symbol") %in% colnames(annot))) {
    stop("Platform annotation needs columns: probe_id,gene_symbol")
  }
  symbols <- annot |>
    select(probe_id, gene_symbol) |>
    rename(!!probe_col := probe_id)
} else {
  stop("No gene symbol column found and no --platform-annotation provided.")
}

merged <- expr_long |>
  left_join(symbols, by = probe_col) |>
  inner_join(map_clean, by = "sample_id")

write_csv(merged, file.path(out_dir, "csv", "merged_expression_metadata_preview.csv"))

bat <- merged |>
  filter(str_to_upper(tissue_group) == "BAT", str_to_upper(temperature_group) %in% c("COLD", "CONTROL"))

if (nrow(bat) == 0) {
  stop("No BAT samples with temperature_group in {COLD, CONTROL}. Fix sample map and rerun.")
}

bat_samples <- bat |>
  distinct(sample_id, temperature_group) |>
  arrange(temperature_group, sample_id)

write_csv(bat_samples, file.path(out_dir, "csv", "bat_samples_used.csv"))

log_lines <- c(log_lines,
               sprintf("BAT analysis samples: %d", nrow(bat_samples)),
               sprintf("BAT group counts: %s",
                       paste(capture.output(print(table(bat_samples$temperature_group))), collapse = " ")))

cldn1 <- bat |>
  filter(str_to_upper(gene_symbol) == "CLDN1")

if (nrow(cldn1) == 0) {
  stop("CLDN1 not found. Check annotation / mapping.")
}

cldn1_vals <- cldn1 |>
  select(all_of(c(probe_col, "gene_symbol", "sample_id", "temperature_group", "expression")))
write_csv(cldn1_vals, file.path(out_dir, "csv", "cldn1_sample_values.csv"))

cldn1_summary <- cldn1 |>
  group_by(temperature_group) |>
  summarise(
    count = sum(!is.na(expression)),
    mean = mean(expression, na.rm = TRUE),
    median = median(expression, na.rm = TRUE),
    sd = sd(expression, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(cldn1_summary, file.path(out_dir, "csv", "cldn1_group_summary.csv"))

cold_vals <- cldn1 |>
  filter(str_to_upper(temperature_group) == "COLD") |>
  pull(expression) |>
  na.omit()
ctrl_vals <- cldn1 |>
  filter(str_to_upper(temperature_group) == "CONTROL") |>
  pull(expression) |>
  na.omit()

if (length(cold_vals) >= 2 && length(ctrl_vals) >= 2) {
  t_res <- t.test(cold_vals, ctrl_vals)
  t_stat <- unname(t_res$statistic)
  p_value <- t_res$p.value
} else {
  t_stat <- NA_real_
  p_value <- NA_real_
  log_lines <- c(log_lines, "WARNING: <2 values in one group for CLDN1; t-test not computed.")
}

cldn1_stats <- tibble(
  cold_n = length(cold_vals),
  control_n = length(ctrl_vals),
  cold_mean = ifelse(length(cold_vals) > 0, mean(cold_vals), NA_real_),
  control_mean = ifelse(length(ctrl_vals) > 0, mean(ctrl_vals), NA_real_),
  log2_fc_cold_minus_control = ifelse(length(cold_vals) > 0 && length(ctrl_vals) > 0,
                                      mean(cold_vals) - mean(ctrl_vals), NA_real_),
  welch_t_stat = t_stat,
  p_value = p_value
)
write_csv(cldn1_stats, file.path(out_dir, "csv", "cldn1_stats.csv"))

p <- ggplot(cldn1, aes(x = temperature_group, y = expression)) +
  geom_boxplot(fill = "#a6cee3") +
  geom_jitter(width = 0.12, color = "#1f78b4", size = 2) +
  labs(title = "CLDN1 expression in human BAT: COLD vs CONTROL",
       x = "Temperature group", y = "Expression") +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "plots", "cldn1_bat_cold_vs_control.png"), plot = p, width = 7, height = 5, dpi = 300)

all_gene_stats <- bat |>
  filter(!is.na(gene_symbol)) |>
  group_by(gene_symbol) |>
  group_modify(~{
    v_cold <- .x |> filter(str_to_upper(temperature_group) == "COLD") |> pull(expression) |> na.omit()
    v_ctrl <- .x |> filter(str_to_upper(temperature_group) == "CONTROL") |> pull(expression) |> na.omit()
    if (length(v_cold) < 2 || length(v_ctrl) < 2) {
      return(tibble())
    }
    t_res <- t.test(v_cold, v_ctrl)
    tibble(
      cold_mean = mean(v_cold),
      control_mean = mean(v_ctrl),
      log2_fc_cold_minus_control = mean(v_cold) - mean(v_ctrl),
      welch_t_stat = unname(t_res$statistic),
      p_value = t_res$p.value
    )
  }) |>
  ungroup()

if (nrow(all_gene_stats) > 0) {
  all_gene_stats <- all_gene_stats |>
    mutate(fdr_bh = p.adjust(p_value, method = "BH")) |>
    arrange(p_value)
}

write_csv(all_gene_stats, file.path(out_dir, "csv", "bat_cold_vs_control_all_genes.csv"))

log_lines <- c(log_lines, "Output files written under:", normalizePath(out_dir))
writeLines(log_lines, file.path(out_dir, "logs", "run_summary.txt"))

cat("Analysis complete. See run summary:", file.path(out_dir, "logs", "run_summary.txt"), "\n")
