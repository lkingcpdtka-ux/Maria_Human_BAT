# CLDN1 in human BAT after cold exposure (E-MTAB-4031)

This repo now has an **R workflow with two stages**:

1. **Metadata sanity-check stage** (works with the files you already downloaded: `idf` + `sdrf`).
2. **Expression analysis stage** (runs once you also download the processed expression matrix).

---

## Files you currently have

From your screenshot, you downloaded:

- `E-MTAB-4031.idf.txt`
- `E-MTAB-4031.sdrf.txt`

That is enough to run **metadata-only sanity checks** and produce structured reports about how the dataset is organized.

## 1) Install R packages

```r
install.packages(c("optparse", "readr", "dplyr", "tidyr", "stringr", "ggplot2"))
```

## 2) Run metadata-only sanity checks (with your current files)

```bash
Rscript scripts/analyze_cldn1_bat_cold.R \
  --idf "data/E-MTAB-4031/E-MTAB-4031.idf.txt" \
  --sdrf "data/E-MTAB-4031/E-MTAB-4031.sdrf.txt" \
  --metadata-only \
  --out-dir "results/cldn1_bat_cold"
```

This writes:

- `results/cldn1_bat_cold/csv/sdrf_column_summary.csv`
- `results/cldn1_bat_cold/csv/sdrf_tissue_related_preview.csv` (if tissue-related columns exist)
- `results/cldn1_bat_cold/csv/sdrf_temperature_related_preview.csv` (if condition columns exist)
- `results/cldn1_bat_cold/csv/sdrf_file_references.csv` (all file names referenced in SDRF)
- `results/cldn1_bat_cold/csv/idf_file_references.csv` (file-like entries found in IDF)
- `results/cldn1_bat_cold/logs/run_summary.txt`

Use these files to understand data structure and identify what additional data file(s) you still need.

## 3) Download missing expression file(s)

Use `idf_file_references.csv` and `sdrf_file_references.csv` to identify the processed expression matrix filename(s) referenced by the study.

After you download/unzip the processed matrix, run full analysis.

## 4) Run full CLDN1 analysis

```bash
Rscript scripts/analyze_cldn1_bat_cold.R \
  --idf "data/E-MTAB-4031/E-MTAB-4031.idf.txt" \
  --sdrf "data/E-MTAB-4031/E-MTAB-4031.sdrf.txt" \
  --expression "data/E-MTAB-4031/<processed_expression_matrix>.txt" \
  --out-dir "results/cldn1_bat_cold"
```

If gene symbols are not in the expression matrix, add:

```bash
  --platform-annotation "data/E-MTAB-4031/<annotation_file>.txt"
```

If automatic sample grouping is wrong, provide your own mapping:

```bash
  --sample-map "data/E-MTAB-4031/my_sample_map.csv"
```

Manual sample map columns:

- `sample_id`
- `tissue_group` (e.g., `BAT`)
- `temperature_group` (`COLD` or `CONTROL`)

## 5) Full-analysis outputs

- `results/cldn1_bat_cold/plots/`
  - `cldn1_bat_cold_vs_control.png`
- `results/cldn1_bat_cold/csv/`
  - `inferred_or_input_sample_map.csv`
  - `merged_expression_metadata_preview.csv`
  - `bat_samples_used.csv`
  - `cldn1_sample_values.csv`
  - `cldn1_group_summary.csv`
  - `cldn1_stats.csv`
  - `bat_cold_vs_control_all_genes.csv`
- `results/cldn1_bat_cold/logs/`
  - `run_summary.txt`

## Interpretation notes

- `cldn1_stats.csv` reports `log2_fc_cold_minus_control = mean(COLD) - mean(CONTROL)`.
- Positive values suggest higher CLDN1 in cold-exposed BAT.
- Check sample counts and map quality first; small n can make p-values unstable.

## Reference

- Study: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-4031
- Paper: https://www.cell.com/cell-metabolism/fulltext/S1550-4131(16)30185-1
