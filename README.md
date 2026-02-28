# CLDN1 in human BAT after cold exposure (E-MTAB-4031)

This repo contains an **R script** to test whether **CLDN1** is increased in **human BAT** after cold exposure using ArrayExpress study **E-MTAB-4031**.

## What this script does

- Reads a processed expression matrix and SDRF metadata.
- Runs sanity checks so you can understand data structure before analysis.
- Auto-infers sample groups (BAT/WAT and COLD/CONTROL), or accepts a manual sample map.
- Extracts CLDN1 values in BAT and compares COLD vs CONTROL.
- Produces output folders with plots, CSVs, and run logs.

## 1) Download files from ArrayExpress

Study page: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-4031

Download at least:

- Processed expression matrix (`E-MTAB-4031.processed.*.zip`)
- SDRF file (`E-MTAB-4031.sdrf.txt`)

Optional (if gene symbols are missing from expression matrix):

- Platform annotation table with `probe_id` and `gene_symbol`

Suggested local structure:

```bash
mkdir -p data/E-MTAB-4031
# unzip downloaded files into data/E-MTAB-4031
```

## 2) Install R packages

```r
install.packages(c("optparse", "readr", "dplyr", "tidyr", "stringr", "ggplot2"))
```

## 3) Run analysis

```bash
Rscript scripts/analyze_cldn1_bat_cold.R \
  --expression "data/E-MTAB-4031/<your_processed_matrix>.txt" \
  --sdrf "data/E-MTAB-4031/E-MTAB-4031.sdrf.txt" \
  --out-dir "results/cldn1_bat_cold"
```

If expression lacks gene symbols:

```bash
Rscript scripts/analyze_cldn1_bat_cold.R \
  --expression "data/E-MTAB-4031/<your_processed_matrix>.txt" \
  --sdrf "data/E-MTAB-4031/E-MTAB-4031.sdrf.txt" \
  --platform-annotation "data/E-MTAB-4031/<annotation_file>.txt" \
  --out-dir "results/cldn1_bat_cold"
```

## 4) Output folders and files

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

## 5) First-time ArrayExpress sanity workflow

1. Run once using auto-inference.
2. Open `run_summary.txt` and verify:
   - expression and SDRF dimensions
   - detected sample columns
   - detected SDRF sample-ID column
   - BAT and COLD/CONTROL counts
3. Open `inferred_or_input_sample_map.csv` and verify group labels per sample.
4. If labels are wrong, create a manual sample map and rerun with `--sample-map`.

### Manual sample-map format

CSV/TSV with columns:

- `sample_id`
- `tissue_group` (for example `BAT`)
- `temperature_group` (`COLD` or `CONTROL`)

Example:

```csv
sample_id,tissue_group,temperature_group
Sample_01,BAT,COLD
Sample_02,BAT,CONTROL
```

Run with manual map:

```bash
Rscript scripts/analyze_cldn1_bat_cold.R \
  --expression "data/E-MTAB-4031/<your_processed_matrix>.txt" \
  --sdrf "data/E-MTAB-4031/E-MTAB-4031.sdrf.txt" \
  --sample-map "data/E-MTAB-4031/my_sample_map.csv" \
  --out-dir "results/cldn1_bat_cold"
```

## Interpretation notes

- `cldn1_stats.csv` includes `log2_fc_cold_minus_control` as mean(COLD) - mean(CONTROL).
- Positive values suggest higher CLDN1 in cold-exposed BAT.
- Treat p-values carefully if sample sizes are small.

## Reference paper

https://www.cell.com/cell-metabolism/fulltext/S1550-4131(16)30185-1
