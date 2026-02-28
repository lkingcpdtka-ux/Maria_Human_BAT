# CLDN1 in human BAT after cold exposure (E-MTAB-4031)

Yes — for this dataset type, you should generally use **limma** for differential expression.

This workflow supports both:
1. **Metadata sanity checks** (IDF/SDRF only), and
2. **Full expression analysis** with:
   - simple Welch tests (quick check), and
   - **limma** (recommended primary analysis).

## Why your previous run failed

You used `source("...Main_Analysis.R")`, which executes the script body in your interactive session and expects CLI args like `--sdrf`.

The script is now refactored so you can either:
- run with `Rscript ... --args`, or
- `source()` it and call `main(c(...))` explicitly.

## 1) Install packages

CRAN packages:

```r
install.packages(c("optparse", "readr", "dplyr", "tidyr", "stringr", "ggplot2"))
```

Bioconductor (`limma`):

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("limma")
```

## 2) Run metadata-only sanity checks (with your current files)

```bash
Rscript scripts/analyze_cldn1_bat_cold.R \
  --idf "data/E-MTAB-4031/E-MTAB-4031.idf.txt" \
  --sdrf "data/E-MTAB-4031/E-MTAB-4031.sdrf.txt" \
  --metadata-only \
  --out-dir "results/cldn1_bat_cold"
```

Key sanity outputs:
- `sdrf_column_summary.csv`
- `sdrf_tissue_related_preview.csv` (if detected)
- `sdrf_temperature_related_preview.csv` (if detected)
- `sdrf_file_references.csv`
- `idf_file_references.csv`
- `logs/run_summary.txt`

## 3) Run full analysis (after expression matrix download)

```bash
Rscript scripts/analyze_cldn1_bat_cold.R \
  --idf "data/E-MTAB-4031/E-MTAB-4031.idf.txt" \
  --sdrf "data/E-MTAB-4031/E-MTAB-4031.sdrf.txt" \
  --expression "data/E-MTAB-4031/<processed_expression_matrix>.txt" \
  --out-dir "results/cldn1_bat_cold"
```

Optional:
- `--platform-annotation <file>` if gene symbols are absent in expression matrix.
- `--sample-map <file>` to override auto-inferred BAT/COLD labels.

## 4) limma outputs (recommended)

When `limma` is installed and both COLD/CONTROL BAT groups are present, the script writes:
- `bat_cold_vs_control_limma_all_genes.csv`
- `cldn1_limma_stats.csv`

These are your preferred inferential outputs over simple per-gene t-tests.

## 5) If you want to use `source()` interactively

```r
source("scripts/analyze_cldn1_bat_cold.R")
main(c(
  "--idf", "data/E-MTAB-4031/E-MTAB-4031.idf.txt",
  "--sdrf", "data/E-MTAB-4031/E-MTAB-4031.sdrf.txt",
  "--metadata-only",
  "--out-dir", "results/cldn1_bat_cold"
))
```

## References
- Study: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-4031
- Paper: https://www.cell.com/cell-metabolism/fulltext/S1550-4131(16)30185-1
