#!/usr/bin/env python3
"""
03_compare_mof_vs_gwat.py
--------------------------
Compare MOF heterozygous (Mof+/-) visceral WAT bulk RNA-seq data from
GSE162653 to user's male gWAT adipocyte knockdown dataset.

This script:
  1. Loads processed DEG lists from both datasets
  2. Identifies overlapping DEGs (concordant and discordant)
  3. Performs gene set overlap analysis (Fisher's exact, hypergeometric)
  4. Generates comparison plots (overlap Venn, rank-rank, scatter, GSEA-style)
  5. Runs pathway enrichment on shared gene signatures
  6. Exports comparison tables

Input files:
  - MOF het DEG list:  results/mof_qc/deseq2_*.csv   (from 02_qc_mof_data.py)
                    OR  data/GSE162653/SuppData3_DEGs.xlsx  (from paper)
  - User gWAT DEGs:    User-provided file (CSV/TSV/XLSX with gene, log2FC, padj)

Usage:
    python3 scripts/03_compare_mof_vs_gwat.py \
        --mof-degs results/mof_qc/deseq2_Mofplusminus_vs_Mofplusplus.csv \
        --gwat-degs <path_to_your_gwat_degs.csv> \
        --out-dir results/comparison

    # Or use supplementary data directly:
    python3 scripts/03_compare_mof_vs_gwat.py \
        --mof-degs data/GSE162653/SuppData3_DEGs.xlsx \
        --gwat-degs <path_to_your_gwat_degs.csv> \
        --out-dir results/comparison
"""

import argparse
import os
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def load_deg_table(filepath, dataset_name="dataset"):
    """Load a DEG table from CSV, TSV, or XLSX format.

    Expects columns that can be mapped to: gene, log2FoldChange, padj.
    Returns a standardized DataFrame.
    """
    print(f"  Loading {dataset_name} DEGs from: {filepath}")

    if filepath.endswith(".xlsx") or filepath.endswith(".xls"):
        # Try to read first sheet; if multiple, let user know
        xls = pd.ExcelFile(filepath)
        print(f"    Excel sheets: {xls.sheet_names}")
        df = pd.read_excel(filepath, sheet_name=0)
    elif filepath.endswith(".csv"):
        df = pd.read_csv(filepath, index_col=0)
    else:
        # Try tab-separated
        df = pd.read_csv(filepath, sep="\t", index_col=0)

    print(f"    Raw shape: {df.shape}")
    print(f"    Columns: {list(df.columns)}")

    # Standardize column names
    col_map = {}
    for col in df.columns:
        col_lower = col.lower().strip()
        if col_lower in ("log2foldchange", "log2fc", "logfc", "lfc", "log2_fold_change"):
            col_map[col] = "log2FoldChange"
        elif col_lower in ("padj", "fdr", "adj.p.val", "p.adjust", "q_value", "qvalue",
                           "adjusted_pvalue", "adj_pval"):
            col_map[col] = "padj"
        elif col_lower in ("pvalue", "pval", "p.value", "p_value"):
            col_map[col] = "pvalue"
        elif col_lower in ("gene", "gene_name", "symbol", "gene_symbol", "geneid"):
            col_map[col] = "gene"
        elif col_lower in ("basemean", "averageexpression", "ave.expr", "aveexpr"):
            col_map[col] = "baseMean"

    df = df.rename(columns=col_map)

    # If gene is a column (not index), set it as index
    if "gene" in df.columns:
        df = df.set_index("gene")

    # Ensure required columns
    if "log2FoldChange" not in df.columns:
        print(f"    WARNING: No log2FoldChange column found in {dataset_name}")
        print(f"    Available columns: {list(df.columns)}")
        return None

    if "padj" not in df.columns and "pvalue" in df.columns:
        print("    No padj column; using raw pvalue as padj")
        df["padj"] = df["pvalue"]

    # Clean
    df = df.dropna(subset=["log2FoldChange"])
    df.index = df.index.astype(str).str.strip()
    df.index.name = "gene"

    print(f"    Standardized: {len(df)} genes with log2FC values")
    if "padj" in df.columns:
        sig = (df["padj"] < 0.05).sum()
        print(f"    Significant (padj < 0.05): {sig}")

    return df


def convert_mouse_genes(gene_list):
    """Ensure mouse gene names are in proper format (capitalize first letter).

    Mouse genes: Capitalize first letter, rest lowercase (e.g., Pparg, Glut4).
    Human genes: All uppercase (e.g., PPARG, GLUT4).
    """
    converted = []
    for g in gene_list:
        g = str(g).strip()
        if g.isupper() and len(g) > 1:
            # Likely human format -> convert to mouse
            converted.append(g[0].upper() + g[1:].lower())
        else:
            converted.append(g)
    return converted


# ---------------------------------------------------------------------------
# Comparison functions
# ---------------------------------------------------------------------------

def overlap_analysis(mof_degs, gwat_degs, fdr_thresh=0.05, out_dir="."):
    """Identify overlapping DEGs between the two datasets."""
    print("\n--- Overlap Analysis ---")

    # Significant DEGs in each
    mof_sig = mof_degs[mof_degs["padj"] < fdr_thresh].copy() if "padj" in mof_degs.columns else mof_degs
    gwat_sig = gwat_degs[gwat_degs["padj"] < fdr_thresh].copy() if "padj" in gwat_degs.columns else gwat_degs

    mof_genes = set(mof_sig.index)
    gwat_genes = set(gwat_sig.index)
    all_genes = set(mof_degs.index) | set(gwat_degs.index)
    overlap = mof_genes & gwat_genes

    print(f"  MOF het DEGs (FDR<{fdr_thresh}): {len(mof_genes)}")
    print(f"  gWAT KD DEGs (FDR<{fdr_thresh}): {len(gwat_genes)}")
    print(f"  Overlapping DEGs: {len(overlap)}")
    print(f"  Universe size: {len(all_genes)}")

    # Classify overlap direction
    if overlap:
        overlap_df = pd.DataFrame(index=sorted(overlap))
        overlap_df["mof_log2FC"] = mof_sig.loc[overlap_df.index, "log2FoldChange"].values
        overlap_df["gwat_log2FC"] = gwat_sig.loc[overlap_df.index, "log2FoldChange"].values
        overlap_df["same_direction"] = (
            np.sign(overlap_df["mof_log2FC"]) == np.sign(overlap_df["gwat_log2FC"])
        )
        if "padj" in mof_sig.columns:
            overlap_df["mof_padj"] = mof_sig.loc[overlap_df.index, "padj"].values
        if "padj" in gwat_sig.columns:
            overlap_df["gwat_padj"] = gwat_sig.loc[overlap_df.index, "padj"].values

        concordant = overlap_df["same_direction"].sum()
        discordant = (~overlap_df["same_direction"]).sum()
        print(f"  Concordant (same direction): {concordant}")
        print(f"  Discordant (opposite direction): {discordant}")

        overlap_df.to_csv(os.path.join(out_dir, "overlapping_degs.csv"))
    else:
        overlap_df = pd.DataFrame()

    # Fisher's exact test for enrichment of overlap
    a = len(overlap)  # both
    b = len(mof_genes - overlap)  # MOF only
    c = len(gwat_genes - overlap)  # gWAT only
    d = len(all_genes) - a - b - c  # neither
    contingency = [[a, b], [c, d]]
    odds_ratio, fisher_p = stats.fisher_exact(contingency, alternative="greater")
    print(f"\n  Fisher's exact test for overlap enrichment:")
    print(f"    Odds ratio: {odds_ratio:.2f}")
    print(f"    p-value: {fisher_p:.2e}")

    return mof_genes, gwat_genes, overlap, overlap_df, fisher_p


def plot_overlap(mof_genes, gwat_genes, overlap, fisher_p, out_dir):
    """Generate overlap visualization plots."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("MOF Het vs gWAT KD - DEG Overlap", fontsize=14, fontweight="bold")

    # Venn diagram
    ax = axes[0]
    try:
        v = venn2(
            [mof_genes, gwat_genes],
            set_labels=("MOF Het\n(Mof+/- WAT)", "gWAT KD\n(Your data)"),
            ax=ax,
        )
        ax.set_title(f"Overlap: {len(overlap)} genes\n(Fisher p={fisher_p:.2e})")
    except Exception:
        # Fallback if matplotlib_venn not available
        ax.text(0.5, 0.5, f"MOF: {len(mof_genes)}\ngWAT: {len(gwat_genes)}\nOverlap: {len(overlap)}",
                ha="center", va="center", fontsize=14, transform=ax.transAxes)
        ax.set_title("DEG Overlap")

    # Upset-style bar
    ax = axes[1]
    categories = ["MOF only", "Overlap", "gWAT only"]
    values = [len(mof_genes - overlap), len(overlap), len(gwat_genes - overlap)]
    colors = ["#4C72B0", "#55A868", "#DD8452"]
    ax.bar(categories, values, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_ylabel("Number of DEGs")
    ax.set_title("DEG Distribution")
    for i, v in enumerate(values):
        ax.text(i, v + max(values) * 0.02, str(v), ha="center", fontweight="bold")

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "05_deg_overlap.png"), dpi=200, bbox_inches="tight")
    plt.close()


def plot_fc_comparison(mof_degs, gwat_degs, overlap_df, out_dir):
    """Plot fold-change comparisons between datasets."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle("MOF Het vs gWAT KD - Fold Change Comparison", fontsize=14, fontweight="bold")

    # 1. Scatter of log2FC for overlapping DEGs
    ax = axes[0]
    if len(overlap_df) > 0:
        colors = ["#55A868" if sd else "#C44E52" for sd in overlap_df["same_direction"]]
        ax.scatter(overlap_df["mof_log2FC"], overlap_df["gwat_log2FC"],
                   c=colors, s=30, alpha=0.7, edgecolor="black", linewidth=0.3)
        # Add correlation
        r, p = stats.pearsonr(overlap_df["mof_log2FC"], overlap_df["gwat_log2FC"])
        ax.set_title(f"Overlapping DEGs (r={r:.2f}, p={p:.2e})")
        # Add diagonal
        lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lims, lims, "k--", alpha=0.3)
        ax.axhline(0, color="gray", linewidth=0.5)
        ax.axvline(0, color="gray", linewidth=0.5)
    ax.set_xlabel("MOF Het log2FC (Mof+/- vs Mof+/+)")
    ax.set_ylabel("gWAT KD log2FC")

    # 2. Rank-rank plot (all genes in common)
    ax = axes[1]
    common_genes = sorted(set(mof_degs.index) & set(gwat_degs.index))
    if len(common_genes) > 100:
        mof_ranks = mof_degs.loc[common_genes, "log2FoldChange"].rank()
        gwat_ranks = gwat_degs.loc[common_genes, "log2FoldChange"].rank()
        ax.scatter(mof_ranks, gwat_ranks, s=3, alpha=0.2, c="gray", rasterized=True)
        rho, p = stats.spearmanr(mof_ranks, gwat_ranks)
        ax.set_title(f"Rank-Rank Plot (Spearman rho={rho:.2f})")
    else:
        ax.set_title("Rank-Rank Plot (insufficient common genes)")
    ax.set_xlabel("MOF Het rank")
    ax.set_ylabel("gWAT KD rank")

    # 3. RRHO-style: running enrichment
    ax = axes[2]
    if len(common_genes) > 100:
        # Sort by MOF log2FC and compute running overlap with gWAT top genes
        mof_sorted = mof_degs.loc[common_genes].sort_values("log2FoldChange", ascending=False)
        gwat_top_up = set(gwat_degs.loc[common_genes].nlargest(
            max(int(len(common_genes) * 0.1), 50), "log2FoldChange"
        ).index)
        gwat_top_down = set(gwat_degs.loc[common_genes].nsmallest(
            max(int(len(common_genes) * 0.1), 50), "log2FoldChange"
        ).index)

        running_up = []
        running_down = []
        for i, gene in enumerate(mof_sorted.index, 1):
            running_up.append(len(gwat_top_up & set(mof_sorted.index[:i])) / max(i, 1))
            running_down.append(len(gwat_top_down & set(mof_sorted.index[:i])) / max(i, 1))

        x = np.arange(1, len(mof_sorted) + 1)
        ax.plot(x, running_up, color="#C44E52", label="gWAT Up", linewidth=1)
        ax.plot(x, running_down, color="#4C72B0", label="gWAT Down", linewidth=1)
        ax.axhline(len(gwat_top_up) / len(common_genes), color="gray", linestyle="--", alpha=0.5)
        ax.set_xlabel("MOF genes ranked by log2FC (high → low)")
        ax.set_ylabel("Fraction of gWAT top genes seen")
        ax.set_title("Running Enrichment")
        ax.legend(fontsize=8)
    else:
        ax.set_title("Running Enrichment (insufficient data)")

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "06_fc_comparison.png"), dpi=200, bbox_inches="tight")
    plt.close()


def pathway_overlap(overlap_df, out_dir):
    """Simple GO-like analysis of overlapping genes using keyword matching.

    For proper pathway analysis, use external tools like Enrichr or g:Profiler.
    This function exports gene lists formatted for those tools.
    """
    print("\n--- Preparing Gene Lists for Pathway Analysis ---")

    if len(overlap_df) == 0:
        print("  No overlapping DEGs to analyze.")
        return

    # Export gene lists for external tools
    all_overlap = sorted(overlap_df.index)
    concordant_up = sorted(
        overlap_df[(overlap_df["same_direction"]) & (overlap_df["mof_log2FC"] > 0)].index
    )
    concordant_down = sorted(
        overlap_df[(overlap_df["same_direction"]) & (overlap_df["mof_log2FC"] < 0)].index
    )
    discordant = sorted(overlap_df[~overlap_df["same_direction"]].index)

    gene_lists = {
        "all_overlap": all_overlap,
        "concordant_up": concordant_up,
        "concordant_down": concordant_down,
        "discordant": discordant,
    }

    for name, genes in gene_lists.items():
        if genes:
            path = os.path.join(out_dir, f"genelist_{name}.txt")
            with open(path, "w") as f:
                f.write("\n".join(genes))
            print(f"  {name}: {len(genes)} genes -> {path}")

    print("\n  Upload these gene lists to:")
    print("    - Enrichr: https://maayanlab.cloud/Enrichr/")
    print("    - g:Profiler: https://biit.cs.ut.ee/gprofiler/gost")
    print("    - DAVID: https://david.ncifcrf.gov/")
    print("  for pathway enrichment analysis.")


def export_comparison_summary(mof_degs, gwat_degs, overlap_df, fisher_p, out_dir):
    """Export a comprehensive comparison summary."""
    print("\n--- Exporting Comparison Summary ---")

    # Full merged table
    common_genes = sorted(set(mof_degs.index) & set(gwat_degs.index))
    merged = pd.DataFrame(index=common_genes)
    merged["mof_log2FC"] = mof_degs.loc[common_genes, "log2FoldChange"].values
    if "padj" in mof_degs.columns:
        merged["mof_padj"] = mof_degs.loc[common_genes, "padj"].values
    merged["gwat_log2FC"] = gwat_degs.loc[common_genes, "log2FoldChange"].values
    if "padj" in gwat_degs.columns:
        merged["gwat_padj"] = gwat_degs.loc[common_genes, "padj"].values

    # Classification
    merged["mof_sig"] = merged.get("mof_padj", pd.Series(dtype=float)).lt(0.05) if "mof_padj" in merged else False
    merged["gwat_sig"] = merged.get("gwat_padj", pd.Series(dtype=float)).lt(0.05) if "gwat_padj" in merged else False
    merged["both_sig"] = merged["mof_sig"] & merged["gwat_sig"]
    merged["concordant"] = np.sign(merged["mof_log2FC"]) == np.sign(merged["gwat_log2FC"])

    merged = merged.sort_values("mof_log2FC", ascending=False)
    merged.to_csv(os.path.join(out_dir, "full_comparison_table.csv"))

    # Summary stats
    summary_path = os.path.join(out_dir, "comparison_summary.txt")
    with open(summary_path, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("MOF Het vs gWAT KD Comparison Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"MOF het dataset: GSE162653 (Pessoa Rodrigues et al. 2021)\n")
        f.write(f"  Total genes: {len(mof_degs)}\n")
        if "padj" in mof_degs.columns:
            f.write(f"  DEGs (FDR<0.05): {(mof_degs['padj']<0.05).sum()}\n")
        f.write(f"\ngWAT KD dataset (user data):\n")
        f.write(f"  Total genes: {len(gwat_degs)}\n")
        if "padj" in gwat_degs.columns:
            f.write(f"  DEGs (FDR<0.05): {(gwat_degs['padj']<0.05).sum()}\n")
        f.write(f"\nComparison:\n")
        f.write(f"  Common genes: {len(common_genes)}\n")
        f.write(f"  Overlapping DEGs: {len(overlap_df)}\n")
        if len(overlap_df) > 0:
            f.write(f"  Concordant: {overlap_df['same_direction'].sum()}\n")
            f.write(f"  Discordant: {(~overlap_df['same_direction']).sum()}\n")
        f.write(f"  Fisher's exact p-value: {fisher_p:.2e}\n")

        # Correlation on common genes
        if len(common_genes) > 10:
            r, p = stats.pearsonr(
                mof_degs.loc[common_genes, "log2FoldChange"],
                gwat_degs.loc[common_genes, "log2FoldChange"],
            )
            rho, rho_p = stats.spearmanr(
                mof_degs.loc[common_genes, "log2FoldChange"],
                gwat_degs.loc[common_genes, "log2FoldChange"],
            )
            f.write(f"\n  Pearson r: {r:.3f} (p={p:.2e})\n")
            f.write(f"  Spearman rho: {rho:.3f} (p={rho_p:.2e})\n")

    print(f"  Summary saved to: {summary_path}")
    print(f"  Full table saved to: full_comparison_table.csv")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare MOF het WAT vs user gWAT KD bulk RNA-seq data"
    )
    parser.add_argument(
        "--mof-degs", required=True,
        help="Path to MOF het DEG results (CSV from 02_qc_mof_data.py or Supp Data XLSX)",
    )
    parser.add_argument(
        "--gwat-degs", required=True,
        help="Path to your gWAT KD DEG results (CSV/TSV/XLSX with gene, log2FC, padj)",
    )
    parser.add_argument(
        "--out-dir", default="results/comparison",
        help="Output directory (default: results/comparison)",
    )
    parser.add_argument(
        "--fdr", type=float, default=0.05,
        help="FDR threshold for calling DEGs (default: 0.05)",
    )
    parser.add_argument(
        "--convert-human", action="store_true",
        help="Convert human gene symbols to mouse format (PPARG -> Pparg)",
    )
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print("=" * 70)
    print("MOF Het vs gWAT KD - Cross-dataset Comparison")
    print("=" * 70)

    # Load datasets
    mof_degs = load_deg_table(args.mof_degs, "MOF het")
    gwat_degs = load_deg_table(args.gwat_degs, "gWAT KD")

    if mof_degs is None or gwat_degs is None:
        print("\nERROR: Could not load one or both DEG tables.")
        sys.exit(1)

    # Gene name harmonization
    if args.convert_human:
        gwat_degs.index = convert_mouse_genes(gwat_degs.index)

    common = set(mof_degs.index) & set(gwat_degs.index)
    print(f"\n  Common genes between datasets: {len(common)}")
    if len(common) == 0:
        print("  WARNING: No common genes found! Check gene name format.")
        print(f"  MOF example genes: {list(mof_degs.index[:5])}")
        print(f"  gWAT example genes: {list(gwat_degs.index[:5])}")
        print("  Try --convert-human if one dataset uses human gene names.")
        sys.exit(1)

    # Run comparison
    mof_genes, gwat_genes, overlap, overlap_df, fisher_p = overlap_analysis(
        mof_degs, gwat_degs, args.fdr, args.out_dir
    )

    # Plots
    try:
        plot_overlap(mof_genes, gwat_genes, overlap, fisher_p, args.out_dir)
        print("  Overlap plot saved.")
    except ImportError:
        print("  matplotlib_venn not installed; skipping Venn diagram.")
        print("  Install with: pip install matplotlib-venn")

    plot_fc_comparison(mof_degs, gwat_degs, overlap_df, args.out_dir)
    print("  FC comparison plots saved.")

    # Pathway analysis prep
    pathway_overlap(overlap_df, args.out_dir)

    # Export
    export_comparison_summary(mof_degs, gwat_degs, overlap_df, fisher_p, args.out_dir)

    print("\n" + "=" * 70)
    print("COMPARISON COMPLETE")
    print("=" * 70)
    print(f"Results in: {args.out_dir}/")
    print("  - overlapping_degs.csv")
    print("  - full_comparison_table.csv")
    print("  - comparison_summary.txt")
    print("  - genelist_*.txt (for pathway analysis)")
    print("  - Plots: 05_deg_overlap.png, 06_fc_comparison.png")


if __name__ == "__main__":
    main()
