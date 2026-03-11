#!/usr/bin/env python3
"""
02_qc_mof_data.py
-----------------
Quality control and exploratory analysis of the MOF heterozygous (Mof+/-)
bulk RNA-seq data from GSE162653.

This script:
  1. Loads the counts matrix (featureCounts output from GEO or user-provided)
  2. Loads sample metadata
  3. Performs library-level QC (total counts, detected genes, mapping stats)
  4. Generates QC plots (library size, gene detection, PCA, sample correlation)
  5. Runs differential expression with PyDESeq2
  6. Produces MA plots, volcano plots, and heatmaps of top DEGs
  7. Exports filtered/normalized data for downstream comparison

Input files expected in data/GSE162653/:
  - A counts matrix (tab-separated, genes x samples) — either downloaded from
    GEO or from the paper's supplementary data
  - sample_metadata.csv (from 01_download_mof_data.py or manually curated)

Usage:
    python3 scripts/02_qc_mof_data.py [--data-dir data/GSE162653] [--out-dir results/mof_qc]
"""

import argparse
import glob
import gzip
import os
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

warnings.filterwarnings("ignore", category=FutureWarning)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def find_counts_file(data_dir):
    """Auto-detect the counts matrix file in the data directory."""
    patterns = [
        "*.counts*.txt*",
        "*.count*.tsv*",
        "*featureCounts*",
        "*raw_counts*",
        "*gene_counts*",
        "*.txt",
        "*.tsv",
        "*.csv",
    ]
    candidates = []
    for pat in patterns:
        candidates.extend(glob.glob(os.path.join(data_dir, pat)))

    # Filter out metadata files
    candidates = [
        c for c in candidates
        if "metadata" not in os.path.basename(c).lower()
        and "series_matrix" not in os.path.basename(c).lower()
        and not c.endswith(".xlsx")
    ]

    if not candidates:
        return None

    # Prefer larger files (counts matrices are typically bigger)
    candidates.sort(key=lambda x: os.path.getsize(x), reverse=True)
    return candidates[0]


def load_counts(filepath):
    """Load a counts matrix from various formats."""
    print(f"  Loading counts from: {filepath}")

    if filepath.endswith(".gz"):
        opener = lambda: gzip.open(filepath, "rt")
    else:
        opener = lambda: open(filepath, "r")

    # Detect delimiter
    with opener() as f:
        # Skip comment lines
        header_line = ""
        for line in f:
            if not line.startswith("#"):
                header_line = line
                break

    sep = "\t" if "\t" in header_line else ","

    # Load with pandas
    df = pd.read_csv(filepath, sep=sep, comment="#", index_col=0)

    # If featureCounts format, it may have Chr, Start, End, Strand, Length columns
    meta_cols = ["Chr", "Start", "End", "Strand", "Length", "chr", "start",
                 "end", "strand", "length", "gene_name", "gene_biotype"]
    gene_info_cols = [c for c in df.columns if c in meta_cols]
    if gene_info_cols:
        gene_info = df[gene_info_cols].copy()
        df = df.drop(columns=gene_info_cols)
    else:
        gene_info = None

    # Ensure all remaining columns are numeric (sample columns)
    numeric_df = df.apply(pd.to_numeric, errors="coerce")
    non_numeric = numeric_df.columns[numeric_df.isna().all()]
    if len(non_numeric) > 0:
        print(f"  Warning: Dropping non-numeric columns: {list(non_numeric)}")
        numeric_df = numeric_df.drop(columns=non_numeric)

    print(f"  Loaded counts matrix: {numeric_df.shape[0]} genes x {numeric_df.shape[1]} samples")
    return numeric_df, gene_info


def load_or_create_metadata(data_dir, sample_names):
    """Load sample metadata or infer from sample names."""
    meta_path = os.path.join(data_dir, "sample_metadata.csv")

    if os.path.exists(meta_path):
        meta = pd.read_csv(meta_path)
        print(f"  Loaded metadata: {meta.shape[0]} samples, columns: {list(meta.columns)}")
        return meta

    # Try to infer metadata from sample names
    print("  No metadata file found. Inferring from sample names...")
    meta = pd.DataFrame({"sample": sample_names})

    # Common patterns in MOF study sample names
    genotype = []
    diet = []
    for s in sample_names:
        s_lower = s.lower()
        if "het" in s_lower or "mof+/-" in s_lower or "htz" in s_lower or "+/-" in s:
            genotype.append("Mof+/-")
        elif "wt" in s_lower or "mof+/+" in s_lower or "+/+" in s:
            genotype.append("Mof+/+")
        else:
            genotype.append("unknown")

        if "hfd" in s_lower or "high" in s_lower:
            diet.append("HFD")
        elif "sd" in s_lower or "chow" in s_lower or "normal" in s_lower:
            diet.append("SD")
        else:
            diet.append("unknown")

    meta["genotype"] = genotype
    meta["diet"] = diet
    meta["condition"] = meta["genotype"] + "_" + meta["diet"]

    meta.to_csv(meta_path, index=False)
    print(f"  Inferred metadata saved to: {meta_path}")
    print(f"  Please review and edit if needed!")
    return meta


# ---------------------------------------------------------------------------
# QC Functions
# ---------------------------------------------------------------------------

def library_qc(counts, out_dir):
    """Compute and plot library-level QC metrics."""
    print("\n--- Library-level QC ---")

    # Basic stats
    lib_sizes = counts.sum(axis=0)
    genes_detected = (counts > 0).sum(axis=0)
    genes_gt10 = (counts > 10).sum(axis=0)

    qc_df = pd.DataFrame({
        "sample": counts.columns,
        "total_counts": lib_sizes.values,
        "genes_detected": genes_detected.values,
        "genes_gt10_counts": genes_gt10.values,
        "mean_count_per_gene": counts.mean(axis=0).values,
        "median_count_per_gene": counts.median(axis=0).values,
    })

    # Print summary
    print(f"  Total counts range: {lib_sizes.min():,.0f} - {lib_sizes.max():,.0f}")
    print(f"  Mean library size: {lib_sizes.mean():,.0f}")
    print(f"  Genes detected (>0): {genes_detected.min()} - {genes_detected.max()}")
    print(f"  Genes with >10 counts: {genes_gt10.min()} - {genes_gt10.max()}")

    # Save QC table
    qc_path = os.path.join(out_dir, "library_qc_stats.csv")
    qc_df.to_csv(qc_path, index=False)
    print(f"  QC stats saved to: {qc_path}")

    # --- Plots ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("GSE162653 - MOF Het Bulk RNA-seq QC", fontsize=14, fontweight="bold")

    # 1. Library size barplot
    ax = axes[0, 0]
    colors = ["#4C72B0" if i % 2 == 0 else "#DD8452" for i in range(len(qc_df))]
    ax.bar(range(len(qc_df)), qc_df["total_counts"] / 1e6, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Library Size (millions)")
    ax.set_title("Library Size per Sample")
    ax.set_xticks(range(len(qc_df)))
    ax.set_xticklabels(qc_df["sample"], rotation=45, ha="right", fontsize=7)
    ax.axhline(y=lib_sizes.mean() / 1e6, color="red", linestyle="--", alpha=0.5, label="Mean")
    ax.legend(fontsize=8)

    # 2. Genes detected barplot
    ax = axes[0, 1]
    ax.bar(range(len(qc_df)), qc_df["genes_detected"], color=colors, edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Genes Detected (>0 counts)")
    ax.set_title("Genes Detected per Sample")
    ax.set_xticks(range(len(qc_df)))
    ax.set_xticklabels(qc_df["sample"], rotation=45, ha="right", fontsize=7)

    # 3. Count distribution (boxplot)
    ax = axes[1, 0]
    log_counts = np.log2(counts + 1)
    bp_data = [log_counts[col].values for col in log_counts.columns]
    bp = ax.boxplot(bp_data, tick_labels=counts.columns, patch_artist=True, showfliers=False)
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(colors[i] if i < len(colors) else "#4C72B0")
    ax.set_xlabel("Sample")
    ax.set_ylabel("log2(counts + 1)")
    ax.set_title("Count Distribution per Sample")
    ax.tick_params(axis="x", rotation=45, labelsize=7)

    # 4. Library size vs genes detected
    ax = axes[1, 1]
    ax.scatter(qc_df["total_counts"] / 1e6, qc_df["genes_detected"],
               s=80, edgecolor="black", linewidth=0.5, c="#4C72B0")
    for _, row in qc_df.iterrows():
        ax.annotate(row["sample"], (row["total_counts"] / 1e6, row["genes_detected"]),
                     fontsize=6, ha="center", va="bottom")
    ax.set_xlabel("Library Size (millions)")
    ax.set_ylabel("Genes Detected")
    ax.set_title("Library Size vs Genes Detected")

    plt.tight_layout()
    plot_path = os.path.join(out_dir, "01_library_qc.png")
    plt.savefig(plot_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  QC plot saved to: {plot_path}")

    return qc_df


def sample_similarity(counts, metadata, out_dir):
    """PCA, sample correlation, and hierarchical clustering."""
    print("\n--- Sample Similarity Analysis ---")

    # CPM normalization for visualization
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes, axis=1) * 1e6
    log_cpm = np.log2(cpm + 1)

    # Filter lowly expressed genes for PCA
    keep = (cpm > 1).sum(axis=1) >= 2  # expressed in at least 2 samples
    log_cpm_filt = log_cpm.loc[keep]
    print(f"  Genes passing filter (CPM>1 in >=2 samples): {keep.sum()}")

    # --- PCA ---
    from sklearn.decomposition import PCA

    pca = PCA(n_components=min(5, len(counts.columns)))
    pca_result = pca.fit_transform(log_cpm_filt.T)
    var_explained = pca.explained_variance_ratio_ * 100

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle("GSE162653 - Sample Similarity", fontsize=14, fontweight="bold")

    # PCA scatter
    ax = axes[0]
    # Color by condition if available
    if metadata is not None and "condition" in metadata.columns:
        conditions = metadata["condition"].values
        unique_conds = sorted(set(conditions))
        cmap = plt.cm.Set1
        color_map = {c: cmap(i / max(len(unique_conds) - 1, 1)) for i, c in enumerate(unique_conds)}
        colors = [color_map.get(c, "gray") for c in conditions]
        for cond in unique_conds:
            mask = [c == cond for c in conditions]
            ax.scatter(
                pca_result[mask, 0], pca_result[mask, 1],
                c=[color_map[cond]], label=cond, s=80, edgecolor="black", linewidth=0.5
            )
        ax.legend(fontsize=8, loc="best")
    else:
        ax.scatter(pca_result[:, 0], pca_result[:, 1], s=80, edgecolor="black", linewidth=0.5)

    for i, name in enumerate(counts.columns):
        ax.annotate(name, (pca_result[i, 0], pca_result[i, 1]),
                     fontsize=6, ha="center", va="bottom")
    ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)")
    ax.set_title("PCA of log2-CPM")

    # Scree plot
    ax = axes[1]
    ax.bar(range(1, len(var_explained) + 1), var_explained, color="#4C72B0", edgecolor="black")
    ax.set_xlabel("Principal Component")
    ax.set_ylabel("% Variance Explained")
    ax.set_title("Scree Plot")

    # Sample correlation heatmap
    ax = axes[2]
    corr = log_cpm_filt.corr(method="spearman")
    sns.heatmap(
        corr, ax=ax, cmap="RdYlBu_r", vmin=0.8, vmax=1.0,
        annot=True, fmt=".2f", square=True,
        xticklabels=counts.columns, yticklabels=counts.columns,
        annot_kws={"fontsize": 6},
    )
    ax.set_title("Spearman Correlation")
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    ax.tick_params(axis="y", rotation=0, labelsize=7)

    plt.tight_layout()
    plot_path = os.path.join(out_dir, "02_sample_similarity.png")
    plt.savefig(plot_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Sample similarity plot saved to: {plot_path}")

    # Save PCA coordinates
    pca_df = pd.DataFrame(
        pca_result,
        columns=[f"PC{i+1}" for i in range(pca_result.shape[1])],
        index=counts.columns,
    )
    pca_df.to_csv(os.path.join(out_dir, "pca_coordinates.csv"))

    # Save correlation matrix
    corr.to_csv(os.path.join(out_dir, "sample_correlation_spearman.csv"))

    return log_cpm_filt, corr


def run_deseq2(counts, metadata, contrast_col, contrast_ref, contrast_test, out_dir, extra_tag=""):
    """Run PyDESeq2 differential expression analysis."""
    print(f"\n--- Differential Expression: {contrast_test} vs {contrast_ref} ---")

    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        print("  PyDESeq2 not installed. Skipping DE analysis.")
        print("  Install with: pip install pydeseq2")
        return None

    # Prepare data
    # Filter samples for the contrast
    mask = metadata[contrast_col].isin([contrast_ref, contrast_test])
    meta_sub = metadata.loc[mask].copy()
    meta_sub.index = meta_sub["sample"]

    sample_cols = [s for s in meta_sub["sample"] if s in counts.columns]
    if len(sample_cols) < 4:
        print(f"  Not enough samples ({len(sample_cols)}) for DE analysis. Need >= 4.")
        return None

    counts_sub = counts[sample_cols].copy()
    meta_sub = meta_sub.loc[sample_cols]

    # Filter lowly expressed genes
    min_counts = 10
    keep = (counts_sub >= min_counts).sum(axis=1) >= 2
    counts_sub = counts_sub.loc[keep]
    print(f"  Samples: {len(sample_cols)} ({contrast_test}: {(meta_sub[contrast_col]==contrast_test).sum()}, "
          f"{contrast_ref}: {(meta_sub[contrast_col]==contrast_ref).sum()})")
    print(f"  Genes after filtering: {counts_sub.shape[0]}")

    # Ensure integer counts
    counts_sub = counts_sub.round(0).astype(int)

    # Run DESeq2
    dds = DeseqDataSet(
        counts=counts_sub.T,
        metadata=meta_sub[[contrast_col]],
        design=f"~{contrast_col}",
    )
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=[contrast_col, contrast_test, contrast_ref])
    stat_res.summary()

    results = stat_res.results_df
    results = results.sort_values("padj")

    # Save results
    tag = f"{contrast_test}_vs_{contrast_ref}".replace("+", "plus").replace("/", "").replace("-", "minus").replace(" ", "_")
    if extra_tag:
        tag = f"{tag}_{extra_tag}"
    res_path = os.path.join(out_dir, f"deseq2_{tag}.csv")
    results.to_csv(res_path)
    print(f"  DESeq2 results saved to: {res_path}")

    # Summary
    sig = results[results["padj"] < 0.05]
    up = sig[sig["log2FoldChange"] > 0]
    down = sig[sig["log2FoldChange"] < 0]
    print(f"  Significant DEGs (FDR < 0.05): {len(sig)}")
    print(f"    Upregulated in {contrast_test}: {len(up)}")
    print(f"    Downregulated in {contrast_test}: {len(down)}")

    # --- Plots ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"DESeq2: {contrast_test} vs {contrast_ref}", fontsize=14, fontweight="bold")

    # MA plot
    ax = axes[0]
    results["sig"] = results["padj"] < 0.05
    ax.scatter(
        results.loc[~results["sig"], "baseMean"],
        results.loc[~results["sig"], "log2FoldChange"],
        c="gray", alpha=0.3, s=5, rasterized=True,
    )
    ax.scatter(
        results.loc[results["sig"], "baseMean"],
        results.loc[results["sig"], "log2FoldChange"],
        c="red", alpha=0.5, s=8, rasterized=True,
    )
    ax.set_xscale("log")
    ax.set_xlabel("Mean Expression (baseMean)")
    ax.set_ylabel("log2 Fold Change")
    ax.set_title("MA Plot")
    ax.axhline(0, color="black", linewidth=0.5)

    # Volcano plot
    ax = axes[1]
    neg_log10_pval = -np.log10(results["padj"].clip(lower=1e-300))
    ax.scatter(
        results.loc[~results["sig"], "log2FoldChange"],
        neg_log10_pval[~results["sig"]],
        c="gray", alpha=0.3, s=5, rasterized=True,
    )
    ax.scatter(
        results.loc[results["sig"], "log2FoldChange"],
        neg_log10_pval[results["sig"]],
        c="red", alpha=0.5, s=8, rasterized=True,
    )
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.set_title("Volcano Plot")
    ax.axvline(0, color="black", linewidth=0.5)
    ax.axhline(-np.log10(0.05), color="blue", linewidth=0.5, linestyle="--")

    # Top DEGs heatmap
    ax = axes[2]
    top_n = min(30, len(sig))
    if top_n > 0:
        top_genes = sig.index[:top_n]
        # Use normalized counts from DESeq2
        norm_counts = pd.DataFrame(
            dds.layers["normed_counts"],
            index=dds.obs_names,
            columns=dds.var_names,
        )
        heatmap_data = norm_counts[top_genes].T
        # Z-score normalize
        z_scores = heatmap_data.subtract(heatmap_data.mean(axis=1), axis=0).div(
            heatmap_data.std(axis=1), axis=0
        )
        sns.heatmap(
            z_scores, ax=ax, cmap="RdBu_r", center=0,
            xticklabels=True, yticklabels=True,
            cbar_kws={"label": "Z-score"},
        )
        ax.set_title(f"Top {top_n} DEGs (Z-score)")
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.tick_params(axis="y", rotation=0, labelsize=6)
    else:
        ax.text(0.5, 0.5, "No significant DEGs", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Top DEGs Heatmap")

    plt.tight_layout()
    plot_path = os.path.join(out_dir, f"03_deseq2_{tag}.png")
    plt.savefig(plot_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  DE plots saved to: {plot_path}")

    return results


def gene_expression_qc(counts, out_dir):
    """Gene-level QC: expression distribution, biotype summary, etc."""
    print("\n--- Gene Expression QC ---")

    total_per_gene = counts.sum(axis=1)
    mean_per_gene = counts.mean(axis=1)
    detected_in = (counts > 0).sum(axis=1)

    # Expression bins
    not_expressed = (total_per_gene == 0).sum()
    low_expr = ((total_per_gene > 0) & (total_per_gene < 10)).sum()
    mid_expr = ((total_per_gene >= 10) & (total_per_gene < 1000)).sum()
    high_expr = (total_per_gene >= 1000).sum()

    print(f"  Total genes: {len(counts)}")
    print(f"  Not expressed (0 counts): {not_expressed}")
    print(f"  Low expression (<10 total): {low_expr}")
    print(f"  Medium expression (10-1000): {mid_expr}")
    print(f"  High expression (>=1000): {high_expr}")

    # Plot gene expression distribution
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("Gene Expression QC", fontsize=14, fontweight="bold")

    # Distribution of mean expression
    ax = axes[0]
    log_mean = np.log10(mean_per_gene[mean_per_gene > 0])
    ax.hist(log_mean, bins=50, color="#4C72B0", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("log10(mean counts)")
    ax.set_ylabel("Number of genes")
    ax.set_title("Gene Expression Distribution")
    ax.axvline(np.log10(10), color="red", linestyle="--", label="Mean=10")
    ax.legend()

    # Number of samples detecting each gene
    ax = axes[1]
    ax.hist(detected_in, bins=range(0, len(counts.columns) + 2),
            color="#DD8452", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Number of samples")
    ax.set_ylabel("Number of genes")
    ax.set_title("Gene Detection Across Samples")

    # Cumulative expression
    ax = axes[2]
    sorted_expr = total_per_gene.sort_values(ascending=False)
    cumsum = sorted_expr.cumsum() / sorted_expr.sum() * 100
    ax.plot(range(len(cumsum)), cumsum.values, color="#4C72B0", linewidth=1.5)
    ax.set_xlabel("Gene rank (by expression)")
    ax.set_ylabel("Cumulative % of total counts")
    ax.set_title("Cumulative Expression")
    # Mark how many genes account for 50% and 90%
    idx_50 = np.searchsorted(cumsum.values, 50)
    idx_90 = np.searchsorted(cumsum.values, 90)
    ax.axhline(50, color="red", linestyle="--", alpha=0.5)
    ax.axhline(90, color="red", linestyle="--", alpha=0.5)
    ax.axvline(idx_50, color="red", linestyle="--", alpha=0.5)
    ax.axvline(idx_90, color="red", linestyle="--", alpha=0.5)
    ax.text(idx_50, 52, f"{idx_50} genes", fontsize=8, color="red")
    ax.text(idx_90, 92, f"{idx_90} genes", fontsize=8, color="red")

    plt.tight_layout()
    plot_path = os.path.join(out_dir, "04_gene_expression_qc.png")
    plt.savefig(plot_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Gene expression QC plot saved to: {plot_path}")


def export_processed_data(counts, metadata, out_dir):
    """Export CPM-normalized and filtered data for downstream use."""
    print("\n--- Exporting Processed Data ---")

    # CPM normalization
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes, axis=1) * 1e6
    log_cpm = np.log2(cpm + 1)

    # Filter lowly expressed genes
    keep = (cpm > 1).sum(axis=1) >= 2
    cpm_filt = cpm.loc[keep]
    log_cpm_filt = log_cpm.loc[keep]
    counts_filt = counts.loc[keep]

    # Save
    counts_filt.to_csv(os.path.join(out_dir, "counts_filtered.csv"))
    cpm_filt.to_csv(os.path.join(out_dir, "cpm_filtered.csv"))
    log_cpm_filt.to_csv(os.path.join(out_dir, "log2cpm_filtered.csv"))

    print(f"  Exported {len(counts_filt)} genes (filtered from {len(counts)} total)")
    print(f"  Files: counts_filtered.csv, cpm_filtered.csv, log2cpm_filtered.csv")

    return counts_filt, cpm_filt, log_cpm_filt


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="QC and exploratory analysis of MOF het bulk RNA-seq (GSE162653)"
    )
    parser.add_argument(
        "--data-dir", default="data/GSE162653",
        help="Directory containing counts and metadata (default: data/GSE162653)",
    )
    parser.add_argument(
        "--out-dir", default="results/mof_qc",
        help="Output directory for results (default: results/mof_qc)",
    )
    parser.add_argument(
        "--counts", default=None,
        help="Path to counts matrix file (auto-detected if not provided)",
    )
    parser.add_argument(
        "--skip-de", action="store_true",
        help="Skip differential expression analysis",
    )
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print("=" * 70)
    print("MOF Het Bulk RNA-seq QC (GSE162653)")
    print("=" * 70)

    # --- Load counts ---
    counts_path = args.counts or find_counts_file(args.data_dir)
    if counts_path is None:
        print(
            "\nERROR: No counts file found in", args.data_dir,
            "\nPlease either:",
            "\n  1. Run 01_download_mof_data.py first",
            "\n  2. Place counts matrix in", args.data_dir,
            "\n  3. Specify --counts <path>",
        )
        sys.exit(1)

    counts, gene_info = load_counts(counts_path)

    # --- Load metadata ---
    metadata = load_or_create_metadata(args.data_dir, list(counts.columns))

    # --- QC ---
    qc_df = library_qc(counts, args.out_dir)
    gene_expression_qc(counts, args.out_dir)
    log_cpm_filt, corr = sample_similarity(counts, metadata, args.out_dir)

    # --- DE analysis ---
    if not args.skip_de:
        # Try genotype comparison on SD samples
        if "genotype" in metadata.columns and "diet" in metadata.columns:
            for diet_val in metadata["diet"].unique():
                if diet_val == "unknown":
                    continue
                meta_diet = metadata[metadata["diet"] == diet_val].copy()
                meta_diet = meta_diet.reset_index(drop=True)
                if len(meta_diet["genotype"].unique()) >= 2:
                    run_deseq2(
                        counts, meta_diet, "genotype", "Mof+/+", "Mof+/-", args.out_dir,
                        extra_tag=diet_val,
                    )
        elif "condition" in metadata.columns:
            conditions = metadata["condition"].unique()
            if len(conditions) >= 2:
                run_deseq2(
                    counts, metadata, "condition", conditions[0], conditions[1], args.out_dir
                )

    # --- Export processed data ---
    export_processed_data(counts, metadata, args.out_dir)

    print("\n" + "=" * 70)
    print("QC COMPLETE")
    print("=" * 70)
    print(f"Results saved to: {args.out_dir}")
    print(f"  - library_qc_stats.csv")
    print(f"  - pca_coordinates.csv")
    print(f"  - sample_correlation_spearman.csv")
    print(f"  - counts_filtered.csv / cpm_filtered.csv / log2cpm_filtered.csv")
    if not args.skip_de:
        print(f"  - deseq2_*.csv (DE results)")
    print(f"  - Plots: 01-04_*.png")


if __name__ == "__main__":
    main()
