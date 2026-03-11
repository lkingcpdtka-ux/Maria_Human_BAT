#!/usr/bin/env python3
"""
Generate synthetic test data to validate the QC and comparison pipeline.
This creates realistic-looking bulk RNA-seq count matrices and DEG tables
mimicking the MOF het WAT and user gWAT KD datasets.

Usage:
    python3 scripts/generate_test_data.py
"""

import os
import numpy as np
import pandas as pd

np.random.seed(42)

# ---------------------------------------------------------------------------
# Generate MOF het counts matrix (mimicking GSE162653 visceral WAT)
# ---------------------------------------------------------------------------
n_genes = 20000
n_samples = 16  # 4 conditions x 4 replicates

# Sample names: Mof+/+ SD, Mof+/+ HFD, Mof+/- SD, Mof+/- HFD
sample_names = []
genotypes = []
diets = []
for gt in ["WT", "Het"]:
    for diet in ["SD", "HFD"]:
        for rep in range(1, 5):
            name = f"Mof_{gt}_{diet}_rep{rep}"
            sample_names.append(name)
            genotypes.append("Mof+/+" if gt == "WT" else "Mof+/-")
            diets.append(diet)

# Gene names (mouse format)
gene_prefixes = ["Pparg", "Glut4", "Fabp4", "Adipoq", "Lep", "Cebpa", "Ucp1",
                  "Prdm16", "Cidea", "Elovl3", "Cox7a1", "Dio2"]
gene_names = gene_prefixes.copy()
for i in range(n_genes - len(gene_prefixes)):
    gene_names.append(f"Gene{i+1:05d}")

# Base expression (log-normal distribution)
base_expr = np.random.lognormal(mean=3, sigma=2, size=n_genes)
base_expr = np.clip(base_expr, 0, 50000)

# Generate counts with biological variation
counts = np.zeros((n_genes, n_samples), dtype=int)
for j in range(n_samples):
    # Add sample-specific noise
    sample_factor = np.random.normal(1, 0.1)
    lib_size_factor = np.random.uniform(0.8, 1.5)

    for i in range(n_genes):
        mu = base_expr[i] * sample_factor * lib_size_factor

        # Add genotype effect for specific genes (MOF targets)
        if genotypes[j] == "Mof+/-" and i < 500:
            mu *= np.random.uniform(0.3, 0.8)  # downregulated in het
        elif genotypes[j] == "Mof+/-" and 500 <= i < 800:
            mu *= np.random.uniform(1.2, 2.0)  # upregulated in het

        # Add diet effect
        if diets[j] == "HFD" and 200 <= i < 600:
            mu *= np.random.uniform(1.3, 2.5)

        counts[i, j] = max(0, int(np.random.negative_binomial(
            n=max(1, int(mu / 5)), p=0.5
        )))

counts_df = pd.DataFrame(counts, index=gene_names, columns=sample_names)
counts_df.index.name = "Geneid"

# Save counts
os.makedirs("data/GSE162653", exist_ok=True)
counts_df.to_csv("data/GSE162653/GSE162653_counts.txt", sep="\t")
print(f"MOF counts: {counts_df.shape} -> data/GSE162653/GSE162653_counts.txt")

# Save metadata
meta_df = pd.DataFrame({
    "sample": sample_names,
    "genotype": genotypes,
    "diet": diets,
    "condition": [f"{g}_{d}" for g, d in zip(genotypes, diets)],
})
meta_df.to_csv("data/GSE162653/sample_metadata.csv", index=False)
print(f"MOF metadata: {meta_df.shape} -> data/GSE162653/sample_metadata.csv")

# ---------------------------------------------------------------------------
# Generate gWAT KD DEG table (mimicking user's adipocyte knockdown data)
# ---------------------------------------------------------------------------
# Use same gene universe but different effects
gwat_genes = gene_names.copy()
log2fc = np.random.normal(0, 0.3, size=n_genes)

# Add strong effects for some genes (KD targets)
# Some overlap with MOF targets, some unique
log2fc[:300] += np.random.uniform(-1.5, -0.5, size=300)  # down in KD (overlap w/ MOF)
log2fc[400:600] += np.random.uniform(0.5, 1.5, size=200)  # up in KD
log2fc[800:1000] += np.random.uniform(-1.0, -0.3, size=200)  # unique to KD

# P-values (correlated with effect size)
pvalues = np.exp(-np.abs(log2fc) * np.random.uniform(3, 10, size=n_genes))
padj = np.minimum(pvalues * n_genes / np.arange(1, n_genes + 1), 1.0)
padj = np.sort(padj)  # Ensure monotonicity
np.random.shuffle(padj)  # Then shuffle back

# Recalculate properly using BH method
ranked_p = np.argsort(pvalues)
padj_bh = np.ones(n_genes)
for rank, idx in enumerate(ranked_p, 1):
    padj_bh[idx] = min(pvalues[idx] * n_genes / rank, 1.0)
# Ensure monotonicity
for i in range(len(ranked_p) - 2, -1, -1):
    idx = ranked_p[i]
    next_idx = ranked_p[i + 1]
    padj_bh[idx] = min(padj_bh[idx], padj_bh[next_idx])

gwat_deg_df = pd.DataFrame({
    "log2FoldChange": log2fc,
    "pvalue": pvalues,
    "padj": padj_bh,
    "baseMean": base_expr * np.random.uniform(0.5, 1.5, size=n_genes),
}, index=gwat_genes)
gwat_deg_df.index.name = "gene"

os.makedirs("data/gwat_kd", exist_ok=True)
gwat_deg_df.to_csv("data/gwat_kd/gwat_male_degs.csv")
print(f"gWAT DEGs: {gwat_deg_df.shape} -> data/gwat_kd/gwat_male_degs.csv")
print(f"  Significant (padj<0.05): {(gwat_deg_df['padj']<0.05).sum()}")

print("\nTest data generated successfully!")
print("Now run:")
print("  python3 scripts/02_qc_mof_data.py")
print("  python3 scripts/03_compare_mof_vs_gwat.py --mof-degs results/mof_qc/deseq2_*.csv --gwat-degs data/gwat_kd/gwat_male_degs.csv")
