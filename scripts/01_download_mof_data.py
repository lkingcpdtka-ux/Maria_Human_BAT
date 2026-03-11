#!/usr/bin/env python3
"""
01_download_mof_data.py
-----------------------
Download and organize bulk RNA-seq data from GSE162653:
  "Histone H4 lysine 16 acetylation controls central carbon metabolism
   and diet-induced obesity in mice" (Pessoa Rodrigues et al., Nat Commun 2021)

This script downloads the GEO series matrix and supplementary count files
for the MOF heterozygous (Mof+/-) visceral WAT bulk RNA-seq dataset.

Usage:
    python3 scripts/01_download_mof_data.py [--out-dir data/GSE162653]
"""

import argparse
import gzip
import io
import os
import sys
import time

import pandas as pd
import requests


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
GEO_ACC = "GSE162653"
FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162653"
SERIES_MATRIX_URL = f"{FTP_BASE}/matrix/GSE162653_series_matrix.txt.gz"
SUPPL_BASE = f"{FTP_BASE}/suppl/"

# Direct supplementary file URLs (typical naming for featureCounts output)
SUPPL_FILES = [
    f"{FTP_BASE}/suppl/GSE162653_featureCounts_counts.txt.gz",
    f"{FTP_BASE}/suppl/GSE162653_counts.txt.gz",
    f"{FTP_BASE}/suppl/GSE162653_raw_counts.txt.gz",
]

# Nature supplementary data (Supplementary Data 3 = DEG list)
NATURE_SUPPL_BASE = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs41467-021-26277-w/MediaObjects/"
)
NATURE_SUPPL_FILES = {
    "SuppData3_DEGs.xlsx": "41467_2021_26277_MOESM6_ESM.xlsx",
    "SuppData4_ChIPseq.xlsx": "41467_2021_26277_MOESM7_ESM.xlsx",
}


def download_file(url, dest_path, retries=3, backoff=2):
    """Download a file with retries and exponential backoff."""
    for attempt in range(retries):
        try:
            print(f"  Downloading: {url}")
            resp = requests.get(url, stream=True, timeout=120)
            resp.raise_for_status()
            with open(dest_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    f.write(chunk)
            size_mb = os.path.getsize(dest_path) / 1e6
            print(f"  -> Saved to {dest_path} ({size_mb:.1f} MB)")
            return True
        except Exception as e:
            wait = backoff ** (attempt + 1)
            print(f"  Attempt {attempt + 1}/{retries} failed: {e}")
            if attempt < retries - 1:
                print(f"  Retrying in {wait}s...")
                time.sleep(wait)
    return False


def download_series_matrix(out_dir):
    """Download and parse the GEO series matrix file."""
    gz_path = os.path.join(out_dir, "GSE162653_series_matrix.txt.gz")
    txt_path = os.path.join(out_dir, "GSE162653_series_matrix.txt")

    if os.path.exists(txt_path):
        print(f"  Series matrix already exists: {txt_path}")
        return txt_path

    if download_file(SERIES_MATRIX_URL, gz_path):
        # Decompress
        with gzip.open(gz_path, "rt") as gz_in:
            with open(txt_path, "w") as txt_out:
                txt_out.write(gz_in.read())
        print(f"  Decompressed to: {txt_path}")
        return txt_path

    return None


def download_suppl_counts(out_dir):
    """Try to download supplementary count files from GEO."""
    for url in SUPPL_FILES:
        fname = url.split("/")[-1]
        dest = os.path.join(out_dir, fname)
        if os.path.exists(dest):
            print(f"  File already exists: {dest}")
            return dest
        if download_file(url, dest):
            return dest
    return None


def download_nature_suppl(out_dir):
    """Download supplementary data files from Nature article."""
    results = {}
    for local_name, remote_name in NATURE_SUPPL_FILES.items():
        dest = os.path.join(out_dir, local_name)
        if os.path.exists(dest):
            print(f"  File already exists: {dest}")
            results[local_name] = dest
            continue
        url = NATURE_SUPPL_BASE + remote_name
        if download_file(url, dest):
            results[local_name] = dest
    return results


def parse_series_matrix(matrix_path):
    """Parse sample metadata from series matrix file."""
    metadata = {}
    data_lines = []
    in_table = False

    with open(matrix_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("!Series_") or line.startswith("!series_"):
                continue
            if line.startswith("!Sample_"):
                key = line.split("\t")[0].replace("!Sample_", "")
                vals = line.split("\t")[1:]
                vals = [v.strip('"') for v in vals]
                if key not in metadata:
                    metadata[key] = []
                metadata[key].append(vals)
            if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
                in_table = True
                data_lines.append(line)
            elif in_table and not line.startswith("!"):
                data_lines.append(line)

    # Build sample info dataframe
    sample_info = {}
    if "geo_accession" in metadata:
        gsm_ids = metadata["geo_accession"][0]
        sample_info["gsm_id"] = gsm_ids
    if "title" in metadata:
        sample_info["title"] = metadata["title"][0]
    if "source_name_ch1" in metadata:
        sample_info["source"] = metadata["source_name_ch1"][0]

    # Extract characteristics
    if "characteristics_ch1" in metadata:
        for char_row in metadata["characteristics_ch1"]:
            # Each row is a list of values across samples
            if char_row and ":" in char_row[0]:
                key = char_row[0].split(":")[0].strip().lower().replace(" ", "_")
                vals = [v.split(":", 1)[1].strip() if ":" in v else v for v in char_row]
                sample_info[key] = vals

    if sample_info:
        df = pd.DataFrame(sample_info)
        return df
    return None


def list_geo_suppl_files(out_dir):
    """List available supplementary files from GEO FTP listing."""
    print("\nAttempting to list supplementary files from GEO...")
    try:
        resp = requests.get(SUPPL_BASE, timeout=30)
        if resp.status_code == 200:
            # Parse HTML directory listing for file links
            import re
            files = re.findall(r'href="([^"]+\.gz)"', resp.text)
            if files:
                print("  Available supplementary files:")
                for f in files:
                    print(f"    - {f}")
                return files
    except Exception as e:
        print(f"  Could not list files: {e}")
    return []


def main():
    parser = argparse.ArgumentParser(
        description="Download MOF het bulk RNA-seq data from GSE162653"
    )
    parser.add_argument(
        "--out-dir",
        default="data/GSE162653",
        help="Output directory for downloaded files (default: data/GSE162653)",
    )
    args = parser.parse_args()

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 70)
    print("MOF Het Bulk RNA-seq Data Download (GSE162653)")
    print("=" * 70)
    print(f"Paper: Pessoa Rodrigues et al., Nat Commun 2021")
    print(f"DOI:   10.1038/s41467-021-26277-w")
    print(f"GEO:   {GEO_ACC}")
    print(f"Output: {out_dir}")
    print()

    # Step 1: Download series matrix (sample metadata)
    print("[1/4] Downloading series matrix (sample metadata)...")
    matrix_path = download_series_matrix(out_dir)

    if matrix_path:
        print("\n[2/4] Parsing sample metadata...")
        sample_df = parse_series_matrix(matrix_path)
        if sample_df is not None:
            meta_out = os.path.join(out_dir, "sample_metadata.csv")
            sample_df.to_csv(meta_out, index=False)
            print(f"  Saved sample metadata to: {meta_out}")
            print(f"  Number of samples: {len(sample_df)}")
            print(f"  Columns: {list(sample_df.columns)}")
            print(f"\n  Sample overview:")
            print(sample_df.to_string(index=False))
        else:
            print("  Warning: Could not parse sample metadata from series matrix")
    else:
        print("  Warning: Could not download series matrix")
        print("  This may be due to network restrictions.")

    # Step 3: Try to download count data
    print("\n[3/4] Downloading supplementary count files...")
    suppl_files = list_geo_suppl_files(out_dir)
    counts_path = download_suppl_counts(out_dir)
    if counts_path:
        print(f"  Downloaded counts: {counts_path}")
    else:
        print("  Could not download count files from GEO.")

    # Step 4: Download Nature supplementary data
    print("\n[4/4] Downloading Nature supplementary data files...")
    nature_files = download_nature_suppl(out_dir)
    for name, path in nature_files.items():
        print(f"  Downloaded: {name} -> {path}")

    # Summary
    print("\n" + "=" * 70)
    print("DOWNLOAD SUMMARY")
    print("=" * 70)
    existing = [f for f in os.listdir(out_dir) if not f.startswith(".")]
    print(f"Files in {out_dir}:")
    for f in sorted(existing):
        fpath = os.path.join(out_dir, f)
        size = os.path.getsize(fpath) / 1e6
        print(f"  {f} ({size:.2f} MB)")

    if not existing or (not counts_path and not nature_files):
        print("\n" + "-" * 70)
        print("MANUAL DOWNLOAD INSTRUCTIONS")
        print("-" * 70)
        print("""
If automatic download failed (e.g., due to network restrictions), you can
manually download the data:

1. GEO Counts Matrix:
   Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162653
   Download the supplementary files and place them in: {out_dir}/

2. Nature Supplementary Data (DEG lists):
   Go to: https://www.nature.com/articles/s41467-021-26277-w
   Scroll to 'Supplementary Information'
   Download:
     - Supplementary Data 3 (DEGs) -> save as: {out_dir}/SuppData3_DEGs.xlsx
     - Supplementary Data 4 (ChIP) -> save as: {out_dir}/SuppData4_ChIPseq.xlsx

3. Or use GEO FTP:
   wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162653/suppl/
   Place .gz files in: {out_dir}/

After downloading, run: python3 scripts/02_qc_mof_data.py
""".format(out_dir=out_dir))

    print("\nDone.")


if __name__ == "__main__":
    main()
