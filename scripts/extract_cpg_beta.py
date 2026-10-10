#!/usr/bin/env python3
"""
extract_cpg_beta.py - Extract per-CpG methylation beta values for all cell types
for every selected DMR block.

For each DMR block that has at least one passing assay, extracts beta values
(methylation fraction) for every CpG in the block, for every cell type, for
every sample. Output is a long-format CSV.

Usage:
  python extract_cpg_beta.py \
      --results results_v2.csv \
      --beta-dir /path/to/beta/files \
      --groups /path/to/groups.csv \
      --cpg-index /path/to/CpG.bed.gz \
      --output cpg_beta_values.csv
"""

import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd


MIN_COV = 5


def read_beta_region(beta_path, start_cpg, end_cpg, min_cov=MIN_COV):
    """Read per-CpG methylation fractions from a wgbstools beta file.

    Binary format: 2 x uint8 per CpG (#meth, #covered), 1-based CpG index.
    Byte offset = (start_cpg - 1) * 2.
    Returns (fracs, meth, cov) arrays: fracs is NaN where coverage < min_cov.
    """
    n_cpgs = end_cpg - start_cpg
    if n_cpgs < 1:
        return np.array([]), np.array([]), np.array([])
    offset = (start_cpg - 1) * 2
    try:
        with open(beta_path, "rb") as f:
            f.seek(offset)
            raw = np.frombuffer(f.read(n_cpgs * 2), dtype=np.uint8).reshape(-1, 2)
        meth = raw[:, 0].astype(float)
        cov = raw[:, 1].astype(float)
        fracs = np.where(cov >= min_cov, meth / np.maximum(cov, 1), np.nan)
        return fracs, meth, cov
    except Exception as e:
        print(f"  WARNING: could not read {os.path.basename(beta_path)}: {e}",
              file=sys.stderr)
        return (np.full(n_cpgs, np.nan),
                np.full(n_cpgs, np.nan),
                np.full(n_cpgs, np.nan))


def load_cpg_index(cpg_index_path):
    """Load the wgbstools CpG index file.

    Returns a dict: {chrom: sorted list of (position, cpg_index)}.
    """
    print(f"Loading CpG index from {cpg_index_path}...", file=sys.stderr)
    df = pd.read_csv(cpg_index_path, sep="\t", header=None,
                     names=["chr", "loc", "site"],
                     compression="gzip" if cpg_index_path.endswith(".gz") else None)
    cpg_dict = {}
    for chrom, group in df.groupby("chr"):
        cpg_dict[chrom] = list(zip(group["loc"].values, group["site"].values))
    total = sum(len(v) for v in cpg_dict.values())
    print(f"  Loaded {total:,} CpGs across {len(cpg_dict)} chromosomes", file=sys.stderr)
    return cpg_dict


def genomic_to_cpg_positions(cpg_dict, chrom, dmr_start, dmr_end):
    """Convert genomic coordinates to list of (cpg_index, genomic_position).

    Returns list of (site, position) tuples for CpGs within the DMR range.
    """
    if chrom not in cpg_dict:
        return []
    return [(s, p) for p, s in cpg_dict[chrom] if dmr_start <= p <= dmr_end]


def load_groups(groups_csv, beta_dir):
    """Load groups CSV and map cell types to beta file paths.

    Returns dict: {cell_type: [(beta_path, sample_name), ...]}.
    """
    groups = pd.read_csv(groups_csv)
    sample_col = [c for c in groups.columns
                  if any(k in c.lower() for k in ["sample", "file", "gsm"])][0]
    group_col = [c for c in groups.columns
                 if any(k in c.lower() for k in ["group", "type", "cell"])][0]
    print(f"Groups: sample='{sample_col}', group='{group_col}'", file=sys.stderr)

    ct_to_betas = {}
    for _, row in groups.iterrows():
        ct = str(row[group_col])
        stem = str(row[sample_col])
        matches = (glob.glob(os.path.join(beta_dir, f"*{stem}*.beta")) or
                   glob.glob(os.path.join(beta_dir, f"{stem.split('_')[0]}*.beta")))
        if matches:
            ct_to_betas.setdefault(ct, []).append((matches[0], stem))

    print(f"Cell types with beta files: {len(ct_to_betas)}", file=sys.stderr)
    for ct, betas in sorted(ct_to_betas.items()):
        print(f"  {ct}: {len(betas)} samples", file=sys.stderr)

    return ct_to_betas


def main():
    parser = argparse.ArgumentParser(
        description="Extract per-CpG beta values for all cell types for selected DMR blocks")
    parser.add_argument("--results", required=True,
                        help="ePCR screening results CSV (to identify selected DMR blocks)")
    parser.add_argument("--beta-dir", required=True,
                        help="Directory containing wgbstools .beta files")
    parser.add_argument("--groups", required=True,
                        help="CSV mapping sample names to cell types")
    parser.add_argument("--cpg-index", required=True,
                        help="wgbstools CpG index file (CpG.bed or CpG.bed.gz)")
    parser.add_argument("--output", required=True,
                        help="Output CSV path")
    args = parser.parse_args()

    print("=" * 70)
    print("Per-CpG Beta Value Extraction", file=sys.stderr)
    print("=" * 70)

    # Load screening results and identify selected DMR blocks
    print(f"\nLoading screening results from {args.results}...", file=sys.stderr)
    results = pd.read_csv(args.results)
    passing = results[results["verdict"] == "PASS"].copy()
    print(f"  {len(passing):,} passing primers", file=sys.stderr)

    # Get unique DMR blocks from passing primers
    dmr_cols = ["cell_type", "chrom", "dmr_start", "dmr_end"]
    for c in dmr_cols:
        if c not in passing.columns:
            print(f"ERROR: column '{c}' not found in results CSV", file=sys.stderr)
            sys.exit(1)
    dmr_blocks = passing[dmr_cols].drop_duplicates().reset_index(drop=True)
    print(f"  {len(dmr_blocks):,} unique DMR blocks with passing assays", file=sys.stderr)

    # Load CpG index
    cpg_dict = load_cpg_index(args.cpg_index)

    # Load groups
    ct_to_betas = load_groups(args.groups, args.beta_dir)
    all_cell_types = sorted(ct_to_betas.keys())

    # Extract beta values
    rows = []
    n_blocks = len(dmr_blocks)

    for i, row in dmr_blocks.iterrows():
        target_ct = row["cell_type"]
        chrom = row["chrom"]
        dmr_start = int(row["dmr_start"])
        dmr_end = int(row["dmr_end"])
        dmr_id = f"{target_ct}_{int(row.get('rank', 0)):04d}"

        if (i + 1) % 50 == 0:
            print(f"  Processing block {i+1}/{n_blocks}...", file=sys.stderr)

        # Get CpG positions within this DMR
        cpg_positions = genomic_to_cpg_positions(cpg_dict, chrom, dmr_start, dmr_end)
        if not cpg_positions:
            continue

        start_cpg = min(s for s, p in cpg_positions)
        end_cpg = max(s for s, p in cpg_positions) + 1

        # Extract beta values for every cell type
        for ct in all_cell_types:
            is_target = (ct == target_ct)
            for beta_path, sample_name in ct_to_betas[ct]:
                fracs, meth, cov = read_beta_region(beta_path, start_cpg, end_cpg)
                for j, (cpg_idx, cpg_pos) in enumerate(cpg_positions):
                    rows.append({
                        "dmr_id": dmr_id,
                        "cell_type_target": target_ct,
                        "chrom": chrom,
                        "dmr_start": dmr_start,
                        "dmr_end": dmr_end,
                        "cpg_index": cpg_idx,
                        "cpg_pos": cpg_pos,
                        "cell_type": ct,
                        "sample_name": sample_name,
                        "beta_value": round(float(fracs[j]), 4) if not np.isnan(fracs[j]) else "",
                        "n_meth": int(meth[j]) if not np.isnan(meth[j]) else "",
                        "n_cov": int(cov[j]) if not np.isnan(cov[j]) else "",
                        "is_target": is_target,
                    })

    # Write output
    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.output, index=False)
    print(f"\nWrote {len(out_df):,} rows to {args.output}", file=sys.stderr)
    print(f"  {n_blocks} DMR blocks x {len(all_cell_types)} cell types", file=sys.stderr)


if __name__ == "__main__":
    main()
