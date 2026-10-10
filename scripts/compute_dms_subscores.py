#!/usr/bin/env python3
"""
compute_dms_subscores.py - Compute DMS (Differential Methylation Suitability)
A-E sub-scores for all DMR blocks in the atlas.

DMS_score = mean(A, B, C, D) * E

  A = 1 - max(per-CpG mean methylation in target samples)
  B = min(per-CpG mean methylation in background samples)
  C = max(0, 1 - 2 * mean(per-CpG std across target samples))
  D = max(0, 1 - 2 * mean(per-CpG std across background samples))
  E = mean(fraction CpGs with data in >=50% target samples,
           fraction CpGs with data in >=50% background samples)

Usage:
  python compute_dms_subscores.py \
      --atlas master_primer_assays_atlas.xlsx \
      --beta-dir /path/to/beta/files \
      --groups /path/to/groups.csv \
      --cpg-index /path/to/CpG.bed.gz \
      --output dms_subscores.csv
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
    Returns np.array of float, NaN where coverage < min_cov.
    """
    n_cpgs = end_cpg - start_cpg
    if n_cpgs < 1:
        return np.array([])
    offset = (start_cpg - 1) * 2
    try:
        with open(beta_path, "rb") as f:
            f.seek(offset)
            raw = np.frombuffer(f.read(n_cpgs * 2), dtype=np.uint8).reshape(-1, 2)
        meth = raw[:, 0].astype(float)
        cov = raw[:, 1].astype(float)
        fracs = np.where(cov >= min_cov, meth / np.maximum(cov, 1), np.nan)
        return fracs
    except Exception as e:
        print(f"  WARNING: could not read {os.path.basename(beta_path)}: {e}",
              file=sys.stderr)
        return np.full(n_cpgs, np.nan)


def cleanliness_score(tg_matrix, bg_matrix):
    """Compute DMS score and A-E sub-scores.

    tg_matrix: shape (n_target_samples, n_cpgs)
    bg_matrix: shape (n_background_samples, n_cpgs)
    """
    tg_cpg_mean = np.nanmean(tg_matrix, axis=0)
    bg_cpg_mean = np.nanmean(bg_matrix, axis=0)

    A = max(0.0, float(1.0 - np.nanmax(tg_cpg_mean))) if not np.all(np.isnan(tg_cpg_mean)) else 0.0
    B = max(0.0, float(np.nanmin(bg_cpg_mean))) if not np.all(np.isnan(bg_cpg_mean)) else 0.0
    C = max(0.0, 1.0 - 2.0 * float(np.nanmean(np.nanstd(tg_matrix, axis=0))))
    D = max(0.0, 1.0 - 2.0 * float(np.nanmean(np.nanstd(bg_matrix, axis=0))))

    tg_cov = float(np.mean(np.mean(~np.isnan(tg_matrix), axis=0) >= 0.5))
    bg_cov = float(np.mean(np.mean(~np.isnan(bg_matrix), axis=0) >= 0.5))
    E = (tg_cov + bg_cov) / 2.0

    score = np.mean([A, B, C, D]) * E
    return (round(score, 4), round(A, 3), round(B, 3),
            round(C, 3), round(D, 3), round(E, 3))


def load_cpg_index(cpg_index_path):
    """Load the wgbstools CpG index file.

    Expected format: tab-separated, columns: chr, loc (genomic position), site (1-based CpG index).
    May be gzipped (CpG.bed.gz) or plain text (CpG.bed).
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


def genomic_to_cpg_range(cpg_dict, chrom, dmr_start, dmr_end):
    """Convert genomic coordinates to CpG index range.

    Returns (start_cpg, end_cpg) where start_cpg is 1-based and end_cpg is exclusive.
    Returns None if no CpGs found in the range.
    """
    if chrom not in cpg_dict:
        return None
    positions = cpg_dict[chrom]
    sites = [s for p, s in positions if dmr_start <= p <= dmr_end]
    if not sites:
        return None
    return min(sites), max(sites) + 1


def load_groups(groups_csv, beta_dir):
    """Load groups CSV and map cell types to beta file paths.

    Returns dict: {cell_type: [beta_file_paths]}.
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
            ct_to_betas.setdefault(ct, []).append(matches[0])

    print(f"Cell types with beta files: {len(ct_to_betas)}", file=sys.stderr)
    for ct, betas in sorted(ct_to_betas.items()):
        print(f"  {ct}: {len(betas)} samples", file=sys.stderr)

    return ct_to_betas


def main():
    parser = argparse.ArgumentParser(
        description="Compute DMS A-E sub-scores for all DMR blocks in the atlas")
    parser.add_argument("--atlas", required=True,
                        help="Path to the atlas xlsx file")
    parser.add_argument("--beta-dir", required=True,
                        help="Directory containing wgbstools .beta files")
    parser.add_argument("--groups", required=True,
                        help="CSV mapping sample names to cell types")
    parser.add_argument("--cpg-index", required=True,
                        help="wgbstools CpG index file (CpG.bed or CpG.bed.gz)")
    parser.add_argument("--output", required=True,
                        help="Output CSV path")
    parser.add_argument("--sheet", default=None,
                        help="Sheet name in atlas xlsx (default: first sheet)")
    args = parser.parse_args()

    print("=" * 70)
    print("DMS Sub-score Computation", file=sys.stderr)
    print("=" * 70)

    # Load atlas
    print(f"\nLoading atlas from {args.atlas}...", file=sys.stderr)
    atlas = pd.read_excel(args.atlas, sheet_name=args.sheet)
    print(f"  {len(atlas):,} rows, {len(atlas.columns)} columns", file=sys.stderr)

    # Get unique DMR blocks
    dmr_cols = ["cell_type", "rank", "chrom", "dmr_start", "dmr_end"]
    for c in dmr_cols:
        if c not in atlas.columns:
            print(f"ERROR: column '{c}' not found in atlas", file=sys.stderr)
            sys.exit(1)
    dmr_blocks = atlas[dmr_cols].drop_duplicates().reset_index(drop=True)
    print(f"  {len(dmr_blocks):,} unique DMR blocks", file=sys.stderr)

    # Load CpG index
    cpg_dict = load_cpg_index(args.cpg_index)

    # Load groups
    ct_to_betas = load_groups(args.groups, args.beta_dir)
    all_cell_types = sorted(ct_to_betas.keys())

    # Compute sub-scores
    results = []
    n_blocks = len(dmr_blocks)

    for i, row in dmr_blocks.iterrows():
        ct = row["cell_type"]
        chrom = row["chrom"]
        dmr_start = int(row["dmr_start"])
        dmr_end = int(row["dmr_end"])

        if (i + 1) % 100 == 0:
            print(f"  Processing block {i+1}/{n_blocks}...", file=sys.stderr)

        # Convert genomic coordinates to CpG indices
        cpg_range = genomic_to_cpg_range(cpg_dict, chrom, dmr_start, dmr_end)
        if cpg_range is None:
            print(f"  WARNING: no CpGs found for {chrom}:{dmr_start}-{dmr_end}",
                  file=sys.stderr)
            continue
        start_cpg, end_cpg = cpg_range
        n_cpgs = end_cpg - start_cpg

        # Get beta files for target and background
        if ct not in ct_to_betas:
            print(f"  WARNING: no beta files for cell type '{ct}'", file=sys.stderr)
            continue
        tg_betas = ct_to_betas[ct]
        bg_betas = [b for o, betas in ct_to_betas.items() if o != ct for b in betas]
        bg_sample = bg_betas[:30]  # cap background at 30 samples for speed

        # Read beta values
        tg_mat = np.array([read_beta_region(bp, start_cpg, end_cpg) for bp in tg_betas])
        bg_mat = np.array([read_beta_region(bp, start_cpg, end_cpg) for bp in bg_sample])

        # Compute scores
        score, A, B, C, D, E = cleanliness_score(tg_mat, bg_mat)

        results.append({
            "cell_type": ct,
            "rank": int(row["rank"]),
            "chrom": chrom,
            "dmr_start": dmr_start,
            "dmr_end": dmr_end,
            "start_cpg": start_cpg,
            "end_cpg": end_cpg,
            "n_cpgs": n_cpgs,
            "dms_score": score,
            "score_A_tg_near0": A,
            "score_B_bg_near1": B,
            "score_C_tg_consist": C,
            "score_D_bg_consist": D,
            "score_E_coverage": E,
        })

    # Write output
    out_df = pd.DataFrame(results)
    out_df.to_csv(args.output, index=False)
    print(f"\nWrote {len(out_df)} rows to {args.output}", file=sys.stderr)
    print(f"  Columns: {list(out_df.columns)}", file=sys.stderr)


if __name__ == "__main__":
    main()
