#!/usr/bin/env python3
"""
annotate_primers.py - Annotate passing primers with SNP, repeat, and segmental
duplication overlap data from UCSC hg19 reference tracks.

Computes values for 5 annotation columns:
  - common_variant_score: count of common SNPs (MAF >= 1%) overlapping primer binding sites
  - snp_count: total known SNPs overlapping primer binding sites
  - segdup_overlap: boolean, True if primer overlaps a segmental duplication
  - repeat_overlap: boolean, True if primer overlaps a repetitive element
  - repeat_family: name of the repetitive element family (if any)

Downloads UCSC tracks (one-time, cached locally):
  - snp151Common: Common SNPs (MAF >= 1%)
  - snp151: All known SNPs
  - rmsk: RepeatMasker repetitive elements
  - genomicSuperDups: Segmental duplications

Primer binding coordinates are found by searching the bisulfite-converted genome
within the DMR block (plus flanking sequence) for each primer sequence.

Usage:
  python annotate_primers.py \
      --atlas master_primer_assays_atlas.xlsx \
      --results results_v2.csv \
      --genome-dir hg19_bs \
      --output primer_annotations.csv \
      [--cache-dir ucsc_cache]
"""

import argparse
import gzip
import os
import struct
import sys
import urllib.request
import numpy as np
import pandas as pd

BASE_CODES = {"A": 0, "C": 1, "G": 2, "T": 3}
CODE_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}


# ---------------------------------------------------------------------------
# Genome loading (numpy .npy format, same as ePCR tool)
# ---------------------------------------------------------------------------

def load_genome(genome_dir, chrom):
    """Load a chromosome as a numpy uint8 array from .npy file."""
    path = os.path.join(genome_dir, f"{chrom}.npy")
    if not os.path.isfile(path):
        return None
    return np.load(path)


# ---------------------------------------------------------------------------
# Bisulfite conversion (string-level, for primer position finding)
# ---------------------------------------------------------------------------

def convert_ct(seq):
    """CT conversion: all C -> T. Returns converted string."""
    return seq.replace("C", "T")


def convert_ct_m(seq):
    """CT_M conversion: CpG C preserved, non-CpG C -> T. Returns converted string."""
    result = list(seq)
    for i in range(len(seq)):
        if seq[i] == "C":
            if i + 1 < len(seq) and seq[i + 1] == "G":
                continue  # CpG C preserved
            result[i] = "T"
    return "".join(result)


def convert_ga(seq):
    """GA conversion: all G -> A. Returns converted string."""
    return seq.replace("G", "A")


def convert_ga_m(seq):
    """GA_M conversion: CpG G preserved, non-CpG G -> A. Returns converted string."""
    result = list(seq)
    for i in range(len(seq)):
        if seq[i] == "G":
            if i > 0 and seq[i - 1] == "C":
                continue  # CpG G preserved
            result[i] = "A"
    return "".join(result)


def reverse_complement(seq):
    """Reverse complement of a DNA string."""
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def genome_to_string(genome_arr, start, end):
    """Extract a substring from a numpy uint8 genome array."""
    start = max(0, start)
    end = min(len(genome_arr), end)
    return "".join(CODE_TO_BASE.get(int(c), "N") for c in genome_arr[start:end])


def find_primer_position(genome_seq, primer_seq, conversions):
    """Search for a primer in converted versions of a genomic sequence.

    Trys each conversion function in order. Returns (start, end) in genomic
    coordinates (0-based, relative to genome_seq start), or None if not found.
    """
    for conv_fn in conversions:
        converted = conv_fn(genome_seq)
        pos = converted.find(primer_seq)
        if pos >= 0:
            return pos, pos + len(primer_seq)
    return None


def find_primer_binding_sites(genome_arr, chrom, dmr_start, dmr_end,
                              left_primer, right_primer, template_used):
    """Find genomic binding coordinates for both primers.

    The ePCR tool searches both FP and RC(RP) against CT/CT_M converted templates.
    For sense assays (S*): FP is in CT alphabet, RP is in GA alphabet.
    For antisense assays (A*): FP is in GA alphabet, RP is in CT alphabet.

    Strategy: search FP and RC(RP) in CT and CT_M converted forward strand.
    Also try GA and GA_M for antisense primers that might not match CT.

    Returns (fp_start, fp_end, rp_start, rp_end) in 1-based genomic coordinates,
    or None for any primer not found.
    """
    if genome_arr is None:
        return None, None, None, None

    # Search window: DMR block + 500bp flanking on each side
    flank = 500
    win_start = max(0, dmr_start - flank)
    win_end = min(len(genome_arr), dmr_end + flank)
    genome_seq = genome_to_string(genome_arr, win_start, win_end)

    # Determine primer alphabets from template
    template = (template_used or "").strip().upper()
    is_sense = template.startswith("S")
    is_methylated = template.endswith("M")

    # Forward primer search: try CT and CT_M first (works for sense),
    # then GA and GA_M (for antisense)
    fp_conversions = [convert_ct, convert_ct_m, convert_ga, convert_ga_m]
    fp_match = find_primer_position(genome_seq, left_primer.upper(), fp_conversions)

    # Reverse primer: RC(RP) searched in CT/CT_M (works for sense),
    # or RP directly in GA/GA_M (for antisense, RC of CT = GA)
    rp_upper = right_primer.upper()
    rp_rc = reverse_complement(rp_upper)
    rp_conversions = [convert_ct, convert_ct_m, convert_ga, convert_ga_m]
    rp_match = find_primer_position(genome_seq, rp_rc, rp_conversions)

    # Convert matches from window-relative to genomic coordinates (1-based)
    fp_start = fp_end = rp_start = rp_end = None
    if fp_match:
        fp_start = win_start + fp_match[0] + 1
        fp_end = win_start + fp_match[1]
    if rp_match:
        rp_start = win_start + rp_match[0] + 1
        rp_end = win_start + rp_match[1]

    return fp_start, fp_end, rp_start, rp_end


# ---------------------------------------------------------------------------
# UCSC track download and overlap checking
# ---------------------------------------------------------------------------

UCSC_TRACKS = {
    "snp151Common": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151Common.txt.gz",
        "description": "Common SNPs (MAF >= 1%)",
    },
    "snp151": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151.txt.gz",
        "description": "All known SNPs",
    },
    "rmsk": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz",
        "description": "RepeatMasker repetitive elements",
    },
    "genomicSuperDups": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz",
        "description": "Segmental duplications",
    },
}


def download_track(track_name, cache_dir):
    """Download a UCSC track to the cache directory. Returns path to cached file."""
    track_info = UCSC_TRACKS[track_name]
    cache_path = os.path.join(cache_dir, f"{track_name}.txt.gz")

    if os.path.isfile(cache_path) and os.path.getsize(cache_path) > 0:
        print(f"  {track_name}: cached ({os.path.getsize(cache_path) / 1e6:.1f} MB)",
              file=sys.stderr)
        return cache_path

    print(f"  {track_name}: downloading from UCSC...", file=sys.stderr)
    url = track_info["url"]
    urllib.request.urlretrieve(url, cache_path)
    print(f"    saved {os.path.getsize(cache_path) / 1e6:.1f} MB to {cache_path}",
          file=sys.stderr)
    return cache_path


def load_snp_track(filepath):
    """Load a UCSC SNP track (snp151 or snp151Common).

    UCSC snp151 format: tab-separated, columns include:
      chrom, chromStart, chromEnd, name, ...
    (bin, chrom, chromStart, chromEnd, name, score, strand, ...)
    We need columns 1 (chrom), 2 (start), 3 (end).

    Returns a dict: {chrom: [(start, end), ...]} sorted by start.
    """
    intervals = {}
    with gzip.open(filepath, "rt") as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            chrom = fields[1]
            start = int(fields[2])
            end = int(fields[3])
            intervals.setdefault(chrom, []).append((start, end))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def load_rmsk_track(filepath):
    """Load the RepeatMasker track.

    UCSC rmsk.txt.gz format (tab-separated):
      bin, swScore, milliDiv, milliDel, milliIns,
      genoName(chrom), genoStart, genoEnd, genoLeft, strand,
      repName, repClass, repFamily, ...

    We need genoName (5), genoStart (6), genoEnd (7), repFamily (12).

    Returns a dict: {chrom: [(start, end, rep_family), ...]} sorted by start.
    """
    intervals = {}
    with gzip.open(filepath, "rt") as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 13:
                continue
            chrom = fields[5]
            start = int(fields[6])
            end = int(fields[7])
            rep_family = fields[12]
            intervals.setdefault(chrom, []).append((start, end, rep_family))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def load_segdup_track(filepath):
    """Load the segmental duplication track.

    UCSC genomicSuperDups format: tab-separated, columns include:
      bin, chrom, chromStart, chromEnd, name, ...
    We need chrom (1), start (2), end (3).

    Returns a dict: {chrom: [(start, end), ...]} sorted by start.
    """
    intervals = {}
    with gzip.open(filepath, "rt") as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            chrom = fields[1]
            start = int(fields[2])
            end = int(fields[3])
            intervals.setdefault(chrom, []).append((start, end))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def count_overlaps(intervals_dict, chrom, qstart, qend):
    """Count how many intervals overlap [qstart, qend] using binary search.

    intervals_dict[chrom] is a sorted list of (start, end) tuples.
    """
    if chrom not in intervals_dict or qstart is None:
        return 0
    intervals = intervals_dict[chrom]
    count = 0
    for istart, iend in intervals:
        if istart > qend:
            break
        if iend >= qstart:
            count += 1
    return count


def find_repeat_family(intervals_dict, chrom, qstart, qend):
    """Find the repeat family with the largest overlap with [qstart, qend].

    Returns (has_overlap, family_name).
    """
    if chrom not in intervals_dict or qstart is None:
        return False, ""
    intervals = intervals_dict[chrom]
    best_family = ""
    best_overlap = 0
    for istart, iend, rep_family in intervals:
        if istart > qend:
            break
        overlap = min(iend, qend) - max(istart, qstart)
        if overlap > best_overlap:
            best_overlap = overlap
            best_family = rep_family
    return best_overlap > 0, best_family


def has_overlap(intervals_dict, chrom, qstart, qend):
    """Check if any interval overlaps [qstart, qend]."""
    if chrom not in intervals_dict or qstart is None:
        return False
    intervals = intervals_dict[chrom]
    for istart, iend in intervals:
        if istart > qend:
            break
        if iend >= qstart:
            return True
    return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Annotate passing primers with SNP, repeat, and segdup overlap data")
    parser.add_argument("--atlas", required=True,
                        help="Path to the atlas xlsx file")
    parser.add_argument("--results", required=True,
                        help="ePCR screening results CSV (to identify passing primers)")
    parser.add_argument("--genome-dir", required=True,
                        help="Directory containing hg19_bs chromosome .npy files")
    parser.add_argument("--output", required=True,
                        help="Output CSV path")
    parser.add_argument("--cache-dir", default="ucsc_cache",
                        help="Directory to cache downloaded UCSC tracks (default: ucsc_cache)")
    args = parser.parse_args()

    print("=" * 70)
    print("Primer Annotation: SNPs, Repeats, Segmental Duplications", file=sys.stderr)
    print("=" * 70)

    # Download UCSC tracks
    print(f"\nDownloading UCSC hg19 tracks to {args.cache_dir}/...", file=sys.stderr)
    os.makedirs(args.cache_dir, exist_ok=True)

    snp_common_path = download_track("snp151Common", args.cache_dir)
    snp_all_path = download_track("snp151", args.cache_dir)
    rmsk_path = download_track("rmsk", args.cache_dir)
    segdup_path = download_track("genomicSuperDups", args.cache_dir)

    # Load tracks into memory
    print("\nLoading tracks into memory...", file=sys.stderr)
    print("  snp151Common...", file=sys.stderr)
    snp_common = load_snp_track(snp_common_path)
    print(f"    {sum(len(v) for v in snp_common.values()):,} SNPs", file=sys.stderr)

    print("  snp151...", file=sys.stderr)
    snp_all = load_snp_track(snp_all_path)
    print(f"    {sum(len(v) for v in snp_all.values()):,} SNPs", file=sys.stderr)

    print("  rmsk...", file=sys.stderr)
    rmsk = load_rmsk_track(rmsk_path)
    print(f"    {sum(len(v) for v in rmsk.values()):,} repeat elements", file=sys.stderr)

    print("  genomicSuperDups...", file=sys.stderr)
    segdups = load_segdup_track(segdup_path)
    print(f"    {sum(len(v) for v in segdups.values()):,} segmental duplications",
          file=sys.stderr)

    # Load screening results and identify passing primers
    print(f"\nLoading screening results from {args.results}...", file=sys.stderr)
    results = pd.read_csv(args.results)
    passing = results[results["verdict"] == "PASS"].copy()
    print(f"  {len(passing):,} passing primers to annotate", file=sys.stderr)

    # Also load atlas for DMR coordinates (in case results CSV doesn't have them)
    print(f"Loading atlas from {args.atlas}...", file=sys.stderr)
    atlas = pd.read_excel(args.atlas)

    # Merge atlas columns if needed
    for col in ["chrom", "dmr_start", "dmr_end", "left_primer", "right_primer",
                "template_used", "cell_type", "rank"]:
        if col not in passing.columns and col in atlas.columns:
            key = "assay_id" if "assay_id" in passing.columns else None
            if key and key in atlas.columns:
                atlas_subset = atlas[[key, col]].drop_duplicates(key)
                passing = passing.merge(atlas_subset, on=key, how="left")

    # Check required columns
    required = ["chrom", "dmr_start", "dmr_end", "left_primer", "right_primer"]
    for col in required:
        if col not in passing.columns:
            print(f"ERROR: column '{col}' not found in results or atlas", file=sys.stderr)
            sys.exit(1)

    # Get unique chromosomes to load
    chroms_needed = sorted(passing["chrom"].unique())
    print(f"\nChromosomes to process: {chroms_needed}", file=sys.stderr)

    # Process each chromosome
    annotations = []
    genome_cache = {}

    for chrom in chroms_needed:
        genome = load_genome(args.genome_dir, chrom)
        if genome is None:
            print(f"  WARNING: genome file not found for {chrom}, skipping",
                  file=sys.stderr)
            continue
        genome_cache[chrom] = genome

        chrom_primers = passing[passing["chrom"] == chrom]
        print(f"\n  {chrom}: {len(chrom_primers)} primers", file=sys.stderr)

        for _, row in chrom_primers.iterrows():
            assay_id = row.get("assay_id", "")
            seq_id = row.get("seq_id", "")
            template = row.get("template_used", "")
            dmr_start = int(row["dmr_start"])
            dmr_end = int(row["dmr_end"])
            left_primer = str(row["left_primer"])
            right_primer = str(row["right_primer"])

            # Find primer binding sites
            fp_start, fp_end, rp_start, rp_end = find_primer_binding_sites(
                genome, chrom, dmr_start, dmr_end,
                left_primer, right_primer, template)

            # Check overlaps for each primer binding site
            # Common SNPs (MAF >= 1%)
            common_snp_fp = count_overlaps(snp_common, chrom, fp_start, fp_end)
            common_snp_rp = count_overlaps(snp_common, chrom, rp_start, rp_end)
            common_variant_score = common_snp_fp + common_snp_rp

            # All SNPs
            snp_fp = count_overlaps(snp_all, chrom, fp_start, fp_end)
            snp_rp = count_overlaps(snp_all, chrom, rp_start, rp_end)
            snp_count = snp_fp + snp_rp

            # Segmental duplications
            segdup_fp = has_overlap(segdups, chrom, fp_start, fp_end)
            segdup_rp = has_overlap(segdups, chrom, rp_start, rp_end)
            segdup_overlap = segdup_fp or segdup_rp

            # Repeats
            rep_fp, fam_fp = find_repeat_family(rmsk, chrom, fp_start, fp_end)
            rep_rp, fam_rp = find_repeat_family(rmsk, chrom, rp_start, rp_end)
            repeat_overlap = rep_fp or rep_rp

            # Pick the repeat family with the largest overlap
            if rep_fp and rep_rp:
                repeat_family = fam_fp  # both overlap, report forward
            elif rep_fp:
                repeat_family = fam_fp
            elif rep_rp:
                repeat_family = fam_rp
            else:
                repeat_family = ""

            annotations.append({
                "assay_id": assay_id,
                "seq_id": seq_id,
                "template_used": template,
                "chrom": chrom,
                "fp_start": fp_start if fp_start else "",
                "fp_end": fp_end if fp_end else "",
                "rp_start": rp_start if rp_start else "",
                "rp_end": rp_end if rp_end else "",
                "common_variant_score": common_variant_score,
                "snp_count": snp_count,
                "segdup_overlap": segdup_overlap,
                "repeat_overlap": repeat_overlap,
                "repeat_family": repeat_family,
            })

    # Write output
    out_df = pd.DataFrame(annotations)
    out_df.to_csv(args.output, index=False)
    print(f"\nWrote {len(out_df)} annotations to {args.output}", file=sys.stderr)

    # Summary
    if len(out_df) > 0:
        print(f"\nSummary:", file=sys.stderr)
        print(f"  Primers with common SNPs: "
              f"{(out_df['common_variant_score'] > 0).sum()}", file=sys.stderr)
        print(f"  Primers with any SNPs: "
              f"{(out_df['snp_count'] > 0).sum()}", file=sys.stderr)
        print(f"  Primers overlapping segdups: "
              f"{out_df['segdup_overlap'].sum()}", file=sys.stderr)
        print(f"  Primers overlapping repeats: "
              f"{out_df['repeat_overlap'].sum()}", file=sys.stderr)


if __name__ == "__main__":
    main()
