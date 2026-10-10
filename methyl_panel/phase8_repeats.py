#!/usr/bin/env python3
"""
Phase 8: Segmental duplication and repeat overlap screening.

Checks whether primer binding sites AND the full amplicon overlap known
segmental duplications or repetitive elements (RepeatMasker). Both are
INFORMATIONAL flags -- they do not reject primers but annotate them with
overlap status.

Primer-level screening:
  - segdup_overlap: True if either primer binding site overlaps a segdup
  - repeat_overlap: True if either primer binding site overlaps a repeat
  - repeat_family: repeat family/class string for primer overlaps

Amplicon-level screening (when amplicon coords provided):
  - amplicon_segdup_overlap: True if any part of the amplicon overlaps a segdup
  - amplicon_repeat_overlap: True if any part of the amplicon overlaps a repeat
  - amplicon_repeat_family: repeat family/class string for amplicon overlaps

Segmental duplications are relevant because primers in duplicated regions
are more likely to have off-target binding sites (though the genome-wide
ePCR screen largely catches this already).

Repetitive elements (SINE, LINE, LTR, etc.) are relevant because primers
in repeats may bind multiple genomic locations. Amplicon-level repeat
overlap is also relevant because repeats in the insert can cause PCR
artifacts or non-specific amplification.

Uses tabix to query pre-downloaded BED files from UCSC:
  - genomicSuperDups: segmental duplications track
  - rmsk: RepeatMasker track
"""

import bisect
import subprocess
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple


@dataclass
class RepeatResult:
    """Result of repeat/segdup screening for a primer pair."""
    # Primer-level fields (always populated)
    segdup_overlap: bool = False
    repeat_overlap: bool = False
    repeat_family: str = ""
    # Amplicon-level fields (populated when amplicon coords provided)
    amplicon_segdup_overlap: bool = False
    amplicon_repeat_overlap: bool = False
    amplicon_repeat_family: str = ""


def _tabix_query(bed_path: str, chrom: str, start: int, end: int) -> List[str]:
    """Query a tabix-indexed BED file for overlaps with a genomic interval.

    Args:
        bed_path: Path to tabix-indexed BED file (.gz with .tbi)
        chrom: Chromosome name (e.g. "chr7")
        start: Start position (1-based, inclusive)
        end: End position (1-based, inclusive)

    Returns:
        List of matching BED lines (tab-separated)
    """
    region = f"{chrom}:{start}-{end}"
    cmd = ["tabix", bed_path, region]
    result = subprocess.run(cmd, capture_output=True, text=True)

    lines = []
    for line in result.stdout.strip().split("\n"):
        if line and not line.startswith("#"):
            lines.append(line)
    return lines


def _extract_families(hit_lines: List[str]) -> str:
    """Extract repeat family names from BED hit lines.

    BED column 4 (0-indexed: 3) contains the name field.
    For rmsk BED, this is repFamily_repClass (e.g., "Alu_SINE").
    """
    families = set()
    for line in hit_lines:
        parts = line.split("\t")
        if len(parts) >= 4:
            families.add(parts[3])
        elif len(parts) >= 3:
            families.add(parts[2])
    return ";".join(sorted(families))


def screen_primer_pair_repeats(chrom: str,
                               left_genomic_start: int, left_len: int,
                               right_genomic_start: int, right_len: int,
                               segdup_bed: Optional[str] = None,
                               repeat_bed: Optional[str] = None,
                               amplicon_start: Optional[int] = None,
                               amplicon_end: Optional[int] = None) -> RepeatResult:
    """
    Screen a primer pair (and optionally the full amplicon) for segmental
    duplication and repeat overlap.

    Args:
        chrom: Chromosome name
        left_genomic_start: Genomic start of left primer (1-based)
        left_len: Length of left primer
        right_genomic_start: Genomic start of right primer (1-based)
        right_len: Length of right primer
        segdup_bed: Path to tabix-indexed segmental duplications BED
        repeat_bed: Path to tabix-indexed RepeatMasker BED
        amplicon_start: Genomic start of full amplicon (1-based, optional)
        amplicon_end: Genomic end of full amplicon (1-based, optional)

    Returns:
        RepeatResult with primer-level and amplicon-level overlap flags
    """
    result = RepeatResult()

    left_end = left_genomic_start + left_len - 1
    right_end = right_genomic_start + right_len - 1

    # --- Primer-level screening ---
    if segdup_bed:
        try:
            left_hits = _tabix_query(segdup_bed, chrom, left_genomic_start, left_end)
            right_hits = _tabix_query(segdup_bed, chrom, right_genomic_start, right_end)
            result.segdup_overlap = len(left_hits) > 0 or len(right_hits) > 0
        except Exception:
            pass

    if repeat_bed:
        try:
            left_hits = _tabix_query(repeat_bed, chrom, left_genomic_start, left_end)
            right_hits = _tabix_query(repeat_bed, chrom, right_genomic_start, right_end)
            all_hits = left_hits + right_hits
            result.repeat_overlap = len(all_hits) > 0
            if all_hits:
                result.repeat_family = _extract_families(all_hits)
        except Exception:
            pass

    # --- Amplicon-level screening ---
    if amplicon_start is not None and amplicon_end is not None:
        if segdup_bed:
            try:
                amp_hits = _tabix_query(segdup_bed, chrom, amplicon_start, amplicon_end)
                result.amplicon_segdup_overlap = len(amp_hits) > 0
            except Exception:
                pass

        if repeat_bed:
            try:
                amp_hits = _tabix_query(repeat_bed, chrom, amplicon_start, amplicon_end)
                result.amplicon_repeat_overlap = len(amp_hits) > 0
                if amp_hits:
                    result.amplicon_repeat_family = _extract_families(amp_hits)
            except Exception:
                pass

    return result


def _load_bed_chrom(bed_path: str, chrom: str, min_start: int, max_end: int
                    ) -> Tuple[List[Tuple[int, int, str]], List[int]]:
    """Load BED intervals for a chromosome range via a single tabix query.

    Returns (sorted_intervals, prefix_max_end) where:
      - sorted_intervals: [(start, end, name), ...] sorted by start (1-based)
      - prefix_max_end: prefix max of end positions for O(log n) overlap checks
    """
    cmd = ["tabix", bed_path, f"{chrom}:{min_start}-{max_end}"]
    result = subprocess.run(cmd, capture_output=True, text=True)

    intervals = []
    for line in result.stdout.strip().split("\n"):
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        start = int(parts[1]) + 1  # BED 0-based half-open -> 1-based inclusive
        end = int(parts[2])
        name = parts[3] if len(parts) >= 4 else ""
        intervals.append((start, end, name))

    intervals.sort(key=lambda x: x[0])
    prefix_max_end = []
    max_e = 0
    for iv in intervals:
        max_e = max(max_e, iv[1])
        prefix_max_end.append(max_e)

    return intervals, prefix_max_end


def _check_overlap(intervals: List[Tuple[int, int, str]], prefix_max_end: List[int],
                   start: int, end: int) -> Tuple[bool, List[str]]:
    """Check if any interval overlaps [start, end] using binary search + prefix max.

    Returns (has_overlap, list_of_overlapping_names).
    """
    if not intervals:
        return False, []

    starts = [iv[0] for iv in intervals]
    idx = bisect.bisect_right(starts, end)
    if idx == 0:
        return False, []

    # Quick check: if max end in [0, idx) < start, no overlap
    if prefix_max_end[idx - 1] < start:
        return False, []

    # Find all overlapping intervals
    overlapping = []
    for i in range(idx):
        iv_start, iv_end, iv_name = intervals[i]
        if iv_end >= start and iv_start <= end:
            overlapping.append(iv_name)

    return True, overlapping


def batch_screen_repeats(primers_data: List[dict],
                         segdup_bed: Optional[str] = None,
                         repeat_bed: Optional[str] = None,
                         flank: int = 100) -> List[RepeatResult]:
    """Batch screen all primers for segdup and repeat overlaps.

    Queries tabix once per chromosome per BED file instead of once per primer.
    For 10K+ primers this is ~1000x faster than per-primer screening.

    Args:
        primers_data: List of primer dicts from primers.json
        segdup_bed: Path to tabix-indexed segmental duplications BED
        repeat_bed: Path to tabix-indexed RepeatMasker BED
        flank: Flank size used during primer design

    Returns:
        List of RepeatResult objects, one per primer (same order as input)
    """
    # Group primer indices by chromosome
    by_chrom: Dict[str, List[int]] = defaultdict(list)
    for i, p in enumerate(primers_data):
        chrom = p.get("chrom")
        dmr_start = p.get("dmr_start")
        if chrom and dmr_start is not None:
            by_chrom[chrom].append(i)

    results: List[Optional[RepeatResult]] = [None] * len(primers_data)

    for chrom, indices in by_chrom.items():
        # Find min/max coordinates across all amplicons
        min_start = float('inf')
        max_end = 0
        for i in indices:
            p = primers_data[i]
            dmr_start = p["dmr_start"]
            left_genomic_start = dmr_start - flank + p["left_start"]
            right_genomic_end = dmr_start - flank + p["right_start"] + p.get("right_len", 0) - 1
            min_start = min(min_start, left_genomic_start)
            max_end = max(max_end, right_genomic_end)

        min_start = int(min_start)
        max_end = int(max_end)

        # Load BED intervals for this chromosome range (one tabix call each)
        segdup_intervals, segdup_prefix = ([], [])
        if segdup_bed:
            segdup_intervals, segdup_prefix = _load_bed_chrom(segdup_bed, chrom, min_start, max_end)

        repeat_intervals, repeat_prefix = ([], [])
        if repeat_bed:
            repeat_intervals, repeat_prefix = _load_bed_chrom(repeat_bed, chrom, min_start, max_end)

        # Screen each primer using in-memory binary search
        for i in indices:
            p = primers_data[i]
            dmr_start = p["dmr_start"]
            left_genomic_start = dmr_start - flank + p["left_start"]
            right_genomic_start = dmr_start - flank + p["right_start"]
            left_len = p["left_len"]
            right_len = p.get("right_len", len(p.get("right_primer", "")))
            left_end = left_genomic_start + left_len - 1
            right_end = right_genomic_start + right_len - 1
            amplicon_start = left_genomic_start
            amplicon_end = right_end

            result = RepeatResult()

            # Primer-level
            if segdup_bed:
                l_ovlp, _ = _check_overlap(segdup_intervals, segdup_prefix, left_genomic_start, left_end)
                r_ovlp, _ = _check_overlap(segdup_intervals, segdup_prefix, right_genomic_start, right_end)
                result.segdup_overlap = l_ovlp or r_ovlp

            if repeat_bed:
                l_ovlp, l_names = _check_overlap(repeat_intervals, repeat_prefix, left_genomic_start, left_end)
                r_ovlp, r_names = _check_overlap(repeat_intervals, repeat_prefix, right_genomic_start, right_end)
                result.repeat_overlap = l_ovlp or r_ovlp
                all_names = set(l_names + r_names)
                if all_names:
                    result.repeat_family = ";".join(sorted(all_names))

            # Amplicon-level
            if segdup_bed:
                amp_ovlp, _ = _check_overlap(segdup_intervals, segdup_prefix, amplicon_start, amplicon_end)
                result.amplicon_segdup_overlap = amp_ovlp

            if repeat_bed:
                amp_ovlp, amp_names = _check_overlap(repeat_intervals, repeat_prefix, amplicon_start, amplicon_end)
                result.amplicon_repeat_overlap = amp_ovlp
                if amp_names:
                    result.amplicon_repeat_family = ";".join(sorted(set(amp_names)))

            results[i] = result

    # Fill in defaults for primers without chrom/dmr_start
    for i in range(len(results)):
        if results[i] is None:
            results[i] = RepeatResult()

    return results


def main():
    """CLI entry point."""
    import argparse
    parser = argparse.ArgumentParser(
        description="Segmental duplication and repeat overlap screening")
    parser.add_argument("--chrom", required=True, help="Chromosome")
    parser.add_argument("--left-start", type=int, required=True,
                        help="Left primer genomic start (1-based)")
    parser.add_argument("--left-len", type=int, required=True,
                        help="Left primer length")
    parser.add_argument("--right-start", type=int, required=True,
                        help="Right primer genomic start (1-based)")
    parser.add_argument("--right-len", type=int, required=True,
                        help="Right primer length")
    parser.add_argument("--segdup-bed", help="Segmental duplications BED (tabix-indexed)")
    parser.add_argument("--repeat-bed", help="RepeatMasker BED (tabix-indexed)")
    parser.add_argument("--amplicon-start", type=int, help="Amplicon start (1-based)")
    parser.add_argument("--amplicon-end", type=int, help="Amplicon end (1-based)")
    args = parser.parse_args()

    result = screen_primer_pair_repeats(
        args.chrom, args.left_start, args.left_len,
        args.right_start, args.right_len,
        args.segdup_bed, args.repeat_bed,
        args.amplicon_start, args.amplicon_end)

    print(f"Primer segdup overlap: {result.segdup_overlap}")
    print(f"Primer repeat overlap: {result.repeat_overlap}")
    if result.repeat_family:
        print(f"Primer repeat family: {result.repeat_family}")
    print(f"Amplicon segdup overlap: {result.amplicon_segdup_overlap}")
    print(f"Amplicon repeat overlap: {result.amplicon_repeat_overlap}")
    if result.amplicon_repeat_family:
        print(f"Amplicon repeat family: {result.amplicon_repeat_family}")


if __name__ == "__main__":
    main()
