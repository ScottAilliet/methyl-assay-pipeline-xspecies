#!/usr/bin/env python3
"""
Phase 6: Common SNP screening for primers and amplicon.

Screens primers and amplicon against dbSNP to identify common variants that could:
1. Disrupt primer binding (SNP in primer region)
2. Create false methylation calls (SNP at a CpG site -- C->T SNP reads as
   unmethylated after bisulfite conversion regardless of actual methylation)
3. Affect 3' end binding (SNP in last 5 nucleotides -- most critical for extension)

Uses a pre-downloaded dbSNP VCF file (common variants, MAF >= 1%).

Primer-level scoring: 2-digit score where:
- First digit: SNP in primer body (0=no, 1=yes)
- Second digit: SNP in 3' end last 5 nt (0=no, 1=yes)
Score 00 = clean, 10 = body SNP, 01 = 3' SNP, 11 = both

Amplicon-level screening (optional, when amplicon coords + cpg_positions provided):
- amplicon_snp_count: total common SNPs in the full amplicon (primers + insert)
- cpg_snp_count: SNPs falling at measured CpG positions (could confound methylation calls)
- cpg_snp_details: per-SNP detail with ref/alt/maf and is_c_to_t flag

Data source: dbSNP b151 GRCh37 (00-common_all.vcf.gz, MAF >= 1%).
This is the last dbSNP build with native GRCh37 coordinates (released 2018).
For newer variants, gnomAD v2 (GRCh37, 2020) can be used as a supplementary VCF
via the --dbsnp argument (any tabix-indexed VCF with CAF/TOPMED/AF INFO fields works).
"""

import bisect
import subprocess
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Optional, Set, Dict


@dataclass
class SNPResult:
    """Result of SNP screening."""
    # Primer-level fields (always populated)
    common_variant_score: int       # 2-digit score (00, 10, 01, 11)
    has_snp_in_body: bool           # SNP anywhere in primer
    has_snp_in_3prime: bool         # SNP in last 5 nucleotides
    snp_count: int                  # Total SNPs in primer binding sites
    snp_details: List[dict]         # Details of each primer SNP
    # Amplicon-level fields (populated when amplicon coords provided)
    amplicon_snp_count: int = 0     # Total SNPs in full amplicon (primers + insert)
    cpg_snp_count: int = 0          # SNPs at measured CpG positions in amplicon
    cpg_snp_details: List[dict] = field(default_factory=list)


def _query_dbsnp_region(dbsnp_vcf: str, chrom: str, start: int, end: int,
                        maf_threshold: float = 0.01) -> List[dict]:
    """Query dbSNP for common variants in a genomic interval.

    Returns list of dicts with keys: position, ref, alt, maf.
    """
    tabix_chrom = chrom.replace("chr", "") if chrom.startswith("chr") else chrom
    region = f"{tabix_chrom}:{start}-{end}"

    cmd = ["tabix", dbsnp_vcf, region]
    result = subprocess.run(cmd, capture_output=True, text=True)

    snps = []
    for line in result.stdout.strip().split("\n"):
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue

        snp_pos = int(parts[1])
        ref = parts[3]
        alt = parts[4]
        info = parts[7]

        # Extract MAF from INFO field
        # dbSNP b151 uses CAF (common allele frequencies) and TOPMED fields
        # CAF=0.5747,0.4253 means ref freq=0.5747, alt freq=0.4253
        # MAF = 1 - max(CAF values)
        maf = 0.0
        for field_str in info.split(";"):
            if field_str.startswith("CAF="):
                try:
                    freqs = [float(f) for f in field_str.split("=")[1].split(",") if f]
                    if freqs:
                        maf = max(maf, 1.0 - max(freqs), min(freqs))
                except ValueError:
                    pass
            elif field_str.startswith("TOPMED="):
                try:
                    freqs = [float(f) for f in field_str.split("=")[1].split(",") if f]
                    if freqs:
                        maf = max(maf, 1.0 - max(freqs), min(freqs))
                except ValueError:
                    pass
            elif field_str.startswith("AF=") or field_str.startswith("MAF="):
                try:
                    maf = max(maf, float(field_str.split("=")[1]))
                except ValueError:
                    pass
            elif field_str.startswith("GMAF="):
                try:
                    maf = max(maf, float(field_str.split("=")[1]))
                except ValueError:
                    pass

        if maf >= maf_threshold:
            snps.append({
                "position": snp_pos,
                "ref": ref,
                "alt": alt,
                "maf": round(maf, 4),
            })

    return snps


def screen_primer_snp(primer_seq: str, chrom: str, genomic_start: int,
                      dbsnp_vcf: str, maf_threshold: float = 0.01,
                      primer_is_reverse: bool = False) -> SNPResult:
    """
    Screen a primer against dbSNP for common variants.

    Args:
        primer_seq: Primer sequence
        chrom: Chromosome name
        genomic_start: Genomic start position of the primer (1-based)
        dbsnp_vcf: Path to dbSNP VCF file (can be .gz)
        maf_threshold: Minimum MAF to consider a SNP "common"
        primer_is_reverse: If True, primer maps to reverse strand

    Returns:
        SNPResult with primer-level score and details
    """
    primer_len = len(primer_seq)
    genomic_end = genomic_start + primer_len - 1

    raw_snps = _query_dbsnp_region(dbsnp_vcf, chrom, genomic_start, genomic_end, maf_threshold)

    snps = []
    for s in raw_snps:
        pos_in_primer = s["position"] - genomic_start
        if primer_is_reverse:
            pos_in_primer = primer_len - 1 - pos_in_primer

        in_3prime = pos_in_primer >= primer_len - 5

        snps.append({
            "position": s["position"],
            "pos_in_primer": pos_in_primer,
            "ref": s["ref"],
            "alt": s["alt"],
            "maf": s["maf"],
            "in_3prime": in_3prime,
        })

    has_body = len(snps) > 0
    has_3prime = any(s["in_3prime"] for s in snps)
    score = int(has_body) * 10 + int(has_3prime)

    return SNPResult(
        common_variant_score=score,
        has_snp_in_body=has_body,
        has_snp_in_3prime=has_3prime,
        snp_count=len(snps),
        snp_details=snps,
    )


def _screen_amplicon_snp(chrom: str, amplicon_start: int, amplicon_end: int,
                         dbsnp_vcf: str, maf_threshold: float,
                         cpg_positions: Optional[Set[int]] = None) -> tuple:
    """Screen the full amplicon for common SNPs and CpG-affecting SNPs.

    Returns (amplicon_snp_count, cpg_snp_count, cpg_snp_details).
    """
    raw_snps = _query_dbsnp_region(dbsnp_vcf, chrom, amplicon_start, amplicon_end, maf_threshold)

    cpg_set = cpg_positions or set()
    cpg_snps = []
    for s in raw_snps:
        if s["position"] in cpg_set:
            # Determine if this is a C->T SNP (most critical for methylation)
            # On forward strand: REF=C, ALT contains T
            # On reverse strand: REF=G, ALT contains A (complement of C->T)
            alt_alts = s["alt"].split(",")
            is_c_to_t = (s["ref"] == "C" and any("T" in a for a in alt_alts)) or \
                        (s["ref"] == "G" and any("A" in a for a in alt_alts))
            cpg_snps.append({
                "position": s["position"],
                "ref": s["ref"],
                "alt": s["alt"],
                "maf": s["maf"],
                "is_c_to_t": is_c_to_t,
            })

    return len(raw_snps), len(cpg_snps), cpg_snps


def screen_primer_pair_snp(left_primer: str, right_primer: str,
                           chrom: str, left_genomic_start: int,
                           right_genomic_start: int,
                           dbsnp_vcf: str,
                           maf_threshold: float = 0.01,
                           right_primer_len: Optional[int] = None,
                           cpg_positions: Optional[Set[int]] = None) -> SNPResult:
    """
    Screen both primers of a pair against dbSNP.

    When right_primer_len and cpg_positions are provided, also screens the
    full amplicon (left primer start to right primer end) for:
    - amplicon_snp_count: all common SNPs in the amplicon
    - cpg_snp_count: SNPs at measured CpG positions (confound methylation calls)

    Returns the worst-case primer score plus amplicon-level annotations.
    """
    left_result = screen_primer_snp(
        left_primer, chrom, left_genomic_start, dbsnp_vcf, maf_threshold,
        primer_is_reverse=False
    )
    right_result = screen_primer_snp(
        right_primer, chrom, right_genomic_start, dbsnp_vcf, maf_threshold,
        primer_is_reverse=True
    )

    combined_score = max(left_result.common_variant_score, right_result.common_variant_score)
    combined_snps = left_result.snp_details + right_result.snp_details

    result = SNPResult(
        common_variant_score=combined_score,
        has_snp_in_body=left_result.has_snp_in_body or right_result.has_snp_in_body,
        has_snp_in_3prime=left_result.has_snp_in_3prime or right_result.has_snp_in_3prime,
        snp_count=len(combined_snps),
        snp_details=combined_snps,
    )

    # Amplicon-level screening (optional)
    if right_primer_len is not None:
        amplicon_start = left_genomic_start
        amplicon_end = right_genomic_start + right_primer_len - 1
        amp_count, cpg_count, cpg_details = _screen_amplicon_snp(
            chrom, amplicon_start, amplicon_end, dbsnp_vcf, maf_threshold, cpg_positions
        )
        result.amplicon_snp_count = amp_count
        result.cpg_snp_count = cpg_count
        result.cpg_snp_details = cpg_details

    return result


def _find_snps_in_range(sorted_snps: List[dict], sorted_positions: List[int],
                        start: int, end: int) -> List[dict]:
    """Find SNPs within a genomic range using binary search.

    Args:
        sorted_snps: SNP dicts sorted by position
        sorted_positions: Position values (same order as sorted_snps)
        start: Start position (1-based, inclusive)
        end: End position (1-based, inclusive)

    Returns:
        List of SNP dicts within [start, end]
    """
    lo = bisect.bisect_left(sorted_positions, start)
    hi = bisect.bisect_right(sorted_positions, end)
    return sorted_snps[lo:hi]


def batch_screen_snp(primers_data: List[dict], dbsnp_vcf: str,
                     maf_threshold: float = 0.01, flank: int = 100,
                     cpg_positions_by_seqid: Optional[Dict[str, Set[int]]] = None
                     ) -> List[SNPResult]:
    """Batch screen all primers for common SNPs.

    Queries tabix once per chromosome instead of once per primer.
    For 10K+ primers this is ~1000x faster than per-primer screening.

    Args:
        primers_data: List of primer dicts from primers.json
        dbsnp_vcf: Path to tabix-indexed dbSNP VCF
        maf_threshold: Minimum MAF to consider a SNP "common"
        flank: Flank size used during primer design (for genomic coord computation)
        cpg_positions_by_seqid: Dict mapping seq_id to set of CpG positions

    Returns:
        List of SNPResult objects, one per primer (same order as input)
    """
    cpg_positions_by_seqid = cpg_positions_by_seqid or {}

    # Group primer indices by chromosome
    by_chrom: Dict[str, List[int]] = defaultdict(list)
    for i, p in enumerate(primers_data):
        chrom = p.get("chrom")
        dmr_start = p.get("dmr_start")
        if chrom and dmr_start is not None:
            by_chrom[chrom].append(i)

    results: List[Optional[SNPResult]] = [None] * len(primers_data)

    for chrom, indices in by_chrom.items():
        # Find min/max coordinates across all amplicons on this chromosome
        min_start = float('inf')
        max_end = 0
        for i in indices:
            p = primers_data[i]
            dmr_start = p["dmr_start"]
            left_genomic_start = dmr_start - flank + p["left_start"]
            right_genomic_end = dmr_start - flank + p["right_start"] + p.get("right_len", 0) - 1
            min_start = min(min_start, left_genomic_start)
            max_end = max(max_end, right_genomic_end)

        # Single tabix query for the entire chromosome range
        all_snps = _query_dbsnp_region(dbsnp_vcf, chrom, int(min_start), int(max_end), maf_threshold)
        all_snps.sort(key=lambda s: s["position"])
        snp_positions = [s["position"] for s in all_snps]

        # For each primer, find SNPs in its regions via binary search
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

            # SNPs in left primer (forward strand)
            left_raw = _find_snps_in_range(all_snps, snp_positions, left_genomic_start, left_end)
            left_details = []
            for s in left_raw:
                pos_in_primer = s["position"] - left_genomic_start
                left_details.append({
                    "position": s["position"],
                    "pos_in_primer": pos_in_primer,
                    "ref": s["ref"],
                    "alt": s["alt"],
                    "maf": s["maf"],
                    "in_3prime": pos_in_primer >= left_len - 5,
                })

            # SNPs in right primer (reverse strand)
            right_raw = _find_snps_in_range(all_snps, snp_positions, right_genomic_start, right_end)
            right_details = []
            for s in right_raw:
                pos_in_primer = right_len - 1 - (s["position"] - right_genomic_start)
                right_details.append({
                    "position": s["position"],
                    "pos_in_primer": pos_in_primer,
                    "ref": s["ref"],
                    "alt": s["alt"],
                    "maf": s["maf"],
                    "in_3prime": pos_in_primer >= right_len - 5,
                })

            combined_snps = left_details + right_details
            has_body = len(combined_snps) > 0
            has_3prime = any(s["in_3prime"] for s in combined_snps)
            score = int(has_body) * 10 + int(has_3prime)

            # Amplicon-level screening
            amp_snps = _find_snps_in_range(all_snps, snp_positions, amplicon_start, amplicon_end)

            # CpG-SNP detection
            seq_id = p.get("seq_id", "")
            cpg_positions = cpg_positions_by_seqid.get(seq_id, set())
            cpg_snps = []
            for s in amp_snps:
                if s["position"] in cpg_positions:
                    alt_alts = s["alt"].split(",")
                    is_c_to_t = (s["ref"] == "C" and any("T" in a for a in alt_alts)) or \
                                (s["ref"] == "G" and any("A" in a for a in alt_alts))
                    cpg_snps.append({
                        "position": s["position"],
                        "ref": s["ref"],
                        "alt": s["alt"],
                        "maf": s["maf"],
                        "is_c_to_t": is_c_to_t,
                    })

            results[i] = SNPResult(
                common_variant_score=score,
                has_snp_in_body=has_body,
                has_snp_in_3prime=has_3prime,
                snp_count=len(combined_snps),
                snp_details=combined_snps,
                amplicon_snp_count=len(amp_snps),
                cpg_snp_count=len(cpg_snps),
                cpg_snp_details=cpg_snps,
            )

    # Fill in defaults for primers without chrom/dmr_start
    for i in range(len(results)):
        if results[i] is None:
            results[i] = SNPResult(
                common_variant_score=0,
                has_snp_in_body=False,
                has_snp_in_3prime=False,
                snp_count=0,
                snp_details=[],
            )

    return results


def main():
    """CLI entry point."""
    import argparse
    parser = argparse.ArgumentParser(description="Common SNP screening")
    parser.add_argument("--primer", required=True, help="Primer sequence")
    parser.add_argument("--chrom", required=True, help="Chromosome")
    parser.add_argument("--start", type=int, required=True, help="Genomic start (1-based)")
    parser.add_argument("--dbsnp", required=True, help="Path to dbSNP VCF")
    parser.add_argument("--maf", type=float, default=0.01, help="MAF threshold")
    args = parser.parse_args()

    result = screen_primer_snp(args.primer, args.chrom, args.start, args.dbsnp, args.maf)

    print(f"Score: {result.common_variant_score:02d}")
    print(f"SNP in body: {result.has_snp_in_body}")
    print(f"SNP in 3' end: {result.has_snp_in_3prime}")
    print(f"Total SNPs: {result.snp_count}")
    for s in result.snp_details:
        print(f"  pos={s['position']} ref={s['ref']} alt={s['alt']} MAF={s['maf']} 3'={s['in_3prime']}")


if __name__ == "__main__":
    main()
