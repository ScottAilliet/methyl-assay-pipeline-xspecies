#!/usr/bin/env python3
"""
Phase 3: Primer3-based primer design for bisulfite-converted DMR sequences.

Designs MSP (methylation-specific PCR) primers using Primer3 with
Primer3Plus-compatible settings. For each DMR block:

1. Fetch genomic sequence with flanking regions
2. Bisulfite convert to 6 strands
3. Design primers on the appropriate strand:
   - Methylated primers: from SM (sense methylated) or AM (antisense methylated)
   - Unmethylated primers: from SU (sense unmethylated) or AU (antisense unmethylated)
4. Apply bisulfite-specific constraints:
   - At least 2 CpGs per primer, 4 total per pair
   - Terminal CpG at 3' end for methylation discrimination
   - Tm 58-62°C (user-specified, required)
   - Product size 60-150 bp
"""

import os
import sys
import primer3
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple

from methyl_panel.config import Primer3PlusConfig, PipelineConfig
from methyl_panel.phase2_bisulfite_convert import (
    convert_sequence, fetch_genomic_sequence, BisulfiteStrands, find_cpg_positions
)


@dataclass
class PrimerPair:
    """A designed primer pair with all properties."""
    # Identification
    assay_id: str
    seq_id: str
    cell_type_id: str
    template_used: str  # SM, AM, SU, AU

    # Primer sequences
    left_primer: str
    right_primer: str
    left_primer_display: str  # With CpG notation
    right_primer_display: str

    # Positions (0-based, relative to template)
    left_start: int
    left_len: int
    right_start: int
    right_len: int

    # Tm and properties
    left_tm: float
    right_tm: float
    product_size: int

    # CpG counts
    c_total_tail: int   # Total CpGs in both primer 3' tails
    c_total: int        # Total CpGs in amplicon
    left_c_total: int   # CpGs in left primer
    right_c_total: int  # CpGs in right primer
    left_c_tail: int    # CpGs in left primer 3' tail
    right_c_tail: int   # CpGs in right primer 3' tail

    # Mismatch scores (for methylation specificity)
    sense_meth_mismatch_score: int
    sense_unmeth_mismatch_score: int
    anti_meth_mismatch_score: int
    anti_unmeth_mismatch_score: int

    # Primer3 quality metrics
    left_gc_percent: float
    right_gc_percent: float
    left_self_any: float
    right_self_any: float
    left_self_end: float
    right_self_end: float
    pair_compl_any: float
    pair_compl_end: float
    left_end_stability: float
    right_end_stability: float
    penalty: float

    # Filter results (filled by later phases)
    left_structure_mfe: Optional[float] = None
    right_structure_mfe: Optional[float] = None
    # Extended structure metrics (informational)
    left_self_dimer_dg: Optional[float] = None
    right_self_dimer_dg: Optional[float] = None
    heterodimer_dg: Optional[float] = None
    amplicon_mfe: Optional[float] = None
    primer_dimer_prediction: Optional[str] = None
    primer_dimer_end_min_dg: Optional[float] = None
    # Extended dimer metrics (informational)
    heterodimer_full_dg: Optional[float] = None
    left_homodimer_dg: Optional[float] = None
    right_homodimer_dg: Optional[float] = None
    common_variant_score: Optional[int] = None
    snp_count: Optional[int] = None
    # Repeat/segdup flags (informational)
    segdup_overlap: Optional[bool] = None
    repeat_overlap: Optional[bool] = None
    repeat_family: Optional[str] = None
    # Amplicon-level annotations (informational)
    amplicon_snp_count: Optional[int] = None
    cpg_snp_count: Optional[int] = None
    cpg_snp_details: list = None
    amplicon_segdup_overlap: Optional[bool] = None
    amplicon_repeat_overlap: Optional[bool] = None
    amplicon_repeat_family: Optional[str] = None
    mapping_error_note: Optional[str] = None
    mismatch_profile: Optional[str] = None
    mismatch_detail: Optional[str] = None

    # ePCR specificity screening (Step 4b)
    ePCR_verdict: Optional[str] = None
    ePCR_product_count: Optional[int] = None
    ePCR_products: Optional[str] = None
    ePCR_error: Optional[str] = None
    ePCR_elapsed_sec: Optional[float] = None

    # Genomic coordinates (for SNP/repeat screening)
    chrom: Optional[str] = None
    dmr_start: Optional[int] = None
    dmr_end: Optional[int] = None

    # Direction of the parent DMR (hypo or hyper)
    direction: str = "hypo"

    @classmethod
    def from_dict(cls, p: dict) -> "PrimerPair":
        """Create a PrimerPair from a primers.json dict entry.

        Handles all pipeline fields including ePCR, structure, dimer,
        SNP, repeat, and amplicon annotation results. Optional fields use
        .get() so older primers.json files without newer fields still work.
        """
        return cls(
            # Required fields (always present from Step 3)
            assay_id=p["assay_id"], seq_id=p["seq_id"],
            cell_type_id=p["cell_type_id"], template_used=p["template_used"],
            left_primer=p["left_primer"], right_primer=p["right_primer"],
            left_primer_display=p["left_primer_display"],
            right_primer_display=p["right_primer_display"],
            left_start=p["left_start"], left_len=p["left_len"],
            right_start=p["right_start"], right_len=p["right_len"],
            left_tm=p["left_tm"], right_tm=p["right_tm"],
            product_size=p["product_size"],
            c_total_tail=p["c_total_tail"], c_total=p["c_total"],
            left_c_total=p["left_c_total"], right_c_total=p["right_c_total"],
            left_c_tail=p["left_c_tail"], right_c_tail=p["right_c_tail"],
            sense_meth_mismatch_score=p["sense_meth_mismatch_score"],
            sense_unmeth_mismatch_score=p["sense_unmeth_mismatch_score"],
            anti_meth_mismatch_score=p["anti_meth_mismatch_score"],
            anti_unmeth_mismatch_score=p["anti_unmeth_mismatch_score"],
            left_gc_percent=p["left_gc_percent"], right_gc_percent=p["right_gc_percent"],
            # Primer3 self-complementarity fields (not stored in JSON, set to 0)
            left_self_any=0, right_self_any=0, left_self_end=0, right_self_end=0,
            pair_compl_any=0, pair_compl_end=0,
            left_end_stability=0, right_end_stability=0,
            penalty=p["penalty"],
            # Optional screening fields (filled by later steps, may be absent)
            left_structure_mfe=p.get("left_structure_mfe"),
            right_structure_mfe=p.get("right_structure_mfe"),
            left_self_dimer_dg=p.get("left_self_dimer_dg"),
            right_self_dimer_dg=p.get("right_self_dimer_dg"),
            heterodimer_dg=p.get("heterodimer_dg"),
            amplicon_mfe=p.get("amplicon_mfe"),
            primer_dimer_prediction=p.get("primer_dimer_prediction"),
            primer_dimer_end_min_dg=p.get("primer_dimer_end_min_dg"),
            heterodimer_full_dg=p.get("heterodimer_full_dg"),
            left_homodimer_dg=p.get("left_homodimer_dg"),
            right_homodimer_dg=p.get("right_homodimer_dg"),
            # SNP screening (Step 6)
            common_variant_score=p.get("common_variant_score"),
            snp_count=p.get("snp_count"),
            # Repeat/segdup screening (Step 6b)
            segdup_overlap=p.get("segdup_overlap"),
            repeat_overlap=p.get("repeat_overlap"),
            repeat_family=p.get("repeat_family"),
            # Amplicon-level annotations (Step 6/6b)
            amplicon_snp_count=p.get("amplicon_snp_count"),
            cpg_snp_count=p.get("cpg_snp_count"),
            amplicon_segdup_overlap=p.get("amplicon_segdup_overlap"),
            amplicon_repeat_overlap=p.get("amplicon_repeat_overlap"),
            amplicon_repeat_family=p.get("amplicon_repeat_family"),
            # ePCR screening (Step 4)
            mapping_error_note=p.get("mapping_error_note"),
            mismatch_profile=p.get("mismatch_profile"),
            mismatch_detail=p.get("mismatch_detail"),
            ePCR_verdict=p.get("ePCR_verdict"),
            ePCR_product_count=p.get("ePCR_product_count"),
            ePCR_products=p.get("ePCR_products"),
            ePCR_error=p.get("ePCR_error"),
            ePCR_elapsed_sec=p.get("ePCR_elapsed_sec"),
            # Genomic coordinates
            chrom=p.get("chrom"),
            dmr_start=p.get("dmr_start"),
            dmr_end=p.get("dmr_end"),
            direction=p.get("direction", "hypo"),
        )

    @staticmethod
    def load_primers(primers_path: str) -> List["PrimerPair"]:
        """Load primers.json and convert to PrimerPair objects."""
        import json
        with open(primers_path) as f:
            primers_data = json.load(f)
        return [PrimerPair.from_dict(p) for p in primers_data]


def count_cpgs_by_position(cpg_positions: List[int], template_len: int,
                            template_name: str,
                            start: int, length: int,
                            cpg_placement_mode: str = "strict_terminal") -> Tuple[int, int]:
    """
    Count CpG sites in a primer region using known genomic CpG positions.

    Bisulfite conversion is a substitution (C→T), not an indel, so CpG positions
    are identical between the genomic sequence and all bisulfite-converted
    templates. We use the known CpG positions from the genomic sequence rather
    than scanning the converted primer sequence, because:

    - In unmethylated templates (SU/AU), CpG C's become T's, so scanning for
      'CG' misses them. Scanning for 'TG' instead is wrong because non-CpG C's
      also become T's, creating many false 'TG' dinucleotides.
    - Position-based counting is template-independent and always correct.

    For sense strands (SM/SU), cpg_positions are top-strand 0-based positions
    and map directly to template coordinates.

    For antisense strands (AM/AU), the bottom strand is the reverse complement,
    so a top-strand CpG at position p maps to bottom-strand position len-p-2.

    Args:
        cpg_positions: 0-based positions of CpG C's on the top (sense) strand
        template_len: Length of the template sequence
        template_name: "SM", "AM", "SU", or "AU"
        start: 0-based start position of the primer in the template
        length: Length of the primer
        cpg_placement_mode: "strict_terminal" (last 5 nt) or
            "three_prime_half" (last 50% of primer length)

    Returns:
        (total_cpg, tail_cpg) where tail_cpg = CpGs in the 3' discrimination
        region (last 5 nt for strict_terminal, last 50% for three_prime_half)
    """
    # Convert top-strand CpG positions to this template's coordinate system
    if template_name in ("SM", "SU"):
        # Sense strand: positions map directly
        positions = cpg_positions
    else:
        # Antisense strand: reverse-complement mapping
        # A top-strand CpG at position p has its C at p and G at p+1.
        # On the bottom strand (reverse complement), the G becomes the
        # first base and the C becomes the second, at position len-p-2.
        positions = [template_len - p - 2 for p in cpg_positions]

    end = start + length  # exclusive

    # Count CpGs in the full primer region
    total = sum(1 for pos in positions if start <= pos < end)

    # Count CpGs in the 3' discrimination region
    if cpg_placement_mode == "three_prime_half":
        tail_size = max(1, length // 2)
    else:
        tail_size = 5
    tail_start = max(start, end - tail_size)
    tail_cpg = sum(1 for pos in positions if tail_start <= pos < end)

    return total, tail_cpg


def calculate_mismatch_score(primer_seq: str, template_strand: str,
                              primer_start: int) -> int:
    """
    Calculate mismatch score between primer and a template strand.
    Higher score = more mismatches = better discrimination.

    For MSP: a methylated primer should have many mismatches against
    the unmethylated template, and vice versa.
    """
    primer = primer_seq.upper()
    template = template_strand.upper()[primer_start:primer_start + len(primer)]

    if len(template) < len(primer):
        return len(primer)

    mismatches = sum(1 for p, t in zip(primer, template) if p != t)
    return mismatches


def design_primers_for_block(genomic_seq: str, config: Primer3PlusConfig,
                              pipeline_config: PipelineConfig,
                              assay_id: str, seq_id: str,
                              cell_type_id: str,
                              direction: str = "hypo",
                              chrom: str = None,
                              dmr_start: int = None,
                              dmr_end: int = None,
                              strands_override: "BisulfiteStrands" = None) -> List[PrimerPair]:
    """
    Design primers for a DMR block using Primer3.

    Tries both methylated (SM/AM) and unmethylated (SU/AU) templates.

    Args:
        direction: DMR direction ("hypo" or "hyper") — tagged onto each
            PrimerPair for downstream tracking.
        chrom: Chromosome name (e.g. "chr7") — propagated to PrimerPair
            for SNP/repeat screening.
        dmr_start: Genomic start position of the DMR block (1-based) —
            propagated to PrimerPair for computing primer genomic coordinates.
        dmr_end: Genomic end position of the DMR block (1-based).
    """
    # Bisulfite convert
    if strands_override is not None:
        strands = strands_override
    else:
        strands = convert_sequence(genomic_seq)

    # Validate Tm before running
    config.validate_tm()

    # Get Primer3 global args
    global_args = config.to_primer3_global_args()

    all_primers = []

    # Design on each of the 4 bisulfite-converted templates
    for template_name in ["SM", "AM", "SU", "AU"]:
        template_seq = strands.get_strand(template_name)

        # Skip if template is too short
        if len(template_seq) < 100:
            continue

        # Set up sequence args for Primer3
        seq_args = {
            "SEQUENCE_ID": f"{assay_id}_{template_name}",
            "SEQUENCE_TEMPLATE": template_seq,
        }

        # Run Primer3
        try:
            result = primer3.design_primers(seq_args, global_args)
        except Exception as e:
            continue

        n_returned = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
        if n_returned == 0:
            continue

        # Process each returned primer pair
        for i in range(min(n_returned, pipeline_config.max_primers_per_block)):
            left_seq = result.get(f"PRIMER_LEFT_{i}_SEQUENCE", "")
            right_seq = result.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")
            left_start, left_len = result.get(f"PRIMER_LEFT_{i}", (0, 0))
            right_start, right_len = result.get(f"PRIMER_RIGHT_{i}", (0, 0))

            if not left_seq or not right_seq:
                continue

            # Count CpGs using known genomic positions (template-independent)
            template_len = len(template_seq)
            cpg_mode = getattr(pipeline_config, 'cpg_placement_mode',
                               'strict_terminal')
            left_cpg, left_cpg_tail = count_cpgs_by_position(
                strands.cpg_positions, template_len, template_name,
                left_start, left_len, cpg_placement_mode=cpg_mode)
            right_cpg, right_cpg_tail = count_cpgs_by_position(
                strands.cpg_positions, template_len, template_name,
                right_start, right_len, cpg_placement_mode=cpg_mode)

            # Skip if not enough CpGs (bisulfite-specific constraint)
            if (left_cpg + right_cpg) < pipeline_config.min_cpg_pair_total:
                continue
            if left_cpg < pipeline_config.min_cpg_per_primer and right_cpg < pipeline_config.min_cpg_per_primer:
                continue

            # Calculate mismatch scores against all 4 strands
            # For a methylated primer (from SM), it should mismatch SU
            sm = strands.sense_methylated
            su = strands.sense_unmethylated
            am = strands.antisense_methylated
            au = strands.antisense_unmethylated

            # Determine which strand the primer was designed from
            # and calculate mismatch against the opposite methylation state
            if template_name in ("SM", "AM"):
                # Methylated primer: check mismatch against unmethylated
                if template_name == "SM":
                    sm_mismatch = 0  # Perfect match to SM
                    su_mismatch = calculate_mismatch_score(left_seq, su, left_start)
                    am_mismatch = calculate_mismatch_score(left_seq, am, left_start) if template_name == "SM" else 0
                    au_mismatch = calculate_mismatch_score(left_seq, au, left_start)
                else:  # AM
                    sm_mismatch = calculate_mismatch_score(left_seq, sm, left_start)
                    su_mismatch = calculate_mismatch_score(left_seq, su, left_start)
                    am_mismatch = 0
                    au_mismatch = calculate_mismatch_score(left_seq, au, left_start)
            else:  # SU or AU (unmethylated primers)
                if template_name == "SU":
                    sm_mismatch = calculate_mismatch_score(left_seq, sm, left_start)
                    su_mismatch = 0
                    am_mismatch = calculate_mismatch_score(left_seq, am, left_start)
                    au_mismatch = calculate_mismatch_score(left_seq, au, left_start)
                else:  # AU
                    sm_mismatch = calculate_mismatch_score(left_seq, sm, left_start)
                    su_mismatch = calculate_mismatch_score(left_seq, su, left_start)
                    am_mismatch = calculate_mismatch_score(left_seq, am, left_start)
                    au_mismatch = 0

            # Create display sequences with CpG notation (position-based)
            left_display = format_primer_display_by_position(
                left_seq, strands.cpg_positions, template_len,
                template_name, left_start, left_len)
            right_display = format_primer_display_by_position(
                right_seq, strands.cpg_positions, template_len,
                template_name, right_start, right_len)

            # Count total CpGs in amplicon using known genomic positions
            amp_start = left_start
            amp_end = right_start + right_len
            # Convert top-strand CpG positions to this template's coordinate system
            if template_name in ("SM", "SU"):
                amp_positions = strands.cpg_positions
            else:
                amp_positions = [template_len - p - 2 for p in strands.cpg_positions]
            amp_cpgs = sum(1 for pos in amp_positions if amp_start <= pos < amp_end)

            primer_pair = PrimerPair(
                assay_id=assay_id,
                seq_id=seq_id,
                cell_type_id=cell_type_id,
                template_used=template_name,
                left_primer=left_seq,
                right_primer=right_seq,
                left_primer_display=left_display,
                right_primer_display=right_display,
                left_start=left_start,
                left_len=left_len,
                right_start=right_start,
                right_len=right_len,
                left_tm=result.get(f"PRIMER_LEFT_{i}_TM", 0),
                right_tm=result.get(f"PRIMER_RIGHT_{i}_TM", 0),
                product_size=result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0),
                c_total_tail=left_cpg_tail + right_cpg_tail,
                c_total=amp_cpgs,
                left_c_total=left_cpg,
                right_c_total=right_cpg,
                left_c_tail=left_cpg_tail,
                right_c_tail=right_cpg_tail,
                sense_meth_mismatch_score=sm_mismatch,
                sense_unmeth_mismatch_score=su_mismatch,
                anti_meth_mismatch_score=am_mismatch,
                anti_unmeth_mismatch_score=au_mismatch,
                left_gc_percent=result.get(f"PRIMER_LEFT_{i}_GC_PERCENT", 0),
                right_gc_percent=result.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", 0),
                left_self_any=result.get(f"PRIMER_LEFT_{i}_SELF_ANY", 0),
                right_self_any=result.get(f"PRIMER_RIGHT_{i}_SELF_ANY", 0),
                left_self_end=result.get(f"PRIMER_LEFT_{i}_SELF_END", 0),
                right_self_end=result.get(f"PRIMER_RIGHT_{i}_SELF_END", 0),
                pair_compl_any=result.get(f"PRIMER_PAIR_{i}_COMPL_ANY", 0),
                pair_compl_end=result.get(f"PRIMER_PAIR_{i}_COMPL_END", 0),
                left_end_stability=result.get(f"PRIMER_LEFT_{i}_END_STABILITY", 0),
                right_end_stability=result.get(f"PRIMER_RIGHT_{i}_END_STABILITY", 0),
                penalty=result.get(f"PRIMER_PAIR_{i}_PENALTY", 0),
                direction=direction,
                chrom=chrom,
                dmr_start=dmr_start,
                dmr_end=dmr_end,
            )

            all_primers.append(primer_pair)

    return all_primers


def format_primer_display(seq: str) -> str:
    """Format primer sequence with CpG notation (lowercase y for C in CpG).

    .. deprecated:: Uses sequence-based CG scanning which misses CpG sites in
       unmethylated templates (where CpG C→T). Use format_primer_display_by_position
       instead, which uses known genomic CpG positions.
    """
    result = []
    seq_upper = seq.upper()
    for i, base in enumerate(seq_upper):
        if base == 'C' and i + 1 < len(seq_upper) and seq_upper[i + 1] == 'G':
            result.append('y')  # CpG C
        elif base == 'G' and i > 0 and seq_upper[i - 1] == 'C':
            result.append('g')  # CpG G (already shown by preceding y)
        else:
            result.append(base)
    return ''.join(result)


def format_primer_display_by_position(seq: str, cpg_positions: List[int],
                                       template_len: int, template_name: str,
                                       start: int, length: int) -> str:
    """Format primer sequence with CpG notation using known genomic positions.

    Marks CpG C's as 'y' and CpG G's as 'g' (lowercase) to indicate
    methylation-discriminatory positions. Works correctly for ALL templates
    (SM, AM, SU, AU) because it uses genomic CpG positions, not sequence scanning.

    For sense strands (SM/SU), cpg_positions map directly to template coordinates.
    For antisense strands (AM/AU), positions are converted via len-p-2.

    At each CpG position:
    - The C (position p) is marked as 'y'
    - The G (position p+1) is marked as 'g'
    In unmethylated templates, the C has been converted to T, so the display
    shows 'y' at that position even though the actual base is T — this indicates
    "this is a CpG position where methylation discrimination occurs."
    """
    # Convert top-strand CpG positions to this template's coordinate system
    if template_name in ("SM", "SU"):
        positions = cpg_positions
    else:
        positions = [template_len - p - 2 for p in cpg_positions]

    # Build a set of CpG C positions and CpG G positions within this primer
    end = start + length
    cpg_c_positions = set()  # Position of the C in CpG
    cpg_g_positions = set()  # Position of the G in CpG (C+1)

    for pos in positions:
        if start <= pos < end:
            cpg_c_positions.add(pos)
            if pos + 1 < end:
                cpg_g_positions.add(pos + 1)

    # Build display string
    result = []
    seq_upper = seq.upper()
    for i in range(length):
        abs_pos = start + i
        if abs_pos in cpg_c_positions:
            result.append('y')  # CpG C (or T in unmethylated templates)
        elif abs_pos in cpg_g_positions:
            result.append('g')  # CpG G
        else:
            result.append(seq_upper[i])
    return ''.join(result)


def main():
    """CLI entry point for testing."""
    import argparse
    parser = argparse.ArgumentParser(description="Design primers for DMR blocks")
    parser.add_argument("--genome", default="/workspace/wgbs_tools/references/hg19/hg19.fa.gz",
                        help="Path to genome FASTA")
    parser.add_argument("--region", help="Genomic region (chr:start-end)")
    parser.add_argument("--settings", help="Primer3Plus settings file")
    parser.add_argument("--min-tm", type=float, default=58.0)
    parser.add_argument("--opt-tm", type=float, default=60.0)
    parser.add_argument("--max-tm", type=float, default=62.0)
    args = parser.parse_args()

    # Load config
    if args.settings:
        config = Primer3PlusConfig.from_primer3plus_file(args.settings)
    else:
        config = Primer3PlusConfig()

    # Set Tm (required)
    config.primer_min_tm = args.min_tm
    config.primer_opt_tm = args.opt_tm
    config.primer_max_tm = args.max_tm

    pipeline_config = PipelineConfig()

    # Fetch sequence
    if args.region:
        chrom, coords = args.region.split(":")
        start, end = coords.split("-")
        seq = fetch_genomic_sequence(args.genome, chrom, int(start), int(end), flank=50)
    else:
        print("Please provide --region")
        sys.exit(1)

    print(f"Sequence length: {len(seq)}")
    print(f"Tm range: {config.primer_min_tm}-{config.primer_max_tm}°C")

    # Design primers
    primers = design_primers_for_block(
        seq, config, pipeline_config, "test_001", "TEST_0001", "TEST"
    )

    print(f"\nDesigned {len(primers)} primer pairs:")
    for p in primers[:5]:
        print(f"\n  {p.assay_id} (template: {p.template_used})")
        print(f"    Left:  {p.left_primer} (Tm={p.left_tm:.1f}°C, {p.left_c_total}CpGs)")
        print(f"    Right: {p.right_primer} (Tm={p.right_tm:.1f}°C, {p.right_c_total}CpGs)")
        print(f"    Product: {p.product_size}bp, Total CpGs: {p.c_total}")
        print(f"    Penalty: {p.penalty:.2f}")


if __name__ == "__main__":
    main()
