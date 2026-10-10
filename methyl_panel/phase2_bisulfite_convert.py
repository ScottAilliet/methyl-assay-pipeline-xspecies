#!/usr/bin/env python3
"""
Phase 2: Bisulfite conversion of genomic sequences.

Bisulfite conversion of DNA creates 4 new strands from the original 2:
- Original: top strand (sense) + bottom strand (antisense)
- After bisulfite conversion:
  - Sense strand, unmethylated (SU): C→T on top strand
  - Sense strand, methylated (SM): C→T except CpG C's stay C
  - Antisense strand, unmethylated (AU): C→T on bottom strand (read 5'→3')
  - Antisense strand, methylated (AM): C→T except CpG C's stay C

Plus the original unconverted strands:
  - Sense (S): original top strand
  - Antisense (A): original bottom strand (reverse complement)

Total: 6 strands (S, SM, SU, A, AM, AU)

For MSP (methylation-specific PCR):
- Methylated primers: designed from SM or AM (CpG C's preserved)
- Unmethylated primers: designed from SU or AU (all C→T)
- The primer's 3' end should be at a CpG site to discriminate methylation state
"""

from dataclasses import dataclass
from typing import List, Tuple, Optional
import subprocess
import os


@dataclass
class BisulfiteStrands:
    """The 6 bisulfite-converted strands for a genomic region."""
    sense: str           # S: original top strand
    sense_methylated: str    # SM: top strand, CpG C's preserved, other C→T
    sense_unmethylated: str  # SU: top strand, all C→T
    antisense: str       # A: original bottom strand (reverse complement)
    antisense_methylated: str    # AM: bottom strand, CpG C's preserved, other C→T
    antisense_unmethylated: str  # AU: bottom strand, all C→T
    cpg_positions: List[int]     # Positions of CpG C's on the top strand (0-based)

    def get_strand(self, name: str) -> str:
        """Get a strand by name."""
        return {
            "S": self.sense,
            "SM": self.sense_methylated,
            "SU": self.sense_unmethylated,
            "A": self.antisense,
            "AM": self.antisense_methylated,
            "AU": self.antisense_unmethylated,
        }[name]


def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'n'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))


def find_cpg_positions(seq: str) -> List[int]:
    """Find positions of C in CpG contexts on the top strand (0-based)."""
    positions = []
    seq_upper = seq.upper()
    for i in range(len(seq_upper) - 1):
        if seq_upper[i] == 'C' and seq_upper[i + 1] == 'G':
            positions.append(i)
    return positions


def bisulfite_convert_top(seq: str, cpg_positions: List[int],
                          methylated: bool) -> str:
    """
    Bisulfite convert a top (sense) strand.
    - Unmethylated: all C → T
    - Methylated: C → T except at CpG positions (C's in CpG context stay C)
    """
    seq_upper = seq.upper()
    cpg_set = set(cpg_positions)
    result = []
    for i, base in enumerate(seq_upper):
        if base == 'C':
            if methylated and i in cpg_set:
                result.append('C')  # CpG C preserved in methylated
            else:
                result.append('T')  # All other C → T
        else:
            result.append(base)
    return ''.join(result)


def bisulfite_convert_bottom(top_seq: str, cpg_positions: List[int],
                              methylated: bool) -> str:
    """
    Bisulfite convert the bottom (antisense) strand.
    The bottom strand is the reverse complement of the top strand.
    CpG on top strand = GpC on bottom strand (reading 3'→5' on bottom = 5'→3' on top).
    On the bottom strand read 5'→3', the G of CpG becomes the C to preserve.

    For unmethylated: all C → T on bottom strand
    For methylated: C → T except at positions corresponding to CpG G's on top strand
    """
    bottom = reverse_complement(top_seq)
    # CpG positions on top strand: position i has C, i+1 has G
    # On bottom strand (reverse complement), the G at top position i+1
    # becomes C at bottom position len-i-2 (0-based from 5' of bottom)
    # The C at top position i becomes G at bottom position len-i-1
    # So the CpG on bottom strand is at position len-i-2 (C) and len-i-1 (G)
    # Wait, let me think more carefully.

    # Top:    5' ...C G... 3'  (CpG at position i, i+1)
    # Bottom: 3' ...G C... 5'  (reading 3'→5')
    # Bottom read 5'→3': ...C G... (reversed)
    # The C on the bottom strand that pairs with G on top is at position:
    # bottom_5p_index = len(top) - 1 - (i + 1) = len(top) - i - 2

    bottom_cpg_positions = []
    for i in cpg_positions:
        # The G on top at position i+1 pairs with C on bottom
        # Bottom strand 5'→3' position = len - 1 - (i+1) = len - i - 2
        bottom_cpg_positions.append(len(top_seq) - i - 2)

    return bisulfite_convert_top(bottom, bottom_cpg_positions, methylated)


def convert_sequence(genomic_seq: str) -> BisulfiteStrands:
    """
    Convert a genomic sequence to all 6 bisulfite strands.

    Args:
        genomic_seq: Genomic DNA sequence (top/sense strand, 5'→3')

    Returns:
        BisulfiteStrands object with all 6 strands
    """
    seq = genomic_seq.upper()
    cpg_positions = find_cpg_positions(seq)

    # Sense strands
    sense = seq
    sense_methylated = bisulfite_convert_top(seq, cpg_positions, methylated=True)
    sense_unmethylated = bisulfite_convert_top(seq, cpg_positions, methylated=False)

    # Antisense strands
    antisense = reverse_complement(seq)
    antisense_methylated = bisulfite_convert_bottom(seq, cpg_positions, methylated=True)
    antisense_unmethylated = bisulfite_convert_bottom(seq, cpg_positions, methylated=False)

    return BisulfiteStrands(
        sense=sense,
        sense_methylated=sense_methylated,
        sense_unmethylated=sense_unmethylated,
        antisense=antisense,
        antisense_methylated=antisense_methylated,
        antisense_unmethylated=antisense_unmethylated,
        cpg_positions=cpg_positions,
    )


def fetch_genomic_sequence(genome_fasta: str, chrom: str, start: int,
                           end: int, flank: int = 0) -> str:
    """
    Fetch genomic sequence from a FASTA file using samtools faidx.

    Args:
        genome_fasta: Path to genome FASTA (can be .gz)
        chrom: Chromosome name (e.g. chr1)
        start: Start position (1-based, inclusive)
        end: End position (1-based, inclusive)
        flank: Number of bp to add on each side

    Returns:
        Genomic sequence (uppercase)
    """
    fetch_start = max(1, start - flank)
    fetch_end = end + flank
    region = f"{chrom}:{fetch_start}-{fetch_end}"

    cmd = ["samtools", "faidx", genome_fasta, region]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        if genome_fasta.endswith(".gz"):
            raise RuntimeError(
                f"samtools faidx failed on {genome_fasta}: {result.stderr.strip()}\n"
                f"This file appears to be regular gzip (not bgzip). "
                f"samtools faidx cannot seek in regular gzip files.\n"
                f"Fix: gunzip {genome_fasta} and use the uncompressed .fa file, "
                f"or recompress with bgzip: bgzip {genome_fasta[:-3]}"
            )
        raise RuntimeError(f"samtools faidx failed: {result.stderr}")

    lines = result.stdout.strip().split("\n")
    if len(lines) < 2:
        raise RuntimeError(f"No sequence returned for {region}")

    seq = "".join(lines[1:]).upper()
    return seq


def get_cpg_labels_in_amplicon(cpg_positions: List[int], amp_start: int,
                                amp_end: int) -> List[Tuple[str, int]]:
    """
    Get CpG labels (CpG_1, CpG_2, ...) for CpGs within an amplicon.

    Args:
        cpg_positions: Positions of all CpGs in the block (0-based, relative to block start)
        amp_start: Start of amplicon (0-based, relative to block start)
        amp_end: End of amplicon (0-based, relative to block start)

    Returns:
        List of (label, position) tuples for CpGs in the amplicon
    """
    labels = []
    cpg_idx = 0
    for i, pos in enumerate(cpg_positions):
        if amp_start <= pos < amp_end:
            cpg_idx += 1
            labels.append((f"CpG_{cpg_idx}", pos))
    return labels


def main():
    """CLI entry point for testing."""
    import argparse
    parser = argparse.ArgumentParser(description="Bisulfite convert a genomic sequence")
    parser.add_argument("--seq", help="Genomic sequence to convert")
    parser.add_argument("--genome", help="Path to genome FASTA")
    parser.add_argument("--region", help="Genomic region (chr:start-end)")
    args = parser.parse_args()

    if args.seq:
        seq = args.seq
    elif args.genome and args.region:
        chrom, coords = args.region.split(":")
        start, end = coords.split("-")
        seq = fetch_genomic_sequence(args.genome, chrom, int(start), int(end))
    else:
        # Demo sequence
        seq = "ATCGATCGCGATCGATCGCGATCGATCG"
        print(f"Demo sequence: {seq}")

    strands = convert_sequence(seq)
    print(f"\nCpG positions: {strands.cpg_positions}")
    print(f"S:  {strands.sense}")
    print(f"SM: {strands.sense_methylated}")
    print(f"SU: {strands.sense_unmethylated}")
    print(f"A:  {strands.antisense}")
    print(f"AM: {strands.antisense_methylated}")
    print(f"AU: {strands.antisense_unmethylated}")


if __name__ == "__main__":
    main()
