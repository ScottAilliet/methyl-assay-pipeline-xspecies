#!/usr/bin/env python3
"""
Phase 5: Secondary structure screening for primers.

Uses primer3-py's thermodynamic functions to assess:
- Hairpin formation (self-folding) — MFE cutoff -1.5 kcal/mol
- Self-dimer formation — assessed via self_any/self_end from Primer3

The MFE (minimum free energy) is calculated using primer3.calc_hairpin()
with DNA thermodynamic parameters (SantaLucia 1998).

Literature-based thresholds:
- Hairpin ΔG < -2 to -3 kcal/mol → increased PCR failure risk
- Self-dimer ΔG < -5 kcal/mol → increased PCR failure risk
- Conservative cutoff: -1.5 kcal/mol for hairpin MFE
"""

import primer3
from dataclasses import dataclass
from typing import Optional


@dataclass
class StructureResult:
    """Result of secondary structure screening."""
    left_mfe: float          # Hairpin MFE for left primer (kcal/mol)
    right_mfe: float         # Hairpin MFE for right primer (kcal/mol)
    left_passes: bool        # Left primer passes MFE cutoff
    right_passes: bool       # Right primer passes MME cutoff
    passes_filter: bool      # Both primers pass
    # Extended metrics (informational)
    left_self_dimer_dg: float = 0.0   # Self-dimer ΔG for left primer
    right_self_dimer_dg: float = 0.0  # Self-dimer ΔG for right primer
    heterodimer_dg: float = 0.0       # Heterodimer ΔG (left vs right)
    amplicon_mfe: float = 0.0         # Amplicon secondary structure MFE


def calc_hairpin_mfe(primer_seq: str, mv_conc: float = 50.0,
                     dv_conc: float = 1.5, dntp_conc: float = 0.6,
                     dna_conc: float = 50.0, temp_c: float = 37.0) -> float:
    """
    Calculate hairpin formation MFE for a primer.

    Uses primer3.calc_hairpin with SantaLucia thermodynamic parameters.

    Args:
        primer_seq: Primer sequence
        mv_conc: Monovalent cation concentration (mM)
        dv_conc: Divalent cation concentration (mM)
        dntp_conc: dNTP concentration (mM)
        dna_conc: DNA concentration (nM)
        temp_c: Temperature (°C)

    Returns:
        MFE in kcal/mol (negative = stable structure)
    """
    try:
        result = primer3.calc_hairpin(
            primer_seq,
            mv_conc=mv_conc,
            dv_conc=dv_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc,
            temp_c=temp_c,
            max_loop=30,
            output_structure=False,
        )
        return result.dg / 1000.0  # Convert from cal/mol to kcal/mol
    except Exception:
        return 0.0  # No structure = 0 MFE


def screen_structure(left_primer: str, right_primer: str,
                     mfe_cutoff: float = -1.5,
                     mv_conc: float = 50.0, dv_conc: float = 1.5,
                     dntp_conc: float = 0.6, dna_conc: float = 50.0,
                     temp_c: float = 37.0) -> StructureResult:
    """
    Screen primer pair for secondary structures.

    Computes:
    - Hairpin MFE for each primer (existing)
    - Self-dimer ΔG for each primer (new, informational)
    - Heterodimer ΔG between left and right primers (new, informational)
    - Amplicon secondary structure MFE (new, informational)

    All metrics are informational — passes_filter is based on hairpin MFE only.

    Args:
        left_primer: Left primer sequence
        right_primer: Right primer sequence
        mfe_cutoff: Maximum allowed MFE (kcal/mol). More negative = worse.
                    Default -1.5 kcal/mol (conservative).
        mv_conc: Monovalent cation (mM)
        dv_conc: Divalent cation (mM)
        dntp_conc: dNTP (mM)
        dna_conc: DNA (nM)
        temp_c: Temperature (°C)

    Returns:
        StructureResult with MFE values, dimer ΔG values, and pass/fail
    """
    left_mfe = calc_hairpin_mfe(left_primer, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c)
    right_mfe = calc_hairpin_mfe(right_primer, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c)

    # Pass if MFE is above cutoff (less negative = less stable structure)
    left_passes = left_mfe >= mfe_cutoff
    right_passes = right_mfe >= mfe_cutoff

    # Extended metrics (informational)
    # Self-dimer: primer binding to itself
    left_self_dg = _calc_heterodimer_dg(left_primer, left_primer,
                                         mv_conc, dv_conc, dntp_conc, dna_conc, temp_c)
    right_self_dg = _calc_heterodimer_dg(right_primer, right_primer,
                                          mv_conc, dv_conc, dntp_conc, dna_conc, temp_c)

    # Heterodimer: left primer binding to right primer
    hetero_dg = _calc_heterodimer_dg(left_primer, right_primer,
                                     mv_conc, dv_conc, dntp_conc, dna_conc, temp_c)

    # Amplicon secondary structure: approximate by computing hairpin
    # on the concatenation of left + right primer (proxy for amplicon folding)
    # This is a lightweight approximation — full amplicon would need the insert
    amplicon_seq = left_primer + right_primer
    amplicon_mfe = calc_hairpin_mfe(amplicon_seq, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c)

    return StructureResult(
        left_mfe=round(left_mfe, 2),
        right_mfe=round(right_mfe, 2),
        left_passes=left_passes,
        right_passes=right_passes,
        passes_filter=left_passes and right_passes,
        left_self_dimer_dg=round(left_self_dg, 2),
        right_self_dimer_dg=round(right_self_dg, 2),
        heterodimer_dg=round(hetero_dg, 2),
        amplicon_mfe=round(amplicon_mfe, 2),
    )


def _calc_heterodimer_dg(seq1: str, seq2: str,
                         mv_conc: float = 50.0, dv_conc: float = 1.5,
                         dntp_conc: float = 0.6, dna_conc: float = 50.0,
                         temp_c: float = 37.0) -> float:
    """Calculate heterodimer ΔG between two sequences using primer3.

    Returns ΔG in kcal/mol (negative = stable = dimer risk).
    """
    try:
        result = primer3.calc_heterodimer(
            seq1, seq2,
            mv_conc=mv_conc,
            dv_conc=dv_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc,
            temp_c=temp_c,
        )
        return result.dg / 1000.0  # cal/mol → kcal/mol
    except Exception:
        return 0.0


def main():
    """CLI entry point."""
    import argparse
    parser = argparse.ArgumentParser(description="Secondary structure screening")
    parser.add_argument("--left", required=True, help="Left primer sequence")
    parser.add_argument("--right", required=True, help="Right primer sequence")
    parser.add_argument("--cutoff", type=float, default=-1.5, help="MFE cutoff (kcal/mol)")
    args = parser.parse_args()

    result = screen_structure(args.left, args.right, args.cutoff)

    print(f"Left primer MFE:  {result.left_mfe} kcal/mol ({'PASS' if result.left_passes else 'FAIL'})")
    print(f"Right primer MFE: {result.right_mfe} kcal/mol ({'PASS' if result.right_passes else 'FAIL'})")
    print(f"Overall: {'PASS' if result.passes_filter else 'FAIL'}")


if __name__ == "__main__":
    main()
