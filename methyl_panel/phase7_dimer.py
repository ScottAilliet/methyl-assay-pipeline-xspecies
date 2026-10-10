#!/usr/bin/env python3
"""
Phase 7: DimerDetective — primer-dimer prediction.

Thermodynamic scoring method based on nearest-neighbour 3'-end stability
calculated with primer3.calc_end_stability (SantaLucia 1998).

For each primer pair, both heterodimer orientations are evaluated:
1. 3' end of left primer annealed to right primer
2. 3' end of right primer annealed to left primer

The minimum ΔG across both orientations (end_min_dg, kcal/mol) determines
the risk classification.

Risk tiers (from DimerDetective validation on 200 qPCR assays):
- High: end_min_dg ≤ -2.48 kcal/mol (100% specificity, 76/76 dimer-forming)
- Medium: -2.48 < end_min_dg ≤ -0.18 kcal/mol (mixed zone)
- Low: end_min_dg > -0.18 kcal/mol (100% sensitivity, 0 false negatives)

Conservative cutoff: -1.0 kcal/mol (95% sensitivity, 55% specificity)

Caveat: DimerDetective was validated on standard qPCR primers with 40-60% GC.
Bisulfite-converted MSP primers have 10-25% GC, outside the validation range.
Risk labels should be interpreted as indicative, not quantitative.

Thermodynamic parameters:
- mv_conc = 50 mM, dv_conc = 1.5 mM, dntp_conc = 0.6 mM
- dna_conc = 50 nM, temp_c = 37 °C
"""

import primer3
from dataclasses import dataclass
from typing import Tuple


@dataclass
class DimerResult:
    """Result of primer-dimer prediction."""
    end_min_dg: float          # Minimum ΔG across both orientations (kcal/mol)
    left_to_right_dg: float    # ΔG for left 3' end → right primer
    right_to_left_dg: float    # ΔG for right 3' end → left primer
    risk_tier: str             # high, medium, or low
    passes_filter: bool        # Passes the conservative cutoff
    prediction: str            # Human-readable prediction
    # Extended metrics (informational)
    heterodimer_full_dg: float = 0.0   # Full heterodimer ΔG (primer3.calc_heterodimer)
    left_homodimer_dg: float = 0.0     # Left primer homodimer ΔG
    right_homodimer_dg: float = 0.0    # Right primer homodimer ΔG


def calc_end_stability_dg(seq1: str, seq2: str,
                          mv_conc: float = 50.0, dv_conc: float = 1.5,
                          dntp_conc: float = 0.6, dna_conc: float = 50.0,
                          temp_c: float = 37.0) -> float:
    """
    Calculate 3' end stability (ΔG) between two sequences.

    Uses primer3.calc_end_stability which computes the nearest-neighbour
    thermodynamic stability of the 3' end of seq1 annealed to seq2.

    Args:
        seq1: First sequence (3' end is evaluated)
        seq2: Second sequence
        mv_conc: Monovalent cation (mM)
        dv_conc: Divalent cation (mM)
        dntp_conc: dNTP (mM)
        dna_conc: DNA (nM)
        temp_c: Temperature (°C)

    Returns:
        ΔG in kcal/mol (negative = stable = dimer risk)
    """
    try:
        result = primer3.calc_end_stability(
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


def predict_dimer(left_primer: str, right_primer: str,
                  mv_conc: float = 50.0, dv_conc: float = 1.5,
                  dntp_conc: float = 0.6, dna_conc: float = 50.0,
                  temp_c: float = 37.0,
                  cutoff: float = -1.0) -> DimerResult:
    """
    Predict primer-dimer formation risk using DimerDetective method.

    Evaluates both heterodimer orientations:
    1. Left primer 3' end → right primer
    2. Right primer 3' end → left primer

    The minimum ΔG determines the risk classification.

    Args:
        left_primer: Left primer sequence (5'→3')
        right_primer: Right primer sequence (5'→3')
        mv_conc: Monovalent cation (mM) — default 50
        dv_conc: Divalent cation (mM) — default 1.5
        dntp_conc: dNTP (mM) — default 0.6
        dna_conc: DNA (nM) — default 50
        temp_c: Temperature (°C) — default 37
        cutoff: Conservative ΔG cutoff (kcal/mol) — default -1.0

    Returns:
        DimerResult with ΔG values, risk tier, and pass/fail
    """
    # Orientation 1: left 3' end → right primer
    left_to_right = calc_end_stability_dg(
        left_primer, right_primer,
        mv_conc, dv_conc, dntp_conc, dna_conc, temp_c
    )

    # Orientation 2: right 3' end → left primer
    right_to_left = calc_end_stability_dg(
        right_primer, left_primer,
        mv_conc, dv_conc, dntp_conc, dna_conc, temp_c
    )

    # Minimum ΔG (most stable = highest dimer risk)
    end_min_dg = min(left_to_right, right_to_left)

    # Risk classification
    if end_min_dg <= -2.48:
        risk_tier = "high"
        prediction = "High dimer risk (end_min_dg ≤ -2.48 kcal/mol)"
    elif end_min_dg <= -0.18:
        risk_tier = "medium"
        prediction = "Medium dimer risk (-2.48 < end_min_dg ≤ -0.18 kcal/mol)"
    else:
        risk_tier = "low"
        prediction = "Low dimer risk (end_min_dg > -0.18 kcal/mol)"

    # Pass/fail based on conservative cutoff
    passes = end_min_dg >= cutoff

    # Extended metrics (informational)
    # Full heterodimer: captures internal complementarity, not just 3' end
    heterodimer_full = _calc_full_heterodimer_dg(
        left_primer, right_primer,
        mv_conc, dv_conc, dntp_conc, dna_conc, temp_c)
    # Homodimers: each primer binding to itself
    left_homo = _calc_homodimer_dg(
        left_primer, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c)
    right_homo = _calc_homodimer_dg(
        right_primer, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c)

    return DimerResult(
        end_min_dg=round(end_min_dg, 2),
        left_to_right_dg=round(left_to_right, 2),
        right_to_left_dg=round(right_to_left, 2),
        risk_tier=risk_tier,
        passes_filter=passes,
        prediction=prediction,
        heterodimer_full_dg=round(heterodimer_full, 2),
        left_homodimer_dg=round(left_homo, 2),
        right_homodimer_dg=round(right_homo, 2),
    )


def _calc_full_heterodimer_dg(seq1: str, seq2: str,
                              mv_conc: float = 50.0, dv_conc: float = 1.5,
                              dntp_conc: float = 0.6, dna_conc: float = 50.0,
                              temp_c: float = 37.0) -> float:
    """Calculate full heterodimer ΔG using primer3.calc_heterodimer.

    This captures internal complementarity across the full primer length,
    not just the 3' end stability measured by calc_end_stability.

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


def _calc_homodimer_dg(seq: str,
                       mv_conc: float = 50.0, dv_conc: float = 1.5,
                       dntp_conc: float = 0.6, dna_conc: float = 50.0,
                       temp_c: float = 37.0) -> float:
    """Calculate homodimer ΔG using primer3.calc_homodimer.

    A homodimer is a primer binding to an identical copy of itself.

    Returns ΔG in kcal/mol (negative = stable = dimer risk).
    """
    try:
        result = primer3.calc_homodimer(
            seq,
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
    parser = argparse.ArgumentParser(description="DimerDetective primer-dimer prediction")
    parser.add_argument("--left", required=True, help="Left primer sequence (5'→3')")
    parser.add_argument("--right", required=True, help="Right primer sequence (5'→3')")
    parser.add_argument("--cutoff", type=float, default=-1.0, help="ΔG cutoff (kcal/mol)")
    parser.add_argument("--mv", type=float, default=50.0, help="Monovalent cation (mM)")
    parser.add_argument("--dv", type=float, default=1.5, help="Divalent cation (mM)")
    parser.add_argument("--dntp", type=float, default=0.6, help="dNTP (mM)")
    parser.add_argument("--dna", type=float, default=50.0, help="DNA (nM)")
    parser.add_argument("--temp", type=float, default=37.0, help="Temperature (°C)")
    args = parser.parse_args()

    result = predict_dimer(
        args.left, args.right,
        args.mv, args.dv, args.dntp, args.dna, args.temp, args.cutoff
    )

    print(f"Left→Right ΔG:  {result.left_to_right_dg} kcal/mol")
    print(f"Right→Left ΔG:  {result.right_to_left_dg} kcal/mol")
    print(f"end_min_dg:     {result.end_min_dg} kcal/mol")
    print(f"Risk tier:      {result.risk_tier}")
    print(f"Prediction:     {result.prediction}")
    print(f"Passes filter:  {result.passes_filter} (cutoff: {args.cutoff} kcal/mol)")
    print(f"\nNote: DimerDetective validated on 40-60% GC primers.")
    print(f"Bisulfite MSP primers (10-25% GC) are outside validation range.")
    print(f"Interpret results as indicative, not quantitative.")


if __name__ == "__main__":
    main()
