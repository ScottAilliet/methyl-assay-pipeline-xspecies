"""Test bisulfite conversion strand-specificity and ePCR genome array completeness.

Verifies that:
1. bisulfite_convert_bottom produces a sequence different from RC(bisulfite_convert_top),
   confirming that the bottom strand must be included separately in the ePCR genome arrays.
2. The conversion functions produce correct AM/AU strands for known sequences.
"""
import pytest
from methyl_panel.phase2_bisulfite_convert import (
    convert_sequence, bisulfite_convert_top, bisulfite_convert_bottom,
    reverse_complement, find_cpg_positions,
)


class TestStrandSpecificity:
    """Bisulfite conversion is strand-specific: RC(SM) != AM."""

    def test_rc_sm_not_equal_am(self):
        """RC of top methylated (SM) must NOT equal bottom methylated (AM).

        This is the core reason the ePCR genome arrays need both strands: if
        RC(SM) == AM, a single-strand search would find AM
        primers against the top-strand array. But since RC(SM) != AM,
        AM primers cannot be found in a top-strand-only array.
        """
        # Sequence with both CpG and non-CpG C's
        seq = "ACGTCACTGCGATTCAACGT"
        cpg_pos = find_cpg_positions(seq)

        sm = bisulfite_convert_top(seq, cpg_pos, methylated=True)
        am = bisulfite_convert_bottom(seq, cpg_pos, methylated=True)

        assert reverse_complement(sm) != am, (
            "RC(SM) should NOT equal AM -- bisulfite conversion is strand-specific"
        )

    def test_rc_su_not_equal_au(self):
        """RC of top unmethylated (SU) must NOT equal bottom unmethylated (AU)."""
        seq = "ACGTCACTGCGATTCAACGT"
        cpg_pos = find_cpg_positions(seq)

        su = bisulfite_convert_top(seq, cpg_pos, methylated=False)
        au = bisulfite_convert_bottom(seq, cpg_pos, methylated=False)

        assert reverse_complement(su) != au, (
            "RC(SU) should NOT equal AU -- bisulfite conversion is strand-specific"
        )

    def test_am_correct_for_known_sequence(self):
        """Verify AM is correctly computed for a simple sequence.

        Top:    ACGTCA  (CpG at pos 0: C,G; non-CpG C at pos 4)
        Bottom: TGACGT  (RC of top)
        Bottom CpG at pos 3 (C at 3, G at 4) -- corresponds to top CpG at pos 0
        AM:     TGACGT  (C at pos 3 is CpG -> preserved; no other C's on bottom)
        """
        seq = "ACGTCA"
        cpg_pos = find_cpg_positions(seq)  # [0]
        am = bisulfite_convert_bottom(seq, cpg_pos, methylated=True)
        assert am == "TGACGT", f"Expected TGACGT, got {am}"

    def test_au_correct_for_known_sequence(self):
        """Verify AU is correctly computed: all C->T on bottom strand.

        Top:    ACGTCA
        Bottom: TGACGT
        AU:     TGATGT  (C at pos 3 -> T)
        """
        seq = "ACGTCA"
        cpg_pos = find_cpg_positions(seq)
        au = bisulfite_convert_bottom(seq, cpg_pos, methylated=False)
        assert au == "TGATGT", f"Expected TGATGT, got {au}"

    def test_all_six_strands_present(self):
        """convert_sequence should produce all 6 strands."""
        seq = "ACGTCACTGCGATTCAACGT"
        strands = convert_sequence(seq)

        assert len(strands.sense) == len(seq)
        assert len(strands.sense_methylated) == len(seq)
        assert len(strands.sense_unmethylated) == len(seq)
        assert len(strands.antisense) == len(seq)
        assert len(strands.antisense_methylated) == len(seq)
        assert len(strands.antisense_unmethylated) == len(seq)

    def test_strand_specificity_multiple_sequences(self):
        """Verify strand specificity across multiple random-ish sequences."""
        test_seqs = [
            "ACGTCACTGCGATTCAACGT",
            "CCGATCGCGCGTATATAACG",
            "GCGCGCGCGCGCGCGCGCGC",  # all CpG
            "ATATATATATATATATATAT",  # no C at all
            "CCCCCCCCCCCCCCCCCCCC",  # all C, no CpG
        ]

        for seq in test_seqs:
            cpg_pos = find_cpg_positions(seq)
            sm = bisulfite_convert_top(seq, cpg_pos, methylated=True)
            am = bisulfite_convert_bottom(seq, cpg_pos, methylated=True)
            su = bisulfite_convert_top(seq, cpg_pos, methylated=False)
            au = bisulfite_convert_bottom(seq, cpg_pos, methylated=False)

            # RC(SM) should differ from AM unless all C's are in CpG context
            # (in which case SM == original top, AM == original bottom, and
            # RC(SM) == AM == bottom). So we only assert inequality when
            # there are non-CpG C's.
            non_cpg_c_count = sum(1 for i, b in enumerate(seq)
                                  if b == 'C' and i not in set(cpg_pos))
            if non_cpg_c_count > 0:
                assert reverse_complement(sm) != am, (
                    f"RC(SM) != AM failed for seq with non-CpG C's: {seq}"
                )
                assert reverse_complement(su) != au, (
                    f"RC(SU) != AU failed for seq with non-CpG C's: {seq}"
                )
