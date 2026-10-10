#!/usr/bin/env python3
"""
U-assays-style PDF output module.

Generates a PDF with one page per primer pair, matching the format of
U-assays-Scott.pdf. Each page contains:
1. Assay header (assay_id, seq_id, template)
2. 6-strand sequence visualization (S, SM, SU, A, AM, AU)
3. Primer info table (sequences, Tm, positions)
4. Assay characteristics (product size, CpG counts, mismatch scores)
5. Filter results (ePCR, structure, dimer, SNP)
6. DimerDetective calculations footer
"""

import os
from typing import List, Optional
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm, cm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    PageBreak, KeepTogether
)

from methyl_panel.phase3_primer3_design import PrimerPair
from methyl_panel.phase2_bisulfite_convert import BisulfiteStrands


def format_sequence_with_cpg(seq: str, cpg_positions: list) -> str:
    """Format sequence with CpG sites highlighted in bold."""
    result = []
    seq_upper = seq.upper()
    for i, base in enumerate(seq_upper):
        if i in cpg_positions or (i > 0 and (i - 1) in cpg_positions):
            result.append(f'<b><font color="red">{base}</font></b>')
        else:
            result.append(base)
    return ''.join(result)


def generate_primer_page(p: PrimerPair, strands: Optional[BisulfiteStrands] = None) -> list:
    """Generate the content for one primer pair page."""
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle('CustomTitle', parent=styles['Title'],
                                  fontSize=14, spaceAfter=10)
    header_style = ParagraphStyle('CustomHeader', parent=styles['Heading2'],
                                   fontSize=11, spaceAfter=6, textColor=colors.HexColor('#4472C4'))
    normal_style = styles['Normal']
    mono_style = ParagraphStyle('Mono', parent=styles['Code'],
                                 fontSize=8, leading=10)

    elements = []

    # Title
    elements.append(Paragraph(f"{p.assay_id}", title_style))
    elements.append(Paragraph(
        f"seq_id: {p.seq_id} | template: {p.template_used} | cell_type: {p.cell_type_id}",
        normal_style
    ))
    elements.append(Spacer(1, 5*mm))

    # 6-strand sequence visualization
    elements.append(Paragraph("Sequence (6 strands)", header_style))
    if strands:
        strand_data = [
            ["S  (sense, unconverted):", strands.sense[:100] + ("..." if len(strands.sense) > 100 else "")],
            ["SM (sense, methylated):", strands.sense_methylated[:100] + ("..." if len(strands.sense_methylated) > 100 else "")],
            ["SU (sense, unmethylated):", strands.sense_unmethylated[:100] + ("..." if len(strands.sense_unmethylated) > 100 else "")],
            ["A  (antisense, unconverted):", strands.antisense[:100] + ("..." if len(strands.antisense) > 100 else "")],
            ["AM (antisense, methylated):", strands.antisense_methylated[:100] + ("..." if len(strands.antisense_methylated) > 100 else "")],
            ["AU (antisense, unmethylated):", strands.antisense_unmethylated[:100] + ("..." if len(strands.antisense_unmethylated) > 100 else "")],
        ]
        strand_table = Table(strand_data, colWidths=[50*mm, 120*mm])
        strand_table.setStyle(TableStyle([
            ('FONTNAME', (0, 0), (-1, -1), 'Courier'),
            ('FONTSIZE', (0, 0), (-1, -1), 7),
            ('TEXTCOLOR', (0, 0), (0, -1), colors.HexColor('#4472C4')),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
        ]))
        elements.append(strand_table)
    elements.append(Spacer(1, 5*mm))

    # Primer info table
    elements.append(Paragraph("Primer Information", header_style))
    primer_data = [
        ["", "Left Primer", "Right Primer"],
        ["Sequence (5'→3')", p.left_primer, p.right_primer],
        ["Display", p.left_primer_display, p.right_primer_display],
        ["Tm (°C)", f"{p.left_tm:.2f}", f"{p.right_tm:.2f}"],
        ["Length (bp)", str(p.left_len), str(p.right_len)],
        ["GC (%)", f"{p.left_gc_percent:.1f}", f"{p.right_gc_percent:.1f}"],
        ["CpGs (total)", str(p.left_c_total), str(p.right_c_total)],
        ["CpGs (3' tail)", str(p.left_c_tail), str(p.right_c_tail)],
    ]
    primer_table = Table(primer_data, colWidths=[40*mm, 65*mm, 65*mm])
    primer_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#4472C4')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 8),
        ('FONTNAME', (1, 1), (-1, -1), 'Courier'),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F0F0F0')]),
    ]))
    elements.append(primer_table)
    elements.append(Spacer(1, 5*mm))

    # Assay characteristics
    elements.append(Paragraph("Assay Characteristics", header_style))
    char_data = [
        ["Product size (bp)", str(p.product_size)],
        ["Total CpGs in amplicon", str(p.c_total)],
        ["Total CpGs in 3' tails", str(p.c_total_tail)],
        ["Sense meth mismatch", str(p.sense_meth_mismatch_score)],
        ["Sense unmeth mismatch", str(p.sense_unmeth_mismatch_score)],
        ["Anti meth mismatch", str(p.anti_meth_mismatch_score)],
        ["Anti unmeth mismatch", str(p.anti_unmeth_mismatch_score)],
        ["Penalty", f"{p.penalty:.2f}"],
    ]
    char_table = Table(char_data, colWidths=[60*mm, 40*mm])
    char_table.setStyle(TableStyle([
        ('FONTSIZE', (0, 0), (-1, -1), 8),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
        ('ROWBACKGROUNDS', (0, 0), (-1, -1), [colors.white, colors.HexColor('#F0F0F0')]),
    ]))
    elements.append(char_table)
    elements.append(Spacer(1, 5*mm))

    # Filter results
    elements.append(Paragraph("QC Filter Results", header_style))
    filter_data = [
        ["Filter", "Result", "Pass"],
        ["ePCR specificity", str(p.ePCR_verdict or "N/A"),
         "✓" if p.ePCR_verdict == "PASS" else ("✗" if p.ePCR_verdict and p.ePCR_verdict.startswith("FAIL") else "—")],
        ["Left structure MFE", f"{p.left_structure_mfe} kcal/mol" if p.left_structure_mfe is not None else "N/A",
         "✓" if p.left_structure_mfe and p.left_structure_mfe >= -1.5 else "—"],
        ["Right structure MFE", f"{p.right_structure_mfe} kcal/mol" if p.right_structure_mfe is not None else "N/A",
         "✓" if p.right_structure_mfe and p.right_structure_mfe >= -1.5 else "—"],
        ["Dimer prediction", str(p.primer_dimer_prediction or "N/A"),
         "✓" if p.primer_dimer_end_min_dg and p.primer_dimer_end_min_dg >= -1.0 else "—"],
        ["Dimer end_min_dg", f"{p.primer_dimer_end_min_dg} kcal/mol" if p.primer_dimer_end_min_dg is not None else "N/A", ""],
        ["Common variant score", f"{p.common_variant_score:02d}" if p.common_variant_score is not None else "N/A",
         "✓" if p.common_variant_score is not None and p.common_variant_score == 0 else "—"],
        ["Mapping note", str(p.mapping_error_note or "N/A"), ""],
    ]
    filter_table = Table(filter_data, colWidths=[50*mm, 60*mm, 20*mm])
    filter_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#4472C4')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 8),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F0F0F0')]),
        ('ALIGN', (2, 0), (2, -1), 'CENTER'),
    ]))
    elements.append(filter_table)
    elements.append(Spacer(1, 5*mm))

    # DimerDetective footer
    elements.append(Paragraph("DimerDetective Calculations", ParagraphStyle(
        'Footer', parent=styles['Heading4'], fontSize=9,
        textColor=colors.HexColor('#666666'))))
    footer_text = (
        f"Thermodynamic parameters: mv_conc=50mM, dv_conc=1.5mM, "
        f"dntp_conc=0.6mM, dna_conc=50nM, temp_c=37°C. "
        f"Both heterodimer orientations evaluated. "
        f"end_min_dg = min(left→right, right→left). "
        f"Risk tiers: high (≤-2.48), medium (-2.48 to -0.18), low (>-0.18). "
        f"Conservative cutoff: -1.0 kcal/mol. "
        f"Note: Validated on 40-60% GC primers; MSP primers (10-25% GC) "
        f"are outside validation range — interpret as indicative only."
    )
    elements.append(Paragraph(footer_text, ParagraphStyle(
        'FooterText', parent=styles['Normal'], fontSize=7,
        textColor=colors.HexColor('#666666'), leading=9)))

    return elements


def write_pdf(primer_pairs: List[PrimerPair],
              strands_list: Optional[List[BisulfiteStrands]] = None,
              output_path: str = "primers.pdf") -> str:
    """
    Write a U-assays-style PDF with one page per primer pair.

    Args:
        primer_pairs: List of PrimerPair objects
        strands_list: Optional list of BisulfiteStrands for each pair
        output_path: Output PDF path

    Returns:
        Path to the written PDF
    """
    # Write to workspace first if S3-backed
    if output_path.startswith("/mnt/"):
        tmp_path = "/workspace/_tmp_pdf.pdf"
    else:
        tmp_path = output_path

    doc = SimpleDocTemplate(
        tmp_path, pagesize=A4,
        leftMargin=15*mm, rightMargin=15*mm,
        topMargin=15*mm, bottomMargin=15*mm,
    )

    elements = []
    for i, p in enumerate(primer_pairs):
        strands = strands_list[i] if strands_list else None
        page_elements = generate_primer_page(p, strands)
        elements.extend(page_elements)
        if i < len(primer_pairs) - 1:
            elements.append(PageBreak())

    doc.build(elements)

    # Copy to S3-backed path if needed
    if output_path.startswith("/mnt/") and tmp_path != output_path:
        import shutil
        shutil.copy(tmp_path, output_path)
        os.unlink(tmp_path)

    return output_path


def main():
    """CLI entry point."""
    import argparse
    parser = argparse.ArgumentParser(description="Write U-assays-style PDF")
    parser.add_argument("--output", required=True, help="Output PDF path")
    args = parser.parse_args()

    # Create dummy data
    from methyl_panel.phase3_primer3_design import PrimerPair
    p = PrimerPair(
        assay_id="MONO_0001_M_001", seq_id="MONO_0001", cell_type_id="MONO",
        template_used="SM", left_primer="ATTTATTAATTATTAATTAT",
        right_primer="GTTGTTGTTGTTGTTGTTGT",
        left_primer_display="yTTTATTyATTATTyATTAT",
        right_primer_display="GTTGTTGTTGTTGTTGTTGT",
        left_start=0, left_len=20, right_start=80, right_len=20,
        left_tm=59.5, right_tm=60.2, product_size=100,
        c_total_tail=2, c_total=5, left_c_total=2, right_c_total=1,
        left_c_tail=1, right_c_tail=0,
        sense_meth_mismatch_score=0, sense_unmeth_mismatch_score=3,
        anti_meth_mismatch_score=2, anti_unmeth_mismatch_score=5,
        left_gc_percent=30.0, right_gc_percent=35.0,
        left_self_any=2.0, right_self_any=1.5,
        left_self_end=0.5, right_self_end=0.0,
        pair_compl_any=3.0, pair_compl_end=1.0,
        left_end_stability=-2.0, right_end_stability=-1.5,
        penalty=0.85,
        ePCR_verdict="PASS", ePCR_product_count=1, ePCR_products="chr1:1000-1100", ePCR_error="", ePCR_elapsed_sec=0.5,
        left_structure_mfe=-0.8, right_structure_mfe=-1.2,
        primer_dimer_prediction="low", primer_dimer_end_min_dg=-0.5,
        common_variant_score=0, mapping_error_note="Unique mapping",
    )

    write_pdf([p], output_path=args.output)
    print(f"Written to {args.output}")


if __name__ == "__main__":
    main()
