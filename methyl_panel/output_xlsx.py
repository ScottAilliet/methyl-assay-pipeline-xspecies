#!/usr/bin/env python3
"""
U-assays-style XLSX output module with coupled WGBS data.

Generates an Excel file with:
  Sheet 1 "Primer pairs": orderable primer sequences, CpG annotation,
    DMR metadata (gene, direction, delta_means, DMS score, p-values, etc.),
    and amplicon-level CpG summary (count, mean target/background beta).
  Sheet 2 "Amplicon CpGs": per-CpG methylation data for CpGs within
    each amplicon (position, target_mean_beta, background_mean_beta, delta).
  Sheet 3 "DMR CpGs": per-CpG methylation data for ALL CpGs in each
    DMR block (position, target_mean_beta, background_mean_beta, delta,
    in_amplicon flag).

Column groups in Sheet 1:
  Identity:     assay_id, seq_id
  DMR context:  cell_type, chrom, dmr_start, dmr_end, dmr_length_bp, gene,
                direction, dmr_rank, delta_means, delta_quants, delta_maxmin,
                ttest_pval, mwtest_pval, mvalue_ttest,
                target_mean_meth, background_mean_meth, dms_score,
                dmr_num_cpgs, annotation
  Amplicon:     amplicon_start, amplicon_end, cpgs_in_amplicon,
                amp_target_mean_beta, amp_background_mean_beta, amp_delta
  Primers:      template_used, assay, left_primer, right_primer,
                left_primer_cpg_display, right_primer_cpg_display,
                left_tm_C, right_tm_C, product_size_bp, ...
  Screening:    ePCR, structure, dimer, snp, repeat, mapping, mismatch
"""

import os
from typing import List, Optional, Dict, Any
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

from methyl_panel.phase3_primer3_design import PrimerPair
from methyl_panel.phase0_dmr_discovery import CELL_TYPE_TARGETS_BASE


# ---- Column definitions ----

# DMR-level columns (from dmr_blocks.json, looked up by seq_id)
DMR_COLUMNS = [
    "cell_type",
    "chrom",
    "dmr_start",
    "dmr_end",
    "dmr_length_bp",
    "gene",
    "direction",
    "dmr_rank",
    "delta_means",
    "delta_quants",
    "delta_maxmin",
    "ttest_pval",
    "mwtest_pval",
    "mvalue_ttest",
    "target_mean_meth",
    "background_mean_meth",
    "dms_score",
    "dmr_num_cpgs",
    "annotation",
]

# Amplicon-level columns (computed by intersecting amplicon range with cpg_sites)
AMPLICON_COLUMNS = [
    "amplicon_start",
    "amplicon_end",
    "cpgs_in_amplicon",
    "amp_target_mean_beta",
    "amp_background_mean_beta",
    "amp_delta",
    "min_non_target_delta",
    "worst_celltype",
]

# Primer + screening columns (from PrimerPair, same as before)
PRIMER_COLUMNS = [
    "template_used",
    "assay",
    "left_primer",
    "right_primer",
    "left_primer_cpg_display",
    "right_primer_cpg_display",
    "left_tm_C",
    "right_tm_C",
    "product_size_bp",
    "c_total_tail",
    "c_total",
    "left_c_total",
    "right_c_total",
    "left_c_tail",
    "right_c_tail",
    "sense_meth_mismatch_score",
    "sense_unmeth_mismatch_score",
    "anti_meth_mismatch_score",
    "anti_unmeth_mismatch_score",
    "ePCR_verdict",
    "ePCR_product_count",
    "ePCR_products",
    "ePCR_error",
    "ePCR_elapsed_sec",
    "left_structure_mfe_kcal_mol",
    "right_structure_mfe_kcal_mol",
    "left_self_dimer_dg",
    "right_self_dimer_dg",
    "heterodimer_dg",
    "amplicon_mfe",
    "primer_dimer_prediction",
    "primer_dimer_end_min_dg",
    "heterodimer_full_dg",
    "left_homodimer_dg",
    "right_homodimer_dg",
    "common_variant_score",
    "snp_count",
    "segdup_overlap",
    "repeat_overlap",
    "repeat_family",
    "amplicon_snp_count",
    "cpg_snp_count",
    "amplicon_segdup_overlap",
    "amplicon_repeat_overlap",
    "amplicon_repeat_family",
    "mapping_error_note",
    "mismatch_profile",
    "mismatch_detail",
]

# Full column order: identity + DMR context + amplicon + primers + screening
COLUMNS = ["assay_id", "seq_id"] + DMR_COLUMNS + AMPLICON_COLUMNS + PRIMER_COLUMNS

# Per-CpG sheet columns
CPG_COLUMNS = [
    "assay_id",
    "seq_id",
    "cpg_position",
    "target_mean_beta",
    "background_mean_beta",
    "delta_beta",
]

# DMR CpGs sheet: all CpGs in each DMR block (not just amplicon)
DMR_CPG_COLUMNS = [
    "assay_id",
    "seq_id",
    "cpg_label",
    "cpg_position",
    "target_mean_beta",
    "background_mean_beta",
    "delta_beta",
    "in_amplicon",
]


def _build_dmr_lookup(dmr_blocks: Optional[List[dict]]) -> Dict[str, dict]:
    """Build a {seq_id: block_dict} lookup from dmr_blocks list."""
    if not dmr_blocks:
        return {}
    return {b["seq_id"]: b for b in dmr_blocks if "seq_id" in b}


def _compute_amplicon_cpgs(p: PrimerPair, block: Optional[dict]) -> dict:
    """
    Compute amplicon genomic range and CpG statistics.

    Returns a dict with amplicon_start, amplicon_end, cpgs_in_amplicon,
    amp_target_mean_beta, amp_background_mean_beta, amp_delta,
    and _cpg_rows (list of per-CpG dicts for the second sheet).
    """
    result = {
        "amplicon_start": None,
        "amplicon_end": None,
        "cpgs_in_amplicon": None,
        "amp_target_mean_beta": None,
        "amp_background_mean_beta": None,
        "amp_delta": None,
        "_cpg_rows": [],
    }

    if not block or p.dmr_start is None or p.left_start is None:
        return result

    # Amplicon genomic coordinates
    amp_start = p.dmr_start + p.left_start
    amp_end = p.dmr_start + p.right_start + p.right_len
    result["amplicon_start"] = amp_start
    result["amplicon_end"] = amp_end

    cpg_sites = block.get("cpg_sites", [])
    if not cpg_sites:
        return result

    # Filter CpGs within amplicon range
    amp_cpgs = [c for c in cpg_sites if amp_start <= c["position"] <= amp_end]
    result["cpgs_in_amplicon"] = len(amp_cpgs)

    if amp_cpgs:
        target_betas = [c["target_mean_beta"] for c in amp_cpgs]
        bg_betas = [c["background_mean_beta"] for c in amp_cpgs]
        amp_target = sum(target_betas) / len(target_betas)
        amp_bg = sum(bg_betas) / len(bg_betas)
        result["amp_target_mean_beta"] = round(amp_target, 4)
        result["amp_background_mean_beta"] = round(amp_bg, 4)
        result["amp_delta"] = round(abs(amp_target - amp_bg), 4)

        # Build per-CpG rows for the second sheet
        for c in amp_cpgs:
            delta = c.get("delta_beta")
            if delta is None:
                delta = round(c["target_mean_beta"] - c["background_mean_beta"], 4)
            result["_cpg_rows"].append({
                "assay_id": p.assay_id,
                "seq_id": p.seq_id,
                "cpg_position": c["position"],
                "target_mean_beta": c["target_mean_beta"],
                "background_mean_beta": c["background_mean_beta"],
                "delta_beta": delta,
            })

    return result


def _enrich_row_with_dmr(row: dict, p: PrimerPair, block: Optional[dict]) -> tuple:
    """
    Add DMR-level and amplicon-level fields to a primer row.
    Returns (amplicon_cpg_rows, dmr_cpg_rows, celltype_cpg_rows) for
    sheets 2, 3, and 4.
    """
    # DMR-level fields
    if block:
        row["gene"] = block.get("gene")
        row["dmr_rank"] = block.get("rank")
        row["delta_means"] = block.get("delta_means")
        row["delta_quants"] = block.get("delta_quants")
        row["delta_maxmin"] = block.get("delta_maxmin")
        row["ttest_pval"] = block.get("ttest_pval")
        row["mwtest_pval"] = block.get("mwtest_pval")
        row["mvalue_ttest"] = block.get("mvalue_ttest")
        row["target_mean_meth"] = block.get("target_mean_meth")
        row["background_mean_meth"] = block.get("background_mean_meth")
        row["dms_score"] = block.get("dms_score", block.get("cleanliness_score"))
        row["dmr_num_cpgs"] = block.get("num_cpgs")
        row["annotation"] = block.get("annotation")

    # Fields already on PrimerPair but not previously in Excel
    row["cell_type"] = p.cell_type_id
    row["chrom"] = p.chrom
    row["dmr_start"] = p.dmr_start
    row["dmr_end"] = p.dmr_end
    row["dmr_length_bp"] = (p.dmr_end - p.dmr_start) if (p.dmr_start and p.dmr_end) else None
    row["direction"] = p.direction

    # Amplicon-level fields
    amp = _compute_amplicon_cpgs(p, block)
    row["amplicon_start"] = amp["amplicon_start"]
    row["amplicon_end"] = amp["amplicon_end"]
    row["cpgs_in_amplicon"] = amp["cpgs_in_amplicon"]
    row["amp_target_mean_beta"] = amp["amp_target_mean_beta"]
    row["amp_background_mean_beta"] = amp["amp_background_mean_beta"]
    row["amp_delta"] = amp["amp_delta"]

    # Specificity metric: worst-case separation from non-target cell types
    # For each amplicon CpG, per_celltype_beta has {cell_type: mean_beta}.
    # Compute amplicon-level mean per cell type, then find the non-target
    # cell type with the smallest delta from the target.
    row["min_non_target_delta"] = None
    row["worst_celltype"] = None
    if block and p.dmr_start is not None and amp["amplicon_start"] is not None:
        amp_start = amp["amplicon_start"]
        amp_end = amp["amplicon_end"]
        target_ct = p.cell_type_id
        # Expand target to atlas group names (handles merged blood types
        # like CD4T -> [Blood-T-CD4, Blood-T-CenMem-CD4, ...])
        target_groups = set(CELL_TYPE_TARGETS_BASE.get(target_ct, [target_ct]))
        amp_cpgs = [c for c in block.get("cpg_sites", [])
                    if amp_start <= c["position"] <= amp_end
                    and c.get("per_celltype_beta")]
        if amp_cpgs:
            # Collect all cell type names across amplicon CpGs
            all_cts = set()
            for c in amp_cpgs:
                all_cts.update(c["per_celltype_beta"].keys())

            # Amplicon-level mean beta per cell type
            ct_amp_means = {}
            for ct in all_cts:
                vals = [c["per_celltype_beta"][ct] for c in amp_cpgs
                        if ct in c["per_celltype_beta"]]
                if vals:
                    ct_amp_means[ct] = sum(vals) / len(vals)

            # Target amplicon mean (from existing field or from per_celltype_beta)
            target_amp = amp["amp_target_mean_beta"]
            if target_amp is None:
                target_vals = [ct_amp_means[tg] for tg in target_groups
                               if tg in ct_amp_means]
                if target_vals:
                    target_amp = sum(target_vals) / len(target_vals)

            if target_amp is not None:
                # Deltas for non-target cell types (exclude all target groups)
                non_target = {ct: abs(target_amp - m)
                              for ct, m in ct_amp_means.items()
                              if ct not in target_groups}
                if non_target:
                    worst_ct = min(non_target, key=non_target.get)
                    row["min_non_target_delta"] = round(non_target[worst_ct], 4)
                    row["worst_celltype"] = worst_ct

    # Build DMR CpG rows (all CpGs in the block, not just amplicon)
    dmr_cpg_rows = []
    celltype_cpg_rows = []
    if block and p.dmr_start is not None:
        amp_start = amp["amplicon_start"]
        amp_end = amp["amplicon_end"]
        for c in block.get("cpg_sites", []):
            delta = c.get("delta_beta")
            if delta is None:
                delta = round(c["target_mean_beta"] - c["background_mean_beta"], 4)
            in_amp = (amp_start is not None and
                      amp_start <= c["position"] <= amp_end)
            dmr_cpg_rows.append({
                "assay_id": p.assay_id,
                "seq_id": p.seq_id,
                "cpg_label": c.get("label", ""),
                "cpg_position": c["position"],
                "target_mean_beta": c["target_mean_beta"],
                "background_mean_beta": c["background_mean_beta"],
                "delta_beta": delta,
                "in_amplicon": in_amp,
            })

            # Build per-cell-type row (wide format: one row per CpG,
            # one column per atlas cell type)
            pct = c.get("per_celltype_beta", {})
            if pct:
                ct_row = {
                    "assay_id": p.assay_id,
                    "seq_id": p.seq_id,
                    "cpg_label": c.get("label", ""),
                    "cpg_position": c["position"],
                }
                ct_row.update(pct)
                celltype_cpg_rows.append(ct_row)

    return amp["_cpg_rows"], dmr_cpg_rows, celltype_cpg_rows


def primer_pair_to_row(p: PrimerPair) -> dict:
    """Convert a PrimerPair to a row dict with primer + screening fields."""
    # Determine assay type from template
    if p.template_used in ("SM", "AM"):
        assay = "M"
    elif p.template_used in ("SU", "AU"):
        assay = "U"
    else:
        assay = p.template_used

    return {
        "assay_id": p.assay_id,
        "seq_id": p.seq_id,
        "template_used": p.template_used,
        "assay": assay,
        "left_primer": p.left_primer,
        "right_primer": p.right_primer,
        "left_primer_cpg_display": p.left_primer_display,
        "right_primer_cpg_display": p.right_primer_display,
        "left_tm_C": round(p.left_tm, 2),
        "right_tm_C": round(p.right_tm, 2),
        "product_size_bp": p.product_size,
        "c_total_tail": p.c_total_tail,
        "c_total": p.c_total,
        "left_c_total": p.left_c_total,
        "right_c_total": p.right_c_total,
        "left_c_tail": p.left_c_tail,
        "right_c_tail": p.right_c_tail,
        "sense_meth_mismatch_score": p.sense_meth_mismatch_score,
        "sense_unmeth_mismatch_score": p.sense_unmeth_mismatch_score,
        "anti_meth_mismatch_score": p.anti_meth_mismatch_score,
        "anti_unmeth_mismatch_score": p.anti_unmeth_mismatch_score,
        "ePCR_verdict": p.ePCR_verdict,
        "ePCR_product_count": p.ePCR_product_count,
        "ePCR_products": p.ePCR_products,
        "ePCR_error": p.ePCR_error,
        "ePCR_elapsed_sec": p.ePCR_elapsed_sec,
        "left_structure_mfe_kcal_mol": p.left_structure_mfe,
        "right_structure_mfe_kcal_mol": p.right_structure_mfe,
        "left_self_dimer_dg": p.left_self_dimer_dg,
        "right_self_dimer_dg": p.right_self_dimer_dg,
        "heterodimer_dg": p.heterodimer_dg,
        "amplicon_mfe": p.amplicon_mfe,
        "primer_dimer_prediction": p.primer_dimer_prediction,
        "primer_dimer_end_min_dg": p.primer_dimer_end_min_dg,
        "heterodimer_full_dg": p.heterodimer_full_dg,
        "left_homodimer_dg": p.left_homodimer_dg,
        "right_homodimer_dg": p.right_homodimer_dg,
        "common_variant_score": p.common_variant_score,
        "snp_count": p.snp_count,
        "segdup_overlap": p.segdup_overlap,
        "repeat_overlap": p.repeat_overlap,
        "repeat_family": p.repeat_family,
        "amplicon_snp_count": p.amplicon_snp_count,
        "cpg_snp_count": p.cpg_snp_count,
        "amplicon_segdup_overlap": p.amplicon_segdup_overlap,
        "amplicon_repeat_overlap": p.amplicon_repeat_overlap,
        "amplicon_repeat_family": p.amplicon_repeat_family,
        "mapping_error_note": p.mapping_error_note,
        "mismatch_profile": p.mismatch_profile,
        "mismatch_detail": p.mismatch_detail,
    }


# ---- Styling helpers ----

_HEADER_FONT = Font(bold=True, size=10, color="FFFFFF")
_HEADER_FILL = PatternFill(start_color="4472C4", end_color="4472C4", fill_type="solid")
_DMR_HEADER_FILL = PatternFill(start_color="2E75B6", end_color="2E75B6", fill_type="solid")
_AMP_HEADER_FILL = PatternFill(start_color="548235", end_color="548235", fill_type="solid")
_DMR_CPG_HEADER_FILL = PatternFill(start_color="7030A0", end_color="7030A0", fill_type="solid")
_CPG_CELLTYPE_HEADER_FILL = PatternFill(start_color="2E9BA0", end_color="2E9BA0", fill_type="solid")
_THIN_BORDER = Border(
    left=Side(style='thin'), right=Side(style='thin'),
    top=Side(style='thin'), bottom=Side(style='thin')
)
_PASS_FILL = PatternFill(start_color="C6EFCE", end_color="C6EFCE", fill_type="solid")
_FAIL_FILL = PatternFill(start_color="FFC7CE", end_color="FFC7CE", fill_type="solid")

_COL_WIDTHS = {
    "assay_id": 18, "seq_id": 18,
    "cell_type": 20, "chrom": 8, "dmr_start": 10, "dmr_end": 10,
    "dmr_length_bp": 12, "gene": 12, "direction": 8, "dmr_rank": 8,
    "delta_means": 10, "target_mean_meth": 12, "background_mean_meth": 12,
    "dms_score": 12, "dmr_num_cpgs": 10, "annotation": 14,
    "delta_quants": 10, "delta_maxmin": 10,
    "ttest_pval": 14, "mwtest_pval": 14, "mvalue_ttest": 14,
    "amplicon_start": 10, "amplicon_end": 10, "cpgs_in_amplicon": 10,
    "amp_target_mean_beta": 14, "amp_background_mean_beta": 14, "amp_delta": 10,
    "min_non_target_delta": 14, "worst_celltype": 20,
    "template_used": 12, "assay": 8,
    "left_primer": 28, "right_primer": 28,
    "left_primer_cpg_display": 28, "right_primer_cpg_display": 28,
    "left_tm_C": 10, "right_tm_C": 10, "product_size_bp": 12,
    "c_total_tail": 10, "c_total": 8, "left_c_total": 10, "right_c_total": 10,
    "left_c_tail": 8, "right_c_tail": 8,
    "sense_meth_mismatch_score": 12, "sense_unmeth_mismatch_score": 12,
    "anti_meth_mismatch_score": 12, "anti_unmeth_mismatch_score": 12,
    "left_structure_mfe_kcal_mol": 12, "right_structure_mfe_kcal_mol": 12,
    "primer_dimer_prediction": 20, "primer_dimer_end_min_dg": 12,
    "common_variant_score": 10, "mapping_error_note": 30,
    "mismatch_profile": 12, "mismatch_detail": 50,
    "amplicon_snp_count": 10, "cpg_snp_count": 10,
    "amplicon_segdup_overlap": 12, "amplicon_repeat_overlap": 12,
    "amplicon_repeat_family": 18,
    # CpG sheet
    "cpg_position": 12, "target_mean_beta": 14, "background_mean_beta": 14,
    "delta_beta": 10,
    # DMR CpG sheet
    "cpg_label": 10, "in_amplicon": 10,
}


def _write_sheet(ws, columns, rows, header_fill=_HEADER_FILL,
                 lightweight=False):
    """Write headers and data rows to a worksheet with styling.

    Args:
        ws: Worksheet to write to
        columns: List of column names
        rows: List of row dicts
        header_fill: PatternFill for header row
        lightweight: If True, skip per-cell border/alignment for data rows.
                     Use for large sheets (>10K rows or >50 cols) where
                     openpyxl per-cell styling is prohibitively slow.
    """
    # Headers (always styled)
    for col_idx, col_name in enumerate(columns, 1):
        cell = ws.cell(row=1, column=col_idx, value=col_name)
        cell.font = _HEADER_FONT
        cell.fill = header_fill
        cell.alignment = Alignment(horizontal='center', text_rotation=45)
        cell.border = _THIN_BORDER

    # Data
    for row_idx, row in enumerate(rows, 2):
        for col_idx, col_name in enumerate(columns, 1):
            val = row.get(col_name)
            cell = ws.cell(row=row_idx, column=col_idx, value=val)

            if not lightweight:
                cell.border = _THIN_BORDER
                cell.alignment = Alignment(horizontal='center')

                # Color-code ePCR verdict on the main sheet
                if col_name == "ePCR_verdict" and val is not None:
                    cell.fill = _PASS_FILL if val == "PASS" else _FAIL_FILL

    # Auto-filter and freeze
    ws.auto_filter.ref = ws.dimensions
    ws.freeze_panes = "A2"

    # Column widths
    for col_idx, col_name in enumerate(columns, 1):
        ws.column_dimensions[get_column_letter(col_idx)].width = _COL_WIDTHS.get(col_name, 15)


def write_xlsx(primer_pairs: List[PrimerPair], output_path: str,
               dmr_blocks: Optional[List[dict]] = None,
               sheet_name: str = "Primer pairs") -> str:
    """
    Write primer pairs to an Excel file with coupled WGBS data.

    Sheet 1 "Primer pairs": one row per primer pair, with DMR metadata
    and amplicon CpG summary columns.
    Sheet 2 "Amplicon CpGs": one row per CpG within each amplicon, with
    per-CpG target/background beta values.

    Args:
        primer_pairs: List of PrimerPair objects
        output_path: Output file path
        dmr_blocks: Optional list of DMR block dicts (from dmr_blocks.json).
                    When provided, DMR and amplicon columns are populated.
        sheet_name: Excel sheet name for the main sheet

    Returns:
        Path to the written file
    """
    dmr_lookup = _build_dmr_lookup(dmr_blocks)

    # Build enriched rows + collect per-CpG rows
    rows = []
    all_cpg_rows = []
    all_dmr_cpg_rows = []
    all_celltype_cpg_rows = []
    for p in primer_pairs:
        row = primer_pair_to_row(p)
        block = dmr_lookup.get(p.seq_id)
        cpg_rows, dmr_cpg_rows, celltype_cpg_rows = _enrich_row_with_dmr(row, p, block)
        rows.append(row)
        all_cpg_rows.extend(cpg_rows)
        all_dmr_cpg_rows.extend(dmr_cpg_rows)
        all_celltype_cpg_rows.extend(celltype_cpg_rows)

    # Create workbook
    wb = Workbook()
    wb.remove(wb.active)

    # Sheet 1: Primer pairs
    ws1 = wb.create_sheet(sheet_name)
    _write_sheet(ws1, COLUMNS, rows)

    # Sheet 2: Amplicon CpGs (CpGs within each amplicon)
    if all_cpg_rows:
        ws2 = wb.create_sheet("Amplicon CpGs")
        _write_sheet(ws2, CPG_COLUMNS, all_cpg_rows, header_fill=_AMP_HEADER_FILL)

    # Sheet 3: DMR CpGs (all CpGs in each DMR block)
    if all_dmr_cpg_rows:
        ws3 = wb.create_sheet("DMR CpGs")
        _write_sheet(ws3, DMR_CPG_COLUMNS, all_dmr_cpg_rows,
                     header_fill=_DMR_CPG_HEADER_FILL, lightweight=True)

    # Sheet 4: CpG CellType Beta (per-CpG per-cell-type mean beta, wide format)
    if all_celltype_cpg_rows:
        # Determine cell type columns dynamically from union of all keys
        celltype_keys = set()
        for r in all_celltype_cpg_rows:
            celltype_keys.update(r.keys())
        meta_cols = {"assay_id", "seq_id", "cpg_label", "cpg_position"}
        sorted_celltypes = sorted(celltype_keys - meta_cols)
        cpg_celltype_columns = ["assay_id", "seq_id", "cpg_label", "cpg_position"] + sorted_celltypes

        ws4 = wb.create_sheet("CpG CellType Beta")
        _write_sheet(ws4, cpg_celltype_columns, all_celltype_cpg_rows,
                     header_fill=_CPG_CELLTYPE_HEADER_FILL, lightweight=True)

    # Save (workspace first for S3-backed paths)
    if output_path.startswith("/mnt/"):
        tmp_path = "/workspace/_tmp_xlsx.xlsx"
        wb.save(tmp_path)
        import shutil
        shutil.copy(tmp_path, output_path)
        os.unlink(tmp_path)
    else:
        wb.save(output_path)

    return output_path


def main():
    """CLI entry point for testing."""
    import argparse
    parser = argparse.ArgumentParser(description="Write U-assays-style XLSX")
    parser.add_argument("--output", required=True, help="Output XLSX path")
    parser.add_argument("--primers", help="primers.json path")
    parser.add_argument("--dmr-blocks", help="dmr_blocks.json path")
    args = parser.parse_args()

    import json
    from methyl_panel.phase3_primer3_design import PrimerPair

    # Load primers
    with open(args.primers) as f:
        primers_data = json.load(f)

    # Load DMR blocks
    dmr_blocks = None
    if args.dmr_blocks:
        with open(args.dmr_blocks) as f:
            dmr_blocks = json.load(f)

    # Convert to PrimerPair objects
    pairs = []
    for p in primers_data:
        pairs.append(PrimerPair(
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
            left_self_any=0, right_self_any=0, left_self_end=0, right_self_end=0,
            pair_compl_any=0, pair_compl_end=0,
            left_end_stability=0, right_end_stability=0,
            penalty=p["penalty"],
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
            common_variant_score=p.get("common_variant_score"),
            snp_count=p.get("snp_count"),
            segdup_overlap=p.get("segdup_overlap"),
            repeat_overlap=p.get("repeat_overlap"),
            repeat_family=p.get("repeat_family"),
            amplicon_snp_count=p.get("amplicon_snp_count"),
            cpg_snp_count=p.get("cpg_snp_count"),
            cpg_snp_details=p.get("cpg_snp_details"),
            amplicon_segdup_overlap=p.get("amplicon_segdup_overlap"),
            amplicon_repeat_overlap=p.get("amplicon_repeat_overlap"),
            amplicon_repeat_family=p.get("amplicon_repeat_family"),
            mapping_error_note=p.get("mapping_error_note"),
            mismatch_profile=p.get("mismatch_profile"),
            mismatch_detail=p.get("mismatch_detail"),
            chrom=p.get("chrom"),
            dmr_start=p.get("dmr_start"),
            dmr_end=p.get("dmr_end"),
            direction=p.get("direction", "hypo"),
        ))

    write_xlsx(pairs, args.output, dmr_blocks=dmr_blocks)
    print(f"Written to {args.output} ({len(pairs)} primers, {len(dmr_blocks) if dmr_blocks else 0} DMR blocks)")


if __name__ == "__main__":
    main()
