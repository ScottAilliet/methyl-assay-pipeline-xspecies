#!/usr/bin/env python3
"""
Generate a master Excel file combining all primer assays across all cell types.

Scans results/*/ directories for primers.json and dmr_blocks.json,
joins them on seq_id, and writes a single XLSX with multiple sheets:

  - All_Primers:    Every primer pair (~4,335 rows), enriched with DMR-level info
  - Best_Per_DMR:   Top primer per DMR block (~760 rows), ranked by ePCR pass + dimer + penalty
  - Summary:        Per cell-type-config counts (blocks, primers, ePCR pass/fail, dimer)
  - DMR_Blocks:     One row per DMR block with primer count
  - Per cell type:  Individual sheets for each cell-type config

Usage:
    python3 scripts/generate_master_excel.py [--results-dir results/] [--output master_primer_assays.xlsx]
"""

import os
import sys
import json
import argparse
from typing import Dict, List
from collections import defaultdict

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter


# ── Column definitions ──────────────────────────────────────────────────────

DMR_COLUMNS = [
    "cell_type", "source", "direction", "rank", "chrom", "dmr_start", "dmr_end",
    "gene", "annotation", "delta_means", "dms_score", "cleanliness_score",
    "target_mean_meth", "background_mean_meth", "num_cpgs", "block_len",
]

PRIMER_COLUMNS = [
    "assay_id", "seq_id", "template_used", "assay",
    "left_primer", "right_primer",
    "left_primer_cpg_display", "right_primer_cpg_display",
    "left_tm_C", "right_tm_C", "product_size_bp",
    "c_total", "left_c_total", "right_c_total",
    "c_total_tail", "left_c_tail", "right_c_tail",
    "penalty", "left_gc_percent", "right_gc_percent",
    "sense_meth_mismatch_score", "sense_unmeth_mismatch_score",
    "anti_meth_mismatch_score", "anti_unmeth_mismatch_score",
    "ePCR_verdict", "ePCR_product_count", "ePCR_products", "ePCR_error", "ePCR_elapsed_sec",
    "mismatch_profile", "mismatch_detail", "mapping_error_note",
    "left_structure_mfe_kcal_mol", "right_structure_mfe_kcal_mol",
    "left_self_dimer_dg", "right_self_dimer_dg", "heterodimer_dg",
    "amplicon_mfe",
    "primer_dimer_prediction", "primer_dimer_end_min_dg",
    "heterodimer_full_dg", "left_homodimer_dg", "right_homodimer_dg",
    "common_variant_score", "snp_count",
    "segdup_overlap", "repeat_overlap", "repeat_family",
    "amplicon_snp_count", "cpg_snp_count",
    "amplicon_segdup_overlap", "amplicon_repeat_overlap", "amplicon_repeat_family",
]

ALL_COLUMNS = DMR_COLUMNS + PRIMER_COLUMNS

COL_WIDTHS = {
    "cell_type": 10, "source": 10, "direction": 8, "rank": 6,
    "chrom": 8, "dmr_start": 12, "dmr_end": 12,
    "gene": 12, "annotation": 15, "delta_means": 10,
    "dms_score": 10, "cleanliness_score": 10, "target_mean_meth": 12, "background_mean_meth": 14,
    "num_cpgs": 8, "block_len": 9,
    "assay_id": 18, "seq_id": 14, "template_used": 12, "assay": 6,
    "left_primer": 28, "right_primer": 28,
    "left_primer_cpg_display": 28, "right_primer_cpg_display": 28,
    "left_tm_C": 10, "right_tm_C": 10, "product_size_bp": 12,
    "c_total": 8, "left_c_total": 10, "right_c_total": 10,
    "c_total_tail": 10, "left_c_tail": 8, "right_c_tail": 8,
    "penalty": 8, "left_gc_percent": 10, "right_gc_percent": 10,
    "sense_meth_mismatch_score": 12, "sense_unmeth_mismatch_score": 12,
    "anti_meth_mismatch_score": 12, "anti_unmeth_mismatch_score": 12,
    "ePCR_verdict": 10, "ePCR_product_count": 10, "ePCR_products": 30, "ePCR_error": 20, "ePCR_elapsed_sec": 10,
    "mismatch_profile": 12, "mismatch_detail": 50, "mapping_error_note": 30,
    "left_structure_mfe_kcal_mol": 12, "right_structure_mfe_kcal_mol": 12,
    "left_self_dimer_dg": 12, "right_self_dimer_dg": 12, "heterodimer_dg": 12,
    "amplicon_mfe": 12,
    "primer_dimer_prediction": 20, "primer_dimer_end_min_dg": 12,
    "heterodimer_full_dg": 12, "left_homodimer_dg": 12, "right_homodimer_dg": 12,
    "common_variant_score": 10, "snp_count": 8,
    "segdup_overlap": 10, "repeat_overlap": 10, "repeat_family": 16,
    "amplicon_snp_count": 10, "cpg_snp_count": 10,
    "amplicon_segdup_overlap": 12, "amplicon_repeat_overlap": 12,
    "amplicon_repeat_family": 18,
}

SUMMARY_COLUMNS = [
    "cell_type", "source", "direction", "num_blocks", "num_primers",
    "num_pass_epcr", "num_fail_epcr", "num_epcr_none",
    "num_low_dimer", "num_medium_dimer", "num_high_dimer",
]

DMR_BLOCK_COLUMNS = [
    "cell_type", "source", "direction", "seq_id", "rank",
    "chrom", "dmr_start", "dmr_end", "gene", "annotation",
    "delta_means", "dms_score", "cleanliness_score", "target_mean_meth",
    "background_mean_meth", "num_cpgs", "block_len",
    "num_primers", "has_epcr_pass",
]


# ── Data loading ────────────────────────────────────────────────────────────

def _source_label(dir_name: str) -> str:
    if dir_name.endswith("_atlas"):
        return "atlas"
    elif dir_name.endswith("_relaxed"):
        return "relaxed"
    return "v2.2.8"


def _cell_type_from_dir(dir_name: str) -> str:
    for suffix in ("_atlas", "_relaxed"):
        if dir_name.endswith(suffix):
            return dir_name[:-len(suffix)]
    return dir_name


def load_results_dir(results_dir: str) -> List[dict]:
    """Scan all subdirectories, load primers.json + dmr_blocks.json, join on seq_id."""
    all_rows = []

    for entry in sorted(os.listdir(results_dir)):
        dir_path = os.path.join(results_dir, entry)
        if not os.path.isdir(dir_path):
            continue

        primers_path = os.path.join(dir_path, "primers.json")
        blocks_path = os.path.join(dir_path, "dmr_blocks.json")

        if not os.path.exists(primers_path):
            continue

        with open(primers_path) as f:
            primers = json.load(f)
        if not primers:
            continue

        block_lookup: Dict[str, dict] = {}
        if os.path.exists(blocks_path):
            with open(blocks_path) as f:
                blocks = json.load(f)
            for b in blocks:
                block_lookup[b["seq_id"]] = b

        cell_type = _cell_type_from_dir(entry)
        source = _source_label(entry)

        for p in primers:
            row = _build_primer_row(p, block_lookup, cell_type, source)
            all_rows.append(row)

        print(f"  {entry}: {len(primers)} primers, {len(block_lookup)} blocks")

    return all_rows


def _build_primer_row(p: dict, block_lookup: Dict[str, dict],
                      cell_type: str, source: str) -> dict:
    seq_id = p.get("seq_id", "")
    block = block_lookup.get(seq_id, {})

    template = p.get("template_used", "")
    if template in ("SM", "AM"):
        assay = "M"
    elif template in ("SU", "AU"):
        assay = "U"
    else:
        assay = template

    return {
        # DMR-level
        "cell_type": cell_type, "source": source,
        "direction": p.get("direction") or block.get("direction", "hypo"),
        "rank": block.get("rank"), "chrom": block.get("chrom"),
        "dmr_start": block.get("start"), "dmr_end": block.get("end"),
        "gene": block.get("gene"), "annotation": block.get("annotation"),
        "delta_means": block.get("delta_means"),
        "dms_score": block.get("dms_score"),
        "cleanliness_score": block.get("cleanliness_score"),
        "target_mean_meth": block.get("target_mean_meth"),
        "background_mean_meth": block.get("background_mean_meth"),
        "num_cpgs": block.get("num_cpgs"), "block_len": block.get("block_len"),
        # Primer-level
        "assay_id": p.get("assay_id"), "seq_id": seq_id,
        "template_used": template, "assay": assay,
        "left_primer": p.get("left_primer"), "right_primer": p.get("right_primer"),
        "left_primer_cpg_display": p.get("left_primer_display"),
        "right_primer_cpg_display": p.get("right_primer_display"),
        "left_tm_C": round(p["left_tm"], 2) if p.get("left_tm") is not None else None,
        "right_tm_C": round(p["right_tm"], 2) if p.get("right_tm") is not None else None,
        "product_size_bp": p.get("product_size"),
        "c_total": p.get("c_total"), "left_c_total": p.get("left_c_total"),
        "right_c_total": p.get("right_c_total"),
        "c_total_tail": p.get("c_total_tail"),
        "left_c_tail": p.get("left_c_tail"), "right_c_tail": p.get("right_c_tail"),
        "penalty": p.get("penalty"),
        "left_gc_percent": p.get("left_gc_percent"),
        "right_gc_percent": p.get("right_gc_percent"),
        "sense_meth_mismatch_score": p.get("sense_meth_mismatch_score"),
        "sense_unmeth_mismatch_score": p.get("sense_unmeth_mismatch_score"),
        "anti_meth_mismatch_score": p.get("anti_meth_mismatch_score"),
        "anti_unmeth_mismatch_score": p.get("anti_unmeth_mismatch_score"),
        "ePCR_verdict": p.get("ePCR_verdict"),
        "ePCR_product_count": p.get("ePCR_product_count"),
        "ePCR_products": p.get("ePCR_products"),
        "ePCR_error": p.get("ePCR_error"),
        "ePCR_elapsed_sec": p.get("ePCR_elapsed_sec"),
        "mismatch_profile": p.get("mismatch_profile"),
        "mismatch_detail": p.get("mismatch_detail"),
        "mapping_error_note": p.get("mapping_error_note"),
        "left_structure_mfe_kcal_mol": p.get("left_structure_mfe"),
        "right_structure_mfe_kcal_mol": p.get("right_structure_mfe"),
        "left_self_dimer_dg": p.get("left_self_dimer_dg"),
        "right_self_dimer_dg": p.get("right_self_dimer_dg"),
        "heterodimer_dg": p.get("heterodimer_dg"),
        "amplicon_mfe": p.get("amplicon_mfe"),
        "primer_dimer_prediction": p.get("primer_dimer_prediction"),
        "primer_dimer_end_min_dg": p.get("primer_dimer_end_min_dg"),
        "heterodimer_full_dg": p.get("heterodimer_full_dg"),
        "left_homodimer_dg": p.get("left_homodimer_dg"),
        "right_homodimer_dg": p.get("right_homodimer_dg"),
        "common_variant_score": p.get("common_variant_score"),
        "snp_count": p.get("snp_count"),
        "segdup_overlap": p.get("segdup_overlap"),
        "repeat_overlap": p.get("repeat_overlap"),
        "repeat_family": p.get("repeat_family"),
        "amplicon_snp_count": p.get("amplicon_snp_count"),
        "cpg_snp_count": p.get("cpg_snp_count"),
        "amplicon_segdup_overlap": p.get("amplicon_segdup_overlap"),
        "amplicon_repeat_overlap": p.get("amplicon_repeat_overlap"),
        "amplicon_repeat_family": p.get("amplicon_repeat_family"),
    }


# ── Best-per-DMR selection ──────────────────────────────────────────────────

def select_best_per_dmr(all_rows: List[dict]) -> List[dict]:
    """
    Best primer per DMR block. Ranking:
      1. ePCR_verdict (PASS > FAIL > None)
      2. primer_dimer_prediction (low > medium > high)
      3. lower penalty
    """
    groups = defaultdict(list)
    for row in all_rows:
        key = (row["cell_type"], row["source"],
               row.get("direction", "hypo"), row["seq_id"])
        groups[key].append(row)

    dimer_rank = {"low": 0, "medium": 1, "high": 2}
    epcr_rank = {"PASS": 0, "FAIL": 1, None: 2}

    best_rows = []
    for key, rows in groups.items():
        best = min(rows, key=lambda r: (
            epcr_rank.get(r.get("ePCR_verdict"), 2),
            dimer_rank.get(r.get("primer_dimer_prediction"), 3),
            r.get("penalty") or 999,
        ))
        best_rows.append(best)

    best_rows.sort(key=lambda r: (
        r.get("cell_type", ""), r.get("source", ""), r.get("rank") or 999,
    ))
    return best_rows


# ── Summary + DMR blocks ────────────────────────────────────────────────────

def build_summary(all_rows: List[dict]) -> List[dict]:
    stats = defaultdict(lambda: {
        "num_primers": 0, "num_pass_epcr": 0, "num_fail_epcr": 0,
        "num_epcr_none": 0, "num_low_dimer": 0, "num_medium_dimer": 0,
        "num_high_dimer": 0, "seq_ids": set(),
    })

    for r in all_rows:
        key = (r["cell_type"], r["source"], r.get("direction", "hypo"))
        s = stats[key]
        s["num_primers"] += 1
        s["seq_ids"].add(r["seq_id"])

        verdict = r.get("ePCR_verdict")
        if verdict == "PASS":
            s["num_pass_epcr"] += 1
        elif verdict and verdict.startswith("FAIL"):
            s["num_fail_epcr"] += 1
        else:
            s["num_epcr_none"] += 1

        dim = r.get("primer_dimer_prediction")
        if dim == "low":
            s["num_low_dimer"] += 1
        elif dim == "medium":
            s["num_medium_dimer"] += 1
        elif dim == "high":
            s["num_high_dimer"] += 1

    rows = []
    for (cell_type, source, direction), s in sorted(stats.items()):
        rows.append({
            "cell_type": cell_type, "source": source,
            "direction": direction,
            "num_blocks": len(s["seq_ids"]), "num_primers": s["num_primers"],
            "num_pass_epcr": s["num_pass_epcr"],
            "num_fail_epcr": s["num_fail_epcr"],
            "num_epcr_none": s["num_epcr_none"],
            "num_low_dimer": s["num_low_dimer"],
            "num_medium_dimer": s["num_medium_dimer"],
            "num_high_dimer": s["num_high_dimer"],
        })
    return rows


def build_dmr_blocks(all_rows: List[dict]) -> List[dict]:
    block_stats = defaultdict(lambda: {"num_primers": 0, "has_epcr_pass": False})
    block_info = {}

    for r in all_rows:
        key = (r["cell_type"], r["source"], r.get("direction", "hypo"), r["seq_id"])
        block_stats[key]["num_primers"] += 1
        if r.get("ePCR_verdict") == "PASS":
            block_stats[key]["has_epcr_pass"] = True
        if key not in block_info:
            block_info[key] = {
                "cell_type": r["cell_type"], "source": r["source"],
                "direction": r.get("direction", "hypo"),
                "seq_id": r["seq_id"], "rank": r.get("rank"),
                "chrom": r.get("chrom"),
                "dmr_start": r.get("dmr_start"), "dmr_end": r.get("dmr_end"),
                "gene": r.get("gene"), "annotation": r.get("annotation"),
                "delta_means": r.get("delta_means"),
                "dms_score": r.get("dms_score"),
                "cleanliness_score": r.get("cleanliness_score"),
                "target_mean_meth": r.get("target_mean_meth"),
                "background_mean_meth": r.get("background_mean_meth"),
                "num_cpgs": r.get("num_cpgs"), "block_len": r.get("block_len"),
            }

    rows = []
    for key, info in block_info.items():
        row = dict(info)
        row["num_primers"] = block_stats[key]["num_primers"]
        row["has_epcr_pass"] = block_stats[key]["has_epcr_pass"]
        rows.append(row)

    rows.sort(key=lambda r: (
        r.get("cell_type", ""), r.get("source", ""), r.get("rank") or 999,
    ))
    return rows


# ── XLSX writing ────────────────────────────────────────────────────────────

def _style_sheet(ws, columns: List[str], df: pd.DataFrame, lightweight: bool = False):
    header_font = Font(bold=True, size=10, color="FFFFFF")
    header_fill = PatternFill(start_color="4472C4", end_color="4472C4", fill_type="solid")
    thin_border = Border(
        left=Side(style='thin'), right=Side(style='thin'),
        top=Side(style='thin'), bottom=Side(style='thin')
    )
    pass_fill = PatternFill(start_color="C6EFCE", end_color="C6EFCE", fill_type="solid")
    fail_fill = PatternFill(start_color="FFC7CE", end_color="FFC7CE", fill_type="solid")

    # Header row — always styled
    for col_idx, col_name in enumerate(columns, 1):
        cell = ws.cell(row=1, column=col_idx, value=col_name)
        cell.font = header_font
        cell.fill = header_fill
        cell.alignment = Alignment(horizontal='center', text_rotation=45)

    # Data rows
    for row_idx, (_, row_data) in enumerate(df.iterrows(), 2):
        for col_idx, col_name in enumerate(columns, 1):
            val = row_data.get(col_name)
            cell = ws.cell(row=row_idx, column=col_idx, value=val)

            # Skip per-cell border + alignment in lightweight mode (saves ~10x memory)
            if not lightweight:
                cell.border = thin_border
                cell.alignment = Alignment(horizontal='center')

            # Conditional fills always applied (cheap, high-value)
            if col_name == "ePCR_verdict" and val is not None:
                cell.fill = pass_fill if val == "PASS" else fail_fill
            if col_name == "primer_dimer_prediction" and val is not None:
                if val == "high":
                    cell.fill = fail_fill
                elif val == "low":
                    cell.fill = pass_fill

    ws.auto_filter.ref = ws.dimensions
    ws.freeze_panes = "A2"

    for col_idx, col_name in enumerate(columns, 1):
        ws.column_dimensions[get_column_letter(col_idx)].width = COL_WIDTHS.get(col_name, 15)


def write_master_xlsx(all_rows: List[dict], output_path: str,
                      per_celltype: bool = False):
    wb = Workbook()
    wb.remove(wb.active)

    # Use lightweight styling (no per-cell borders/alignment) for large sheets
    # to avoid OOM. Threshold: >5000 rows gets lightweight mode.
    LIGHTWEIGHT_THRESHOLD = 5000

    # Sheet 1: All_Primers
    df_all = pd.DataFrame(all_rows, columns=ALL_COLUMNS)
    lw = len(df_all) > LIGHTWEIGHT_THRESHOLD
    ws = wb.create_sheet("All_Primers")
    _style_sheet(ws, ALL_COLUMNS, df_all, lightweight=lw)
    print(f"  All_Primers: {len(df_all)} rows{' (lightweight)' if lw else ''}")

    # Sheet 2: Best_Per_DMR
    best_rows = select_best_per_dmr(all_rows)
    df_best = pd.DataFrame(best_rows, columns=ALL_COLUMNS)
    lw = len(df_best) > LIGHTWEIGHT_THRESHOLD
    ws = wb.create_sheet("Best_Per_DMR")
    _style_sheet(ws, ALL_COLUMNS, df_best, lightweight=lw)
    print(f"  Best_Per_DMR: {len(df_best)} rows{' (lightweight)' if lw else ''}")

    # Sheet 3: Summary
    summary_rows = build_summary(all_rows)
    df_summary = pd.DataFrame(summary_rows, columns=SUMMARY_COLUMNS)
    ws = wb.create_sheet("Summary")
    _style_sheet(ws, SUMMARY_COLUMNS, df_summary)
    print(f"  Summary: {len(df_summary)} rows")

    # Sheet 4: DMR_Blocks
    dmr_rows = build_dmr_blocks(all_rows)
    df_dmr = pd.DataFrame(dmr_rows, columns=DMR_BLOCK_COLUMNS)
    ws = wb.create_sheet("DMR_Blocks")
    _style_sheet(ws, DMR_BLOCK_COLUMNS, df_dmr)
    print(f"  DMR_Blocks: {len(df_dmr)} rows")

    # Sheets 5-N: Per cell type (opt-in — redundant with All_Primers + auto_filter)
    if per_celltype:
        cell_types = sorted(df_all["cell_type"].unique())
        sources = sorted(df_all["source"].unique())
        directions = sorted(df_all["direction"].dropna().unique()) if "direction" in df_all.columns else ["hypo"]
        for ct in cell_types:
            for src in sources:
                for d in directions:
                    subset = df_all[
                        (df_all["cell_type"] == ct) &
                        (df_all["source"] == src) &
                        (df_all.get("direction", "hypo") == d)
                    ]
                    if len(subset) == 0:
                        continue
                    parts = [ct]
                    if src != "v2.2.8":
                        parts.append(src)
                    if d != "hypo":
                        parts.append(d)
                    sheet_name = "_".join(parts)[:31]
                    ws = wb.create_sheet(sheet_name)
                    lw = len(subset) > LIGHTWEIGHT_THRESHOLD
                    _style_sheet(ws, ALL_COLUMNS, subset, lightweight=lw)
                    print(f"  {sheet_name}: {len(subset)} rows")
    else:
        print(f"  (Per-cell-type sheets skipped — use --per-celltype to enable)")

    wb.save(output_path)
    print(f"\nMaster Excel saved to {output_path}")


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate master Excel combining all primer assays across cell types"
    )
    parser.add_argument("--results-dir", default="results/",
                        help="Directory containing per-cell-type result subdirectories")
    parser.add_argument("--output", default="master_primer_assays.xlsx",
                        help="Output XLSX path")
    parser.add_argument("--per-celltype", action="store_true",
                        help="Include per-cell-type sheets (default: skip — use All_Primers auto_filter instead)")
    args = parser.parse_args()

    if not os.path.isdir(args.results_dir):
        print(f"Error: results directory '{args.results_dir}' not found")
        sys.exit(1)

    print(f"Scanning {args.results_dir} for primer results...")
    all_rows = load_results_dir(args.results_dir)

    if not all_rows:
        print("Error: no primers found. Run the pipeline first.")
        sys.exit(1)

    print(f"\nTotal: {len(all_rows)} primer pairs across all cell types")
    print(f"\nWriting master Excel to {args.output}...")
    write_master_xlsx(all_rows, args.output, per_celltype=args.per_celltype)


if __name__ == "__main__":
    main()
