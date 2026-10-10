#!/bin/bash
set -e

# ============================================================================
# run_full_pipeline.sh
# Consolidated pipeline: git pull, re-select, compute DMS sub-scores,
# annotate primers, extract CpG beta, merge into final atlas.
#
# Does NOT re-run ePCR screening (already done, results in results_v2.csv).
# Only re-runs the selection step with the new "keep all passing primers" logic.
#
# Usage:
#   1. Edit the USER CONFIG section below with your local paths.
#   2. Run: bash run_full_pipeline.sh
# ============================================================================

# --- USER CONFIG (edit these paths) ---
EPCR_DIR=~/epcr_screening
REPO_DIR=~/cell-specific-epigenetic-dpcr-assay-design-pipeline
BETA_DIR=~/path/to/beta/files
GROUPS_CSV=~/path/to/groups.csv
CPG_INDEX_FILE=~/path/to/CpG.bed.gz
# --- END USER CONFIG ---

echo "=== Step 1: Pull latest code ==="
cd "$REPO_DIR"
git pull

echo ""
echo "=== Step 2: Copy scripts to working directory ==="
cp bisearch_epcr.py "$EPCR_DIR/"
cp scripts/annotate_primers.py "$EPCR_DIR/"
cp scripts/compute_dms_subscores.py "$EPCR_DIR/"
cp scripts/extract_cpg_beta.py "$EPCR_DIR/"
cp scripts/merge_screening_results.py "$EPCR_DIR/"
cp scripts/data_dictionary.py "$EPCR_DIR/"

echo ""
echo "=== Step 3: Re-run selection on existing screening results ==="
echo "(Keeps ALL passing primers - no DMR dedup, no cell-type cap)"
cd "$EPCR_DIR"
python bisearch_epcr.py select \
  --input results_v2.csv \
  --output results_v2_selected.csv

echo ""
echo "=== Step 4: Compute DMS sub-scores ==="
python compute_dms_subscores.py \
  --atlas master_primer_assays_atlas.xlsx \
  --beta-dir "$BETA_DIR" \
  --groups "$GROUPS_CSV" \
  --cpg-index "$CPG_INDEX_FILE" \
  --output dms_subscores.csv

echo ""
echo "=== Step 5: Annotate primers (SNPs, repeats, segdups) ==="
python annotate_primers.py \
  --atlas master_primer_assays_atlas.xlsx \
  --results results_v2.csv \
  --genome-dir hg19_bs \
  --output primer_annotations.csv

echo ""
echo "=== Step 6: Extract per-CpG beta values ==="
python extract_cpg_beta.py \
  --results results_v2.csv \
  --beta-dir "$BETA_DIR" \
  --groups "$GROUPS_CSV" \
  --cpg-index "$CPG_INDEX_FILE" \
  --output cpg_beta_values.csv

echo ""
echo "=== Step 7: Merge everything into final atlas ==="
python merge_screening_results.py \
  --atlas master_primer_assays_atlas.xlsx \
  --results results_v2.csv \
  --selected results_v2_selected.csv \
  --dms-subscores dms_subscores.csv \
  --annotations primer_annotations.csv \
  --output master_atlas_screened_v3.xlsx

echo ""
echo "=== Done! ==="
echo "Output files in $EPCR_DIR:"
echo "  master_atlas_screened_v3.xlsx  - Final merged atlas"
echo "  cpg_beta_values.csv            - Per-CpG beta values for all cell types"
echo "  dms_subscores.csv              - DMS A-E sub-scores"
echo "  primer_annotations.csv         - SNP/repeat/segdup annotations"
