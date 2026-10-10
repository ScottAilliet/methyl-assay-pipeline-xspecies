#!/bin/bash
# =============================================================================
# Cell-Specific Epigenetic dPCR Assay Design Pipeline -- Annotate Existing primers.json with SNP/Repeat/Segdup
# =============================================================================
# Runs only pipeline steps 6 (SNP screening) and 61 (repeat/segdup screening)
# on existing primers.json files. Adds:
#
#   Primer-level:
#     - common_variant_score, snp_count        (from dbSNP common variants)
#     - segdup_overlap                         (from UCSC genomicSuperDups)
#     - repeat_overlap, repeat_family          (from UCSC RepeatMasker)
#
#   Amplicon-level (new -- screens full amplicon, not just primer binding sites):
#     - amplicon_snp_count                     (total common SNPs in amplicon)
#     - cpg_snp_count                          (SNPs at measured CpG positions)
#     - cpg_snp_details                        (per-SNP detail with is_c_to_t flag)
#     - amplicon_segdup_overlap                (segdup overlap in amplicon)
#     - amplicon_repeat_overlap                (repeat overlap in amplicon)
#     - amplicon_repeat_family                 (repeat family in amplicon)
#
# This is much faster than re-running the full pipeline (~seconds per cell type
# vs ~5-15 min) because it skips DMR discovery, primer design, and ePCR screening.
#
# Prerequisites:
#   - data/dbsnp_common.vcf.gz (+ .tbi)     -- run ./download_data.sh
#   - data/segdup.bed.gz (+ .tbi)           -- run ./download_data.sh
#   - data/rmsk.bed.gz (+ .tbi)             -- run ./download_data.sh
#   - tabix (part of htslib)
#   - dmr_blocks.json in each results subdirectory (for CpG positions)
#
# Usage:
#   ./annotate_primers.sh                          # annotate all in results_atlas/
#   ./annotate_primers.sh --results-dir results    # specify results directory
#   ./annotate_primers.sh --resume                 # skip already-annotated cell types
# =============================================================================

set -euo pipefail

# --- Defaults ---
RESULTS_DIR="results_atlas"
DBSNP="data/dbsnp_common.vcf.gz"
SEGDUP_BED="data/segdup.bed.gz"
REPEAT_BED="data/rmsk.bed.gz"
RESUME=false
FLANK=100

# --- Parse arguments ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --results-dir)  RESULTS_DIR="$2"; shift 2 ;;
        --resume)       RESUME=true; shift ;;
        --dbsnp)        DBSNP="$2"; shift 2 ;;
        --segdup-bed)   SEGDUP_BED="$2"; shift 2 ;;
        --repeat-bed)   REPEAT_BED="$2"; shift 2 ;;
        --help|-h)
            echo "Usage: ./annotate_primers.sh [--results-dir DIR] [--resume]"
            echo "                         [--dbsnp PATH] [--segdup-bed PATH] [--repeat-bed PATH]"
            exit 0 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# --- Prerequisite checks ---
echo "=========================================="
echo "Primer Annotation (SNP + Repeat + Segdup)"
echo "=========================================="
echo ""

MISSING=0

for f in "$DBSNP" "$SEGDUP_BED" "$REPEAT_BED"; do
    if [ ! -f "$f" ]; then
        echo "MISSING: $f"
        MISSING=$((MISSING + 1))
    fi
done

# Check tabix indices
for f in "$DBSNP" "$SEGDUP_BED" "$REPEAT_BED"; do
    if [ -f "$f" ] && [ ! -f "${f}.tbi" ]; then
        echo "MISSING INDEX: ${f}.tbi"
        MISSING=$((MISSING + 1))
    fi
done

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "ERROR: $MISSING file(s) missing."
    echo "Run ./download_data.sh first to download annotation data."
    exit 1
fi

# Check tabix is available
if ! command -v tabix &>/dev/null; then
    echo "ERROR: tabix not found in PATH."
    echo "Install with: brew install htslib"
    exit 1
fi

echo "  dbSNP:     $DBSNP"
echo "  Segdup:    $SEGDUP_BED"
echo "  Repeat:    $REPEAT_BED"
echo "  Results:   $RESULTS_DIR"
echo "  Resume:    $RESUME"
echo ""

# --- Find all cell type directories with primers.json ---
CELL_TYPES=()
for dir in "$RESULTS_DIR"/*/; do
    if [ -f "${dir}primers.json" ]; then
        CT=$(basename "$dir")
        CELL_TYPES+=("$CT")
    fi
done

TOTAL=${#CELL_TYPES[@]}
if [ $TOTAL -eq 0 ]; then
    echo "ERROR: No primers.json files found in $RESULTS_DIR/"
    exit 1
fi

echo "Found $TOTAL cell types with primers.json"
echo ""

# --- Annotate each cell type ---
COMPLETED=0
SKIPPED=0
FAILED=0
FAILED_LIST=()
IDX=0

for CT in "${CELL_TYPES[@]}"; do
    IDX=$((IDX + 1))
    OUT_DIR="$RESULTS_DIR/$CT"

    # Check resume -- skip if common_variant_score is already populated
    if [ "$RESUME" = true ]; then
        ALREADY=$(python3 -c "
import json, sys
with open('$OUT_DIR/primers.json') as f:
    primers = json.load(f)
if primers and primers[0].get('common_variant_score') is not None:
    print('yes')
else:
    print('no')
" 2>/dev/null || echo "no")
        if [ "$ALREADY" = "yes" ]; then
            echo "[$IDX/$TOTAL] SKIP: $CT (already annotated)"
            SKIPPED=$((SKIPPED + 1))
            continue
        fi
    fi

    echo -n "[$IDX/$TOTAL] $CT ... "

    # Run steps 6 and 61 only
    if python3 -m methyl_panel.pipeline \
        --steps 6,61 \
        --dbsnp "$DBSNP" \
        --segdup-bed "$SEGDUP_BED" \
        --repeat-bed "$REPEAT_BED" \
        --flank $FLANK \
        --output-dir "$OUT_DIR/" 2>&1 | tail -5; then

        # Count annotations for summary
        SUMMARY=$(python3 -c "
import json
with open('$OUT_DIR/primers.json') as f:
    primers = json.load(f)
n = len(primers)
snp = sum(1 for p in primers if p.get('snp_count', 0) and p['snp_count'] > 0)
cpg = sum(1 for p in primers if p.get('cpg_snp_count', 0) and p['cpg_snp_count'] > 0)
seg = sum(1 for p in primers if p.get('segdup_overlap'))
rep = sum(1 for p in primers if p.get('repeat_overlap'))
print(f'{n} primers, {snp} primer SNPs, {cpg} CpG SNPs, {seg} segdup, {rep} repeats')
" 2>/dev/null || echo "done")

        echo "OK: $SUMMARY"
        COMPLETED=$((COMPLETED + 1))
    else
        echo "FAILED"
        FAILED=$((FAILED + 1))
        FAILED_LIST+=("$CT")
    fi
done

# --- Summary ---
echo ""
echo "=========================================="
echo "Annotation complete"
echo "=========================================="
echo "  Annotated: $COMPLETED"
echo "  Skipped:   $SKIPPED"
echo "  Failed:    $FAILED"
if [ $FAILED -gt 0 ]; then
    echo "  Failed list: ${FAILED_LIST[*]}"
fi
echo ""
echo "Next step: regenerate master Excel:"
echo "  python3 scripts/generate_master_excel.py --results-dir $RESULTS_DIR --output master_primer_assays_atlas.xlsx"
