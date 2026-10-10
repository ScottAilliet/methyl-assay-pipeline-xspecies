# BiSearch ePCR - Hash-Based Primer Specificity Screening

Pure-Python reimplementation of the BiSearch ePCR algorithm for bisulfite PCR
primer specificity checking. No external aligner, BLAST, or Linux-only
dependencies required. Runs on macOS with numpy + pandas + openpyxl.

## Algorithm

Implements the 16-mer hash-based ePCR from Tusnady et al. 2005 (NAR 33:e9)
and Aranyi et al. 2006 (BMC Bioinformatics 7:431):

1. **Seed extraction**: 3'-most 16 bases of each primer
2. **Variant generation**: All seed variants within mismatch tolerance
   (mismatch string "1233333333333333" = up to 3 mismatches anywhere in seed)
3. **Hash lookup**: Sorted numpy index + binary search (np.searchsorted)
4. **Vectorized verification**: Full-primer mismatch counting at all seed hits
   simultaneously using numpy fancy indexing
5. **Product pairing**: FP (strand +) paired with RP (strand -) on same
   chromosome, same genome state, product size <= 1000 bp
6. **Deduplication**: Same locus in unmethylated + methylated states counts
   as one product

**Verdict**: PASS = exactly 1 product total; FAIL = 0 or >1 products.

Checks two bisulfite-converted genome states:
- Unmethylated: CT (all C->T) + GA (all G->A)
- Methylated: CT_M (CpG C preserved) + GA_M (CpG G preserved)

The unconverted genome is NOT checked (bisulfite PCR templates are converted).

## Setup (macOS)

```bash
# Install dependencies
pip install numpy pandas openpyxl

# Download and prepare hg19 genome (one-time, ~30 min)
python prepare_genomes.py --output-dir hg19_bs --genome hg19

# This downloads per-chromosome FASTA from UCSC and saves:
#   hg19_bs/chr1.npy ... chr22.npy, chrX.npy, chrY.npy
#   hg19_bs/chrom_info.json
```

## Usage

### Single primer pair

```bash
python bisearch_epcr.py screen \
    --fp AGAGGTTGGAGATTTTGGAGTGAGA \
    --rp TTTGTGTTTTAGTGGTATTTTATAA \
    --genome-dir hg19_bs
```

### Batch screening with early stopping

```bash
python bisearch_epcr.py screen \
    --input primers.csv \
    --output results.csv \
    --genome-dir hg19_bs \
    --max-pass-per-ct 20 \
    --screen-per-ct 200 \
    --verbose
```

Input CSV columns: `assay_id`, `fp_seq` (or `left_primer`), `rp_seq` (or
`right_primer`), `cell_type`, `delta_means`, `chrom`, `dmr_start`, `dmr_end`.

The batch mode:
1. Groups primers by cell_type, sorts by delta_means (descending)
2. Screens top `--screen-per-ct` primers per cell type (default 200)
3. Processes chromosomes one at a time (builds indexes once per chromosome)
4. Selects top `--max-pass-per-ct` passing assays per cell type (default 20)
5. Writes two files: `results.csv` (all screened) and `results_selected.csv`

### Selection from existing results

```bash
python bisearch_epcr.py select \
    --input results.csv \
    --output selected.csv \
    --max-per-celltype 20
```

### Validation

```bash
python validate_epcr.py --genome-dir hg19_bs --chrom chr22
```

Compares hash-based results against brute-force genome scan.

## Performance

Benchmarked on chr22 (51M bp) with 100 synthetic primer pairs:

| Metric | Value |
|--------|-------|
| Per-primer per-chromosome | ~30 ms |
| Speedup vs original | 19.7x |
| Index building (chr22) | ~20s (4 templates) |
| Memory (chr22) | ~1.3 GB (4 indexes + 4 converted genomes) |

Estimated time for full batch (200 primers x 82 cell types x 24 chromosomes):
- Screening: ~3.3h (1 core), ~50min (4 cores)
- Index building: ~32min (one-time per run)
- Total: ~4h (1 core), ~1.5h (4 cores)

Repetitive primers (>200 verified hits on one template) are immediately
marked FAIL, avoiding the O(N^2) product pairing explosion.

## Pipeline Integration

The hash-based engine is the specificity screening module (Step 4) in the
cell-specific epigenetic dPCR assay design pipeline:

- `methyl_panel/bisearch_epcr_core.py`: Core engine (`screen_pair()`, `run_phase4_epcr()`)
- `methyl_panel/prepare_genomes.py`: Genome preparation (`.npy` file generation)

To enable in the pipeline:
1. Run `python methyl_panel/prepare_genomes.py --output-dir data/hg19_bs`
2. Set `--epcr-genome-dir data/hg19_bs/` on the pipeline CLI (or the
   `EPCR_GENOME_DIR` env var)

## Files

| File | Description |
|------|-------------|
| `bisearch_epcr.py` | Standalone CLI tool (screen + select) |
| `prepare_genomes.py` | Download and encode hg19 as .npy files |
| `validate_epcr.py` | Brute-force correctness validation |
| `requirements.txt` | Python dependencies |

## References

- Tusnady et al. 2005. BiSearch: ePCR tool for bisulfite PCR. NAR 33:e9
- Aranyi et al. 2006. BiSearch ePCR application. BMC Bioinformatics 7:431
- Aranyi & Tusnady 2007. BiSearch ePCR. Methods Mol Biol 402:385
