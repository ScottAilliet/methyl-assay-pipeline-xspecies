# Cell-Specific Epigenetic dPCR Assay Design Pipeline

A complete pipeline for designing bisulfite-PCR primers for cell-type-specific
methylation dPCR assays, based on the Loyfer et al. (Nature 2023) human
methylation atlas (GSE186458).

## Features

- **DMR discovery**: wgbstools `find_markers` discovers DMRs per cell type from 207 atlas samples
- **Atlas background**: All 82 cell types serve as the default background for DMR discovery and per-CpG DMS scoring; WBC-in-blood assays automatically restrict the background to other blood cell types (see Cell Types)
- **Per-CpG resolution**: Individual CpG labels, positions, and beta values for all 207 samples
- **DMS scoring**: Per-CpG DMS score (>= 0.6) filters poor-quality DMR blocks
- **Primer3Plus-compatible**: All 161 Primer3Plus parameters configurable
- **Hash-based ePCR specificity**: BiSearch 16-mer hash screening against bisulfite-converted genome arrays, no external aligner dependency (Step 4)
- **Secondary structure screening**: Hairpin MFE with -1.5 kcal/mol cutoff
- **Common SNP screening**: dbSNP-based screening with 3' end priority
- **Repeat and segdup annotation**: UCSC track-based screening for repeats and segmental duplications
- **DimerDetective**: Primer-dimer prediction with -1.0 kcal/mol conservative cutoff
- **U-assays output format**: XLSX + PDF (one page per primer pair)
- **Modular execution**: Run any combination of pipeline steps
- **macOS compatible**: Runs locally on Apple Silicon (M1/M2/M3/M4) or Intel

## Quick Start

### One-time setup

```bash
cd ~/cell-specific-epigenetic-dpcr-assay-design-pipeline
git pull origin main
chmod +x setup_local.sh
./setup_local.sh
```

This installs Homebrew deps, Python packages, wgbstools,
downloads ~10 GB of data, and decompresses the genome.
Takes ~2-3 hours. Idempotent — safe to re-run.

See `docs/handbook.md` for the full step-by-step handbook.

### Run the pipeline

```bash
source .venv/bin/activate

# Single cell type with DMR discovery
python -m methyl_panel.pipeline --steps all --discover-dmrs \
    --cell-type MONO \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/MONO/

# Multiple cell types
for ct in MONO BCELL NK GRAN CD3T CD4T CD8T; do
    python -m methyl_panel.pipeline --steps all --discover-dmrs \
        --cell-type $ct \
        --genome data/hg19/hg19.fa.gz \
        --min-tm 58 --opt-tm 60 --max-tm 62 \
        --output-dir results/$ct/
done
```

### Post-pipeline annotation

```bash
chmod +x run_full_pipeline.sh
./run_full_pipeline.sh
```

Runs DMS sub-score computation, UCSC annotation, per-CpG beta extraction,
result merging, and data dictionary generation.

## Pipeline Steps

| Step | ID | Module | Description |
|------|----|--------|-------------|
| 0 | 0 | `phase0_dmr_discovery` | DMR discovery via wgbstools find_markers |
| 1 | 1 | `phase1_dmr_loader` | Load DMR blocks from Excel (legacy alternative to Step 0) |
| 2 | 2 | `phase2_bisulfite_convert` | Bisulfite convert to 6 strands |
| 3 | 3 | `phase3_primer3_design` | Design primers with Primer3 |
| 4 | 4 | `bisearch_epcr_core` | Hash-based ePCR specificity screening |
| 5 | 5 | `phase5_structure` | Secondary structure screening |
| 6 | 6 | `phase6_snp` | Common SNP screening (dbSNP) |
| 6b | 61 | `phase8_repeats` | Repeat and segmental duplication annotation |
| 7 | 7 | `phase7_dimer` | DimerDetective primer-dimer prediction |
| 8 | 8 | `output_xlsx` | XLSX output (~79 columns) |
| 9 | 9 | `output_pdf` | PDF output |

**Default step sequence**: `0, 2, 3, 4, 5, 6, 61, 7, 8, 9` (with `--discover-dmrs`)
or `1, 2, 3, 4, 5, 6, 61, 7, 8, 9` (legacy, loads from Excel via Step 1).

Steps 4 and 6 require downloaded data (.npy genome arrays, dbSNP VCF).
If data is missing, those steps are skipped gracefully (NULL values in output).

## ePCR Screening (Step 4)

Hash-based electronic PCR using the BiSearch 16-mer algorithm against
bisulfite-converted hg19 genome arrays (`.npy` format). No external aligner
or samtools dependency.

- **Input**: Primer sequences from Step 3
- **Output**: `ePCR_verdict` (PASS/FAIL), `ePCR_product_count`, `ePCR_products`, `ePCR_error`, `ePCR_elapsed_sec`
- **Verdict logic**: PASS = single product on intended methylation state; FAIL = no product, wrong state, or multiple products
- **Genome arrays**: Use `prepare_genomes.py` to generate `.npy` files from FASTA
- **CLI**: `--epcr-genome-dir data/hg19_bs/` (default)

See `docs/bisearch_epcr_README.md` for algorithm details.

## Pipeline Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Tm range | 55-61 C (opt 58) | Roche DLC conditions |
| Salt preset | roche_dlc | 100mM Na, 4.5mM Mg |
| Product size | 40-150 bp | dPCR amplicon range |
| Min CpG/block | 5 | Enough CpGs for discrimination |
| CpG pair total | 3 | Relaxed for sparse-CpG DMRs |
| CpG per primer | 1 | Relaxed for sparse-CpG DMRs |
| Direction | both | Hypo + hyper DMRs |
| Top markers | 1000 | Per find_markers |
| Delta means | 0.3 | Min methylation difference |
| Max bg samples | 30 | Subsample background for speed |
| DMS score | >= 0.6 | Per-CpG quality filter |
| ePCR mismatch profile | 1233333333333333 | 3' strict, 5' permissive |

## Cell Types

All 82 cell types are discovered from scratch via `--discover-dmrs` using
wgbstools `find_markers` against the Loyfer et al. atlas beta files
(GSE186458, hg19). The full list is defined in `data/sample_groups.csv`,
which maps 207 atlas samples into 82 target groups.

The 82 types span all major tissues — liver, pancreas, colon, lung, heart,
brain, kidney, breast, prostate — and 14 blood/immune subtypes. Seven
convenience aliases merge related blood subtypes into broader assay targets:

| Alias | Atlas groups merged |
|-------|-------------------|
| MONO | Blood-Monocytes |
| BCELL | Blood-B, Blood-B-Mem |
| NK | Blood-NK |
| GRAN | Blood-Granulocytes |
| CD3T | All 9 T-cell subtypes |
| CD4T | 4 CD4+ T-cell subtypes |
| CD8T | 4 CD8+ T-cell subtypes |

22 of the 82 types have only 1 sample. Results for these should be
interpreted cautiously — there is no within-target variability estimate.

### Background selection

`wgbstools find_markers` identifies regions differentially methylated
between a target cell type and a background set. The default background
is all 82 cell types, which is appropriate for counting tissue-derived
cells in blood (where WBC methylation is the natural background).

To count WBC subtypes in blood, the background is automatically restricted
to other blood cell types (36 samples, 14 `Blood-*` groups) — because blood
is composed of WBCs, the relevant discrimination is between blood subtypes,
not against liver or brain tissue.

## Configuration

The pipeline uses `Primer3PlusConfig` (in `config.py`) which supports all 161
Primer3Plus parameters. Load settings from a Primer3Plus settings file:

```python
from methyl_panel.config import Primer3PlusConfig
config = Primer3PlusConfig.from_primer3plus_file("Primer3Plus_Settings.txt")
```

Or configure programmatically:

```python
config = Primer3PlusConfig()
config.primer_min_tm = 55.0
config.primer_opt_tm = 58.0
config.primer_max_tm = 61.0
config.validate_tm()  # Raises if Tm not set
```

### Salt Presets

- **primer3plus**: 50mM Na, 3mM Mg, 1.2mM dNTP, 250nM primer (SantaLucia correction)
- **roche_dlc**: 100mM Na, 4.5mM Mg, 1.2mM dNTP, 250nM primer (Owczarzy correction)

## Data Source

Loyfer et al., "A DNA methylation atlas of normal human cell types",
Nature 613, 355-364 (2023). GEO accession: GSE186458.

## License

MIT
