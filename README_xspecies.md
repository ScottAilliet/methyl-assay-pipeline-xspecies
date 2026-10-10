# methyl-assay-pipeline-xspecies

Cross-species extension of the cell-specific epigenetic dPCR assay design
pipeline (MSP-style methylation-specific PCR primers with per-CpG
cleanliness scoring). The human pipeline (GSE186458, Loyfer et al. 2023,
18 immune cell types, hg19 + hg38 liftover) lives in the parent project;
this repo adds generic DMR discovery and context-aware assay design for
**Mus musculus**, **Arabidopsis thaliana**, and **Danio rerio**.

## What was added

| Module | Purpose |
|---|---|
| `methyl_panel/phase0_generic.py` | Organism-generic DMR discovery: loaders for mouse Bismark cov.gz, arabidopsis Kawakatsu mc_calls TSV (CG/CHG/CHH, strand-aware), zebrafish score-only bigWigs; per-chromosome beta matrices; block clustering; dms scoring; Block_Summary + PerCpG workbook output. |
| `methyl_panel/plant_design.py` | Context-aware (CG/CHG/CHH) MSP assay design for plants: discriminating sites in ANY context feed the design templates and the ePCR target-state mask (template-consistent mode). |
| `methyl_panel/phase3_primer3_design.py` | Patched: `strands_override` parameter so plant context-aware strands bypass the CpG-only conversion. |
| `methyl_panel/pipeline.py` | Patched step 4: batched 16-mer index ePCR screening (builds CT/CT_M indexes once per chromosome/strand instead of once per primer — ~N_primers-fold faster, byte-identical verdicts). |
| `sample_sheets/` | Ready-to-adapt samples.tsv files for the three organisms (GSE191211, GSE79708, GSE122722). |
| `docs/` | Survey report and database catalog from the methylation-database landscape review. |
| `xspecies_results/` | Delivered artifacts from the original session: per-organism assay workbooks, funnel summary, QC summary, audit-gate table. |

## Design/ePCR semantics (plants)

- Discriminating sites (any context) are kept as C in the methylated
  template K_M; the unmethylated template K_U keeps only CpG C's.
- ePCR target state (T) = CpG kept + ALL DMR-block discriminating sites
  kept as C — exactly matches the design template (the plant analog of the
  mammalian fully-methylated CT_M state). Background state (B) = CpG kept +
  WGBS-informed high-methylation CH sites; U = fully converted.
- PASS = exactly one amplified locus whose best target-state hit has
  strictly fewer mismatches than any background/converted hit at that
  locus. Products are deduplicated by locus across methylation states and
  reported in forward (Watson) genome space.
- Design post-filter `--min-disc-in-primers` (default 1): pairs whose
  primer spans contain no discriminating site cannot discriminate and are
  dropped.

## Mammals

Mouse and zebrafish use the mammalian CpG-default pipeline
(`python -m methyl_panel.pipeline`) with the generic phase-0 discovery
front-end; CH contexts can be enabled via phase 0's `--contexts` flag
(mammals default to CG).

## Provenance note (important)

The original cross-species modules were written and validated in a sandbox
session whose working disk was lost. They were **reconstructed** from the
session fix log, surviving file dumps, and the delivered result artifacts,
then re-validated end-to-end on arabidopsis (GSE79708 root cell types):
the full chain (discovery -> design -> ePCR) runs and produces own-locus,
target-specific ePCR-PASS assays. Rerun block/assay sets are
equivalent-quality but **not identical** to the original session's numbers
(e.g. Columella_hyper: 200 capped blocks with 1192 PerCpG sites vs 1083
originally; 6 primer pairs / 1 PASS vs 8 / 2). The authoritative result
artifacts from the original session are preserved under `xspecies_results/`
and were produced by the pre-loss code.

## Reproducing the arabidopsis validation

```bash
# 1. TAIR10 chromosomes (Ensembl Genomes plants release 57) -> .npy
python prepare_genomes.py --output-dir TAIR10_bs \
    --fasta-dir TAIR10_fa --chromosomes 1 2 3 4 5 Mt Pt

# 2. Discovery (CG+CHG+CHH)
python -m methyl_panel.phase0_generic --organism arabidopsis \
    --samples sample_sheets/arabidopsis_GSE79708_samples.tsv \
    --out dmr_blocks.xlsx --contexts CG,CHG,CHH --direction both

# 3. Context-aware design + ePCR
python -m methyl_panel.plant_design design \
    --dmr-xlsx dmr_blocks.xlsx --cell-type Columella_hyper \
    --fasta-dir TAIR10_fa --output-dir run/Columella_hyper
python -m methyl_panel.plant_design epcr \
    --primers run/Columella_hyper/primers.json --genome-dir TAIR10_bs \
    --samples sample_sheets/arabidopsis_GSE79708_samples.tsv \
    --cell-type Columella --dmr-xlsx dmr_blocks.xlsx \
    --cell-type-id Columella_hyper

# 4. Structure / dimer / report
python -m methyl_panel.pipeline --steps 5,7,8 \
    --output-dir run/Columella_hyper --dmr-xlsx dmr_blocks.xlsx \
    --cell-type Columella_hyper
```

## Data sources

- Human: GSE186458 (Loyfer et al., Nature 2023) — parent project
- Mouse brain atlas WGBS: GSE191211 (WT replicates; KO replicates excluded)
- Arabidopsis root cell types: GSE79708 (Kawakatsu et al. 2016 mc_calls)
- Zebrafish PGC/Soma: GSE122722 (score-only methylation bigWigs)

## Known limitations

- SNP and repeat screens (steps 6/6b) were skipped for non-human organisms.
- Zebrafish bigWigs are score-only: no per-site coverage filtering is
  possible; methylation-state thresholds apply to scores.
- Mouse WGBS is low-coverage (mean ~4-8x): use the >=3x filter, merged
  replicates, and relaxed min-methylation thresholds (CONDITIONAL audit
  decision).
- The human panel PerCpG annotation sheets in the parent deliverables
  inherit an annotation-only coordinate caveat (block coordinates and
  primer sequences are sequence-derived and correct).
