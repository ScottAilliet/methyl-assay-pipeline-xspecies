# Methylation Databases with Cell-Type-Specific CpG Information: Landscape Survey and Integration Ranking

Prepared for the Unified Epigenetic Analysis Pipeline / cell-specific dPCR assay database project.
Date: 2026-10-09. Scope: every public methylation database or resource carrying cell-type-specific (or tissue-specific) CpG-level information, any organism, with a feasibility verdict for feeding the per-CpG cleanliness-scoring assay-design pipeline.

---

## 1. Executive summary

- The human cell-type methylome space is dominated by one pipeline-grade resource (Loyfer et al. 2023, the project's current source) and a small number of complementary WGBS resources — chiefly **BLUEPRINT** (~200 immune-cell WGBS methylomes, GRCh38) and **ENCODE/Roadmap** — plus several array-based cell-type reference panels that are useful for cross-validation but not for primer design.
- **No existing database provides ready-to-use cell-type-specific dPCR assays.** Every resource stores methylation data, not assay designs. The project's pipeline fills a genuine gap; the improvement path is ingesting additional WGBS sources, not replacing the approach.
- Organism coverage of public cell-specific CpG methylation data: **human and mouse are the richest**; **Arabidopsis** has the best plant coverage (6 FACS-sorted root cell types + gametophytic cell types); **zebrafish** has sorted-population data (germline vs. soma); **pig** has 9 sorted immune-cell WGBS methylomes; **rice/maize** have gamete/seed-cell methylomes only. **Rat has none** — only tissue/brain-region WGBS and scATAC — which is why rat was dropped from the build plan.
- Top integration recommendations: (1) lift the existing panel to GRCh38 and cross-validate against BLUEPRINT; (2) add BLUEPRINT immune cell types absent from Loyfer (e.g., eosinophils, basophils, macrophage polarization states); (3) use EWAS Data Hub / IDOL array references as orthogonal validation of marker specificity, not as design input.

---

## 2. Human resources

### 2.1 WGBS, per-CpG, cell-type-resolved (pipeline-grade)

| Resource | Cell-type coverage | Build | Format | Pipeline verdict |
|---|---|---|---|---|
| **Loyfer et al. 2023** (GSE186458) — current source | 39 cell-type groups, 207 deep WGBS methylomes (30×), FACS-sorted from healthy tissues | hg19 | per-CpG beta files + methylation blocks | **Primary source.** Best depth + block structure matches pipeline design. |
| **BLUEPRINT** (blueprint-epigenome.eu; data via EBI FTP/EGA) | >200 WGBS maps of sorted immune cells: monocytes, neutrophils, eosinophils, macrophages (M0/M1/M2), naive CD4+/CD8+ T cells, B cells, progenitors, plus leukaemia samples | GRCh38 | per-CpG methylation bigWig + hypo/hyper-methylated region calls | **Highest-value addition.** Adds immune cell types Loyfer lacks (eosinophils, basophils, M1/M2 macrophages). Requires hg38 handling (see §5 liftover). |
| **ENCODE** WGBS (encodeproject.org) | Primary cells (B cell, CD14+ monocyte, CD4+ memory/naive T, endothelial, epithelial, fibroblast, stem cells) + cell lines | hg38 | per-CpG bigWig/bedGraph | Complementary; smaller depth than Loyfer/BLUEPRINT for most types. |
| **Roadmap Epigenomics** (egg2.wustl.edu/roadmap; Human Epigenome Atlas) | 111 consolidated epigenomes: tissues, primary cells, stem cells | hg19 | WGBS (bis-seq) + 450k arrays | Mostly tissue-level; useful as background methylation reference. |
| **IHEC Epigenome Portal** (epigenomesportal.ca) | Aggregates BLUEPRINT + Roadmap + ENCODE + Canadian/Japanese epigenomes | mixed | standardized metadata + links | Meta-catalogue; use for discovery, not as a data source per se. |
| **MethBase** (smithlabresearch.org; UCSC track hub) | Hundreds of public BS-seq methylomes across organisms incl. human cell types/tissues | multiple | per-CpG `.meth.bw` + coverage `.read.bw`, HMR/PMD/AMR bigBed | Usable per-CpG source; depth varies per sample; good for cross-checking DMRs. |
| **MethBank 4.0** (ngdc.cncb.ac.cn/methbank) | 1,449 single-base methylomes, 23 species, 236 tissues/cell lines; normal human cell-type datasets incorporated 2025-03 | multiple | WIG/BED per-CpG + DMR toolkit | Broad but shallow curation per sample; useful for organism coverage (incl. Arabidopsis projects added 2025-08). |
| **scMethBank** (ngdc.cncb.ac.cn/methbank/scm) | 8,328 single-cell WGBS samples, human + mouse, 29 cell types | hg19/mm10 | per-cell per-CpG | Single-cell sparsity makes it unsuitable for direct primer design; useful for marker validation in rare types. |
| **Human brain snmC atlases** (BICCN; Tian et al. 2023 et al.) | Single-cell methylomes of human brain cell types | hg38 | allc/bigwig | Neuronal mCH context complicates CpG-only design; candidate for a future brain panel. |

### 2.2 Array-based cell-type references (validation, not design input)

| Resource | Content | Verdict |
|---|---|---|
| **MethyMark** (fame.edbc.org/methymark) | Cell-type-specific methylation marks in 42 Roadmap cell types/tissues (MethyMarks, Super/Large-MethyMarks) | Array-era marks; useful orthogonal check of marker specificity. |
| **EWAS Data Hub** (ngdc.cncb.ac.cn/ewas/datahub) | 75,344 normalized array samples; reference profiles for 81 tissues/cell types incl. 25 blood cell types, 6 ancestries, 67 diseases | Best array-based reference for blood cell types; batch-normalized. Validation use. |
| **IDOL / FlowSorted.Blood.EPIC + FlowSorted.BloodExtended.EPIC** (Salas et al. 2018, 2022; GSE167998) | 12 sorted leukocyte subtypes (neutrophils, eosinophils, basophils, monocytes, naive/memory B, naive/memory CD4, naive/memory CD8, Treg, NK) + IDOL-optimized CpG libraries | Directly comparable to the project's blood panel; use to confirm dPCR markers match array-derived deconvolution markers. |
| **UniLIFE / EpiDISH** (Guo et al. 2025) | 19 immune cell types, blood, birth-to-old-age reference panel | Extends validation across age; array-based. |
| **Luo et al. 2023** (Genome Med) 12-immune-cell reference | DNAm reference matrix validated against WGBS of sorted cells (IHEC) | Cross-platform validation set. |

### 2.3 Disease-focused human methylation databases (context, not cell-type-specific)

- **DiseaseMeth 3.0** — 162 diseases, 49,949 samples, gene-centric DMGs; not per-CpG cell-type resolved.
- **MethyCancer** — cancer methylation + gene/mutation integration (2007-era).
- **EWAS Atlas** — 1,090,441 trait–methylation associations, 934 traits; knowledgebase of EWAS findings.
- **MethHC 2.0**, **PubMeth**, **MethDB**, **MethPrimerDB** (primers, not data), **MeInfoText** — legacy resources.
- **MethMotif** — cell-specific TF binding motifs with CpG methylation; niche use for motif-level annotation.
- **NGSmethDB** — single-cytosine methylomes + differential cytosines + methylation segments across species; uniform processing; track-hub access.

---

## 3. Organism coverage: which organisms have public cell-specific methylation data?

| Organism | Cell-specific methylation data | Best resource | Pipeline-grade? |
|---|---|---|---|
| **Human** | Extensive (39+ sorted cell types, deep WGBS) | Loyfer 2023; BLUEPRINT | Yes |
| **Mouse** | Extensive: brain (161–4,673 clusters), sorted neurons, immune populations | CEMBA/snmC-seq2 (Liu et al. 2021, GSE132489); snmC-seq3 (Liu et al. 2023); INTACT neurons (GSE191211); HSC (GSE95697) | Yes (sorted bulk); brain atlas requires pseudo-bulking per-cell allc (225 GB at GEO) |
| **Arabidopsis** | 6 FACS-sorted root meristem cell types + sperm/vegetative nucleus/central cell/egg/microspore | Kawakatsu et al. 2016 (GSE79708); 1001 Epigenomes (accessions, not cell types); AraENCODE (tissues) | Yes — allc-style per-cytosine TSVs, all three contexts (CG/CHG/CHH) |
| **Zebrafish** | Germline (PGC) vs. somatic cells at 4 embryonic stages × 2 reps | Norris et al. 2019 (GSE122722) — per-CpG bigwigs | Yes (processed per-CpG bigwigs, danRer10) |
| **Pig** | 9 sorted immune cell populations (B subsets, 4 T fractions, NK, myeloid, neutrophils); 54,391 cDMRs | Corbett et al. 2022, BMC Genomics | Candidate future organism |
| **Rice** | Central cell, egg cell methylomes only | Anderson et al. (plant gamete methylomes) | Too few cell types |
| **Maize** | Embryo vs. endosperm (seed tissues); UMR maps | Bai et al.; Crisp et al. 2020 | Tissue-level only |
| **Rat** | **None.** Tissue/brain-region WGBS/RRBS and scATAC only | — | **No — excluded from build** |
| **Chimpanzee** | Neuronal vs. non-neuronal RRBS (cortex) | Camp et al. 2018 | RRBS-limited; 2 cell classes |
| **Rabbit/dog/cow/chicken** | Fibroblast WGBS comparisons only | Zhuo et al. 2022 | Not cell-panel grade |

Note: no single-cell DNA methylome has been reported in plants (Nobori 2025 review); plant cell-specific methylation comes from FACS-sorted populations.

---

## 4. Integration ranking for the human pipeline

Ranked by (a) data-format fit (per-CpG beta values), (b) cell types added vs. Loyfer 2023, (c) processing effort.

1. **BLUEPRINT WGBS (GRCh38)** — adds eosinophils, basophils, M1/M2 macrophages, progenitor subsets; deep WGBS; same assay-design logic applies. Requires the hg38 panel (§5). Effort: moderate (download bigWigs, build beta matrices).
2. **hg19→hg38 liftover of the existing panel** — unlocks all GRCh38-native resources and modern annotation. Effort: low (done in this work, see §5).
3. **IDOL/FlowSorted.BloodExtended.EPIC + EWAS Data Hub references** — orthogonal validation of the blood panel's marker specificity on an independent platform (array vs. WGBS vs. dPCR). Effort: low.
4. **ENCODE WGBS** — a few additional primary-cell types (e.g., endothelial, keratinocyte, mesenchymal stem cell) with shallower depth. Effort: moderate.
5. **MethBase human methylomes** — cross-check of DMR boundaries and background methylation. Effort: low (per-CpG bigwigs).
6. **scMethBank** — validation of rare-cell markers at single-cell level. Effort: low, validation only.
7. **Roadmap/MethyMark** — legacy array marks; validation only. Effort: low.
8. **DiseaseMeth/EWAS Atlas/MethyCancer** — context for marker-disease associations in the manuscript; not design input.

---

## 5. hg19 → GRCh38 liftover (performed)

- DMR blocks and assay coordinates from the existing hg19 master panel were lifted with the UCSC hg19→hg38 chain; every lifted amplicon's sequence was re-extracted from GRCh38 and both primer binding sites re-verified in bisulfite-converted space.
- Assays failing re-verification are flagged (`liftover_status`) and excluded from the hg38 panel while remaining valid in hg19. **Both hg19 and hg38 panels are primary deliverables** — hg19 for continuity with the Loyfer source data and prior QC, hg38 for BLUEPRINT/ENCODE integration.
- Outputs: `master_assay_panel_hg38.xlsx` (+ parquet), liftover scripts in the new repo.

---

## 6. Cross-species builds (this work)

Pilot-free strategy per approved plan: acquire all processed per-cytosine data → audit gate → full runs on qualifying cell types.

| Organism | Source | Data | Status |
|---|---|---|---|
| Mouse (mm10) | GSE191211 INTACT neurons | 10 bulk WGBS `.cov.gz` (GABAergic, Glutamatergic, Purkinje × 2 reps WT + MeCP2 KO) | Downloaded; audit + run |
| Mouse (mm10) | CEMBA GSE132489/GSE157886 | per-cell allc only (225 GB); no cluster-level per-CpG release | **Audit-gate exclusion: cross-reference only** (pseudo-bulking deferred) |
| Arabidopsis (TAIR10) | GSE79708 | 9 allc-style TSVs (6 cell types + controls), all contexts | Downloaded; context-aware run |
| Zebrafish (danRer10) | GSE122722 | 16 per-CpG bigwigs (PGC vs. Soma × 4 stages × 2 reps) | Downloaded; run |

Rat: excluded (no cell-type methylome data exists — see §3).

---

## 7. Key references

- Loyfer N et al. A DNA methylation atlas of normal human cell types. Nature 613, 355–364 (2023). PMID 36626871
- Schuyler RP et al. Distinct Trends of DNA Methylation Patterning in the Innate and Adaptive Immune Systems. Cell Rep 17, 2101–2111 (2016). BLUEPRINT WGBS resource.
- Salas LA et al. Enhanced cell deconvolution of peripheral blood using DNA methylation for high-resolution immune profiling. Nat Commun 13, 761 (2022).
- Luo W et al. A meta-analysis of immune-cell fractions at high resolution. Genome Med 15, 61 (2023).
- Guo X et al. Unified high-resolution immune cell fraction estimation in blood tissue from birth to old age (UniLIFE). Nat Commun (2025).
- Liu H, Zhou J, Tian W et al. DNA methylation atlas of the mouse brain at single-cell resolution. Nature 598, 120–128 (2021).
- Liu H et al. Single-cell DNA Methylome and 3D Multi-omic Atlas of the Adult Mouse Brain. Nature 623, 120–131 (2023).
- Kawakatsu T et al. Unique cell-type-specific patterns of DNA methylation in the root meristem. Nat Plants 2, 16058 (2016).
- Kawakatsu T et al. Epigenomic Diversity in a Global Collection of Arabidopsis thaliana Accessions (1001 Epigenomes). Cell 166, 481–491 (2016).
- Wang Z et al. AraENCODE: A comprehensive epigenomic database of Arabidopsis thaliana. Mol Plant 16, 1064–1074 (2023).
- Norris WT et al. Retention of paternal DNA methylome in the developing zebrafish germline. Nat Commun 10, 3054 (2019).
- Jang HJ et al. Epigenetic dynamics shaping melanophore and iridophore cell fate in zebrafish. Nat Commun (2021).
- Corbett RJ et al. Assessment of DNA methylation in porcine immune cells reveals novel regulatory elements. BMC Genomics 23, 760 (2022).
- Zhang M et al. MethBank 4.0. Nucleic Acids Res 51, D208–D216 (2023).
- Zong W et al. scMethBank. Nucleic Acids Res 50, D212–D219 (2022).
- Schultz MD et al. MethBase: an encyclopedia of DNA methylomes. Bioinformatics (2015).
- Xiong Z et al. EWAS Data Hub. Nucleic Acids Res 48, D890–D895 (2020).
- Xing J et al. DiseaseMeth version 3.0. Nucleic Acids Res 50, D1208–D1216 (2022).
- Lebrón R et al. NGSmethDB 2017. Nucleic Acids Res 45, D97–D103 (2017).
- Xu J et al. MethyMark / FAME. (fame.edbc.org/methymark)
- Camp JG et al. Cell Type and Species-specific Patterns in Neuronal and Non-neuronal Methylomes of Human and Chimpanzee Cortices. Mol Biol Evol (2018).
- Zhuo Z et al. A comparative methylome analysis of vertebrate fibroblasts. BMC Biol 20, 71 (2022).
