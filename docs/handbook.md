# Cell-Specific Epigenetic dPCR Assay Design Pipeline — Operator Handbook

**Version:** Pipeline v3.2.0
**Repository:** `cell-specific-epigenetic-dpcr-assay-design-pipeline`
**Last updated:** 2026-08-25

---

## PREFACE

### What This Handbook Is

This handbook is the complete guide for running the Cell-Specific Epigenetic dPCR Assay Design Pipeline primer design pipeline on your Mac. It covers installation, data download, running the pipeline, and understanding the output — written for a lab scientist who is not a programmer.

### What Changed in v3.2.0

This version removes bowtie entirely, making hash-based ePCR the sole specificity screening method (Step 4), adds ePCR verdict columns to the Excel/PDF output, and streamlines the repo by removing obsolete files.

Key changes:

1. **Removed: bowtie v1 specificity screening.** The bowtie-based Step 4 (`phase4_bowtie_specificity.py`) and all bowtie-related infrastructure (index builders, conda env, config fields, output columns) have been removed. Hash-based ePCR is now the sole specificity screening method.

2. **Renumbered: ePCR is now Step 4.** The BiSearch 16-mer hash ePCR engine (`bisearch_epcr_core.py`) is now Step 4 (step ID 4), replacing the old bowtie step. For each primer pair, it searches bisulfite-converted hg19 genome arrays (`.npy` format) for PCR products and produces a template-aware verdict: PASS (single product on intended methylation state) or FAIL (no product, wrong methylation state, or multiple products). Uses `--epcr-genome-dir` (default: `data/hg19_bs/`). If the genome directory is not found, ePCR fields are set to NULL and screening is skipped (graceful degradation). Default step sequence is now `[0, 2, 3, 4, 5, 6, 61, 7, 8, 9]`.

3. **New: ePCR columns in Excel/PDF output.** Five columns in `primer_assays.xlsx`: `ePCR_verdict` (PASS/FAIL), `ePCR_product_count`, `ePCR_products` (formatted coordinates), `ePCR_error`, `ePCR_elapsed_sec`.

4. **Streamlined: removed obsolete files.** Old handbook versions, manuscript files, old runner scripts, legacy scripts, bowtie-related modules and indices, Flask web app, and large data xlsx files. The repo now contains only the current pipeline, ePCR tools, current handbook, tests, and essential docs.

5. **Renamed: CSV files.** `full_atlas_groups.csv` → `sample_groups.csv`, `full_atlas_manifest.csv` → `sample_manifest.csv`. Scripts `generate_master_excel.py` and `export_dmr_excel.py` moved to `scripts/`.

### What Changed in v3.1.0

This version adds a standalone ePCR screening workflow, post-pipeline annotation scripts, a consolidated runner, and renames "cleanliness score" to "DMS score" (Differential Methylation Suitability).

Key changes:

1. **Renamed: "cleanliness score" → "DMS score" (Differential Methylation Suitability).** The per-CpG quality score that quantifies how suitable a candidate DMR is for bisulfite PCR primer design is now called the DMS score throughout the handbook and output files. The source code function is still named `cleanliness_score()` for backward compatibility. The formula is unchanged: `score = mean(A, B, C, D) * E` where A = 1 - max(target CpG mean methylation), B = min(background CpG mean methylation), C = max(0, 1 - 2*mean(target CpG std)), D = max(0, 1 - 2*mean(background CpG std)), E = mean(fraction CpGs with adequate coverage in target and background).

2. **New: standalone ePCR screening tool (`bisearch_epcr.py`).** A standalone bisulfite ePCR screening tool at the repository root, distinct from the pipeline's internal Phase 4. Uses hash-based numpy indexing (no bowtie2 or samtools dependency), searches both methylation states, and produces a 3-way template-aware verdict: PASS / FAIL wrong_genome / FAIL not_methylation_specific. See Chapter 6b for full documentation.

3. **New: post-pipeline annotation scripts.** Five scripts that run after the main pipeline to compute DMS sub-scores, annotate primers with UCSC reference tracks, extract per-CpG beta values, merge results, and generate a data dictionary. Orchestrated by `run_full_pipeline.sh`. See Chapter 6b.

4. **Changed: selection logic keeps ALL passing primers.** The ePCR selection step no longer deduplicates by DMR or caps per cell type. All primers that pass the 3-way verdict are kept. Previously, only one primer per DMR and a maximum of 20 per cell type were retained.

5. **Column changes in final atlas.** Removed `direction` (unreliable), `cleanliness_score` (renamed to `dms_score`), and 8 bowtie2-specific columns. Added ePCR verdict columns, DMS sub-score columns, and UCSC annotation columns. See Chapter 6b.9 for the full column list.

### What Changed in v2.2.7

This version fixes five data quality issues found during analysis of the v2.2.6 results.

Key changes:

1. **Fixed: blocks with empty or partial bg_subgroup_meth bypassed the per-subgroup filter.** The filter code had `if not b.bg_subgroup_meth: filtered.append(b); continue` — blocks with no subgroup data passed automatically. And blocks with partial data (e.g. 4 of 13 subgroups present) were checked only on the present subgroups; the 9 missing ones were silently ignored. Now blocks with empty or partial bg_subgroup_meth are **rejected**. Missing subgroup data means we cannot verify that subgroup is sufficiently methylated — the safe choice is to reject. The rejection report now categorizes rejections as "no subgroup data", "partial subgroup data", or threshold violations.

2. **Added: cleanliness score filter (CS >= 0.6).** Blocks with CS < 0.6 (low coverage, high variability, or no usable CpG data) were previously computed but never filtered out. 6 blocks in v2.2.6 had CS < 0.6 (3 NK with CS=0, 1 GRAN with CS=0.34, 1 MONO with CS=0.53, 1 CD3T hyper with CS=0.37). Now they are rejected at Step 0e1, before the per-subgroup filter runs.

3. **Added: sex chromosome exclusion (chrX, chrY, chrM).** 2 blocks in v2.2.6 were on chrX (KDM6A, SEPT6). chrX has copy number differences between sexes (females 2 copies, males 1) and X-inactivation methylation differences, making dPCR quantification unreliable. chrY is male-only. chrM has different copy number biology. All three are now excluded at Step 0d2, before per-CpG extraction.

4. **Fixed: hyper sheet naming in Excel export.** When hypo and hyper markers for the same cell type were exported to one workbook, the hyper sheets got names like `BCELL1`, `CD3T1` — not obviously "hypermethylated." Now sheets are named `BCELL_hypo` and `BCELL_hyper` (or `MONO_hypo`, `NK_hyper`, etc.).

5. **Added: Direction column to all Excel sheets.** Block-level sheets, the summary sheet, and the Per-CpG detail sheet now include a Direction column showing `U` (hypomethylated) or `M` (hypermethylated), so you can always tell which type of marker a row represents.

### What Changed in v2.2.6

This version corrects the DMR discovery strategy and widens the candidate pool.

Key changes:

1. **Top 1000 markers (was 300):** `--top-markers` default changed from 300 to 1000 in `pipeline.py`, `phase0_dmr_discovery.py`, and `config.py`. find_markers now returns the top 1000 DMRs per cell type, giving the strict 0.70 per-subgroup filter 3x more candidates to work with. With 300 candidates, the filter left only 9 markers for MONO and 6 for NK — too few for a panel. With 1000, more correctly-ranked candidates enter the filter.

2. **Hybrid atlas approach deprecated for hypo markers.** The `--use-full-atlas` flag (added in v2.2.5) is no longer recommended for hypomethylated marker discovery. Testing showed it finds more raw candidates (11,214 for MONO vs 3,804 for NK) but the **ranking** is wrong for blood separation: find_markers sorts by delta against the average of all 204 atlas background samples (liver, brain, colon, granulocytes mixed together). The top 300 by atlas delta are the regions most different from non-blood tissues — NOT the regions most different from blood cells. A region where MONO=0.01, granulocytes=0.55, liver=0.95 gets a huge atlas delta (0.84) and ranks high, but granulocytes at 0.55 fails the 0.70 blood filter. Blood-only find_markers ranks by blood delta instead — the correct metric for what we're filtering on. The `--use-full-atlas` flag remains available for experimentation but is not in the recommended commands. Hypermethylated markers (`--only-hyper`) remain valuable — they search a completely different set of regions.

### What Changed in v2.2.5

This version introduces three improvements to DMR marker discovery: a **hybrid atlas approach** for more candidate markers, **hypermethylated marker support**, and a **strict per-subgroup filter** at 0.70.

Key changes:

1. **Hybrid atlas approach (`--use-full-atlas`):** By default, find_markers uses blood-only background (36 samples) to discover DMRs. This is biologically correct but finds fewer candidates for closely related T-cell subtypes (CD4T found only 17 markers, CD8T only 11). The new `--use-full-atlas` flag uses the full 207-sample atlas as find_markers background — this finds more candidate DMRs because the target cell type is more different from liver/brain/colon than from other blood cells. The per-CpG extraction and per-subgroup filter still use blood-only background, so only markers that are also well-separated from every blood cell type survive. The biological logic: markers that distinguish a cell type from all human cell types are a superset of markers that distinguish it from blood cells only — so using the full atlas for discovery and blood-only for quality control gives the best of both worlds.

2. **Hypermethylated markers (`--only-hyper`):** By default, the pipeline discovers **hypomethylated** markers — regions where the target cell type is unmethylated (< 15%) and all other blood cells are methylated (> 65%). The new `--only-hyper` flag discovers the mirror image: regions where the target is methylated and all other blood cells are unmethylated. This doubles the marker search space. The per-subgroup filter is automatically inverted for hyper markers: instead of requiring every background subgroup to be methylated (> 0.70), it requires every background subgroup to be unmethylated (< 0.30). The cleanliness score is also inverted: component A measures how close the target is to 1.0 (instead of 0.0), and component B measures how close the background is to 0.0 (instead of 1.0). The Excel export shows "Max subgroup methylation" instead of "Min subgroup methylation" for hyper markers. The `direction` field in `dmr_blocks.json` is set to `"M"` for hyper and `"U"` for hypo.

3. **Strict per-subgroup filter (0.70, was 0.50):** The `--min-bg-subgroup-meth` default is raised from 0.50 to 0.70. This means every background blood cell type must be at least 70% methylated (for hypo markers) or at most 30% methylated (for hyper markers) at a DMR region. The previous 0.50 threshold allowed markers where a background cell type was only partially methylated — adequate but not ideal for a dPCR assay where binary discrimination matters. The stricter 0.70 threshold ensures the DMR region itself provides strong binary separation, rather than relying on primer design constraints to compensate for marginal discrimination. Discrimination must happen at the DMR region level.

### What Changed in v2.2.4

This version fixes two bugs in per-CpG methylation extraction and Excel export that affected data integrity.

Key changes:

1. **Fixed: Blood-B-Mem silently dropped from MONO/NK/GRAN backgrounds.** The `max_bg_samples=30` truncation in `extract_percpg_methylation_with_indices()` used a naive `background_betas[:30]` slice. For MONO/NK/GRAN, there are 33 background samples, and the last 3 happened to be both Blood-B-Mem samples — so the entire Blood-B-Mem subgroup disappeared from per-subgroup methylation data, the per-subgroup filter, and the Excel output. Replaced with `_subsample_background()`, a subgroup-aware subsampling function that guarantees at least one sample per subgroup using the largest-remainder allocation method. All 13 blood subgroups are now preserved for every cell type.

2. **Fixed: Excel showed all 14 subgroup columns on every sheet.** `export_dmr_excel.py` collected all subgroup names across all cell types and put them on every per-cell-type sheet. For CD4T, 5 columns were empty (4 CD4 target groups + excluded Blood-T-CD3); for CD3T, 9 columns were empty. Now each per-cell-type sheet only includes columns for subgroups that actually appear in that cell type's `bg_subgroup_meth` data. The "All blocks summary" sheet still shows all subgroups (correct for cross-cell-type comparison).

3. **Audit finding: repo 1 vs repo 2 delta differences explained.** Repo 1 used the full 18-cell-type atlas (45 samples) as find_markers background but accidentally computed per-CpG cleanliness scores on a blood-only background (the `bg_betas[:30]` truncation dropped all non-blood samples, which were at the end of the list). This meant find_markers found markers with large deltas (blood vs liver/colon) while per-CpG scoring used blood-only background. Repo 2 correctly uses blood-only for both find_markers and per-CpG scoring, giving smaller but biologically correct deltas for a WBC assay.

### What Changed in v2.2.3

This version improves DMR marker quality by adding a **per-subgroup background check**, tightening the find_markers mean thresholds, and excluding overlapping cell populations from CD4T/CD8T backgrounds.

Key changes:

1. **Per-subgroup background check:** After find_markers, the pipeline now computes mean methylation for **each background blood cell type independently** (e.g. B cells, NK cells, CD8 T cells, monocytes, granulocytes). Any DMR block where a single background subgroup has mean methylation below the threshold is rejected. This catches the key problem with CD4T/CD8T: a region where CD4T is unmethylated (0.01) and the average background looks fine (0.88) but CD8 T cells are only partially methylated (0.55) — which would cause false positives in a real blood sample containing CD8 cells. The threshold can be adjusted with `--min-bg-subgroup-meth`. The default was 0.50 in v2.2.3 and raised to 0.70 in v2.2.5 — every background subgroup must now be at least 70% methylated (for hypo markers) or at most 30% methylated (for hyper markers). This ensures strong binary discrimination at the DMR region level.

2. **T-CD3 excluded from CD4T/CD8T backgrounds:** Blood-T-CD3 is a pan-T cell population containing both CD4+ and CD8+ T cells. It overlaps with the CD4T and CD8T targets, so it is now excluded from their backgrounds (via `BG_EXCLUDE` in `phase0_dmr_discovery.py`). Without this exclusion, the per-subgroup filter rejected 68/73 CD4T markers because T-CD3 was "partially unmethylated" — it contains CD4+ T cells that are unmethylated at CD4-specific regions. With the exclusion, CD4T background is B cells, NK cells, granulocytes, monocytes, and CD8 T cells only.

3. **Tighter find_markers mean thresholds:** Added `--unmeth_mean_thresh 0.15` (target mean must be below 15%) and `--meth_mean_thresh 0.65` (background mean must be above 65%). These were previously disabled (find_markers defaults of 1.0 and 0.0), allowing mediocre markers with target=0.3/bg=0.6 to pass. The `delta_means` threshold stays at 0.3 (was briefly raised to 0.4 in an earlier version of v2.2.3 but reverted because it found too few markers for CD4T/CD8T — only 73 and 19 respectively).

4. **Excel export updated:** `export_dmr_excel.py` now includes per-subgroup methylation columns (one per blood cell type), a "Min subgroup methylation" column, and a "Worst background subgroup" column, so you can see exactly which cell types are well-separated and which are borderline.

5. **Why blood-only background gives lower deltas:** Switching from the full 207-sample atlas to blood-only background (36 samples) naturally reduces delta values for all cell types. MONO vs liver/brain is a large epigenetic distance (delta ~0.95); MONO vs T cells/B cells is smaller (delta ~0.88). CD4T vs CD8T is the smallest — they are closely related T-cell lineages. **This is expected and correct for a WBC assay**: the deltas that matter are the ones against cells that will actually be in the same tube, not against liver tissue.

### What Changed in v2.2.2

This version changes the default number of DMR candidates per cell type from 200 to **300**, and adds a scipy pre-check to prevent cryptic crashes during DMR discovery.

Key changes:

1. **Top 300 markers (was 200):** `--top-markers` default changed from 200 to 300 in `pipeline.py`, `phase0_dmr_discovery.py`, and `config.py`. find_markers now returns the top 300 DMRs per cell type, giving primer3 more candidates to work with. You can still override with `--top-markers 200` or any other number.
2. **scipy pre-check:** Before running `find_markers`, the pipeline now verifies that scipy is installed in the venv Python. If it's missing (e.g. rolled back by a failed `pip install -e .`), the pipeline prints the exact install command instead of letting wgbstools crash with an opaque traceback.
3. **install_macos.sh:** Added scipy re-verification after the wgbstools install step — automatically reinstalls scipy if it was lost.

### What Changed in v2.2.1

This version restricts the DMR discovery background to **blood cells only**. Previously, find_markers used all 207 atlas samples (including liver, brain, colon, etc.) as background. For a WBC assay, this is biologically incorrect — the DMRs must distinguish the target cell type from other blood cells, not from non-blood tissues. Now only the 36 blood samples (14 blood cell types) are used as background.

Key change:

1. **Blood-only background:** `generate_groups_file()` now filters to `Blood-*` groups only. Each cell type's background is the other blood cell types (e.g. MONO background = B, NK, granulocytes, T cells = 33 samples; CD3T background = B, NK, granulocytes, monocytes = 14 samples).

### What Changed in v2.2.0

This version adds **DMR discovery from scratch** — the pipeline is now fully self-contained. You no longer need a pre-computed DMR Excel file. The new Step 0 uses `wgbstools find_markers` to discover differentially methylated regions directly from the raw 207 beta files downloaded by `download_data.sh`.

Key changes:

1. **New Step 0 (DMR Discovery):** `--discover-dmrs` flag runs `wgbstools find_markers` on raw beta files to find DMRs from scratch. No external DMR file needed.
2. **New module:** `phase0_dmr_discovery.py` — generates groups files, runs find_markers, extracts per-CpG methylation, computes cleanliness scores.
3. **`install_macos.sh` fixed:** wgbstools installed via symlink to `.venv/bin/wgbstools` (not `pip install -e .`, which fails because the package has no `__init__.py`). Added scipy dependency. Removed unnecessary `init_genome` call.
4. **Cell type target mapping:** Each of the 7 cell types maps to specific atlas groups (e.g. CD3T = all T cell subtypes, BCELL = B + B-Memory).
5. **Backward compatible:** Step 1 (load from Excel) still works for users who have a pre-computed DMR file.

### Who This Handbook Is For

You are a doctor, medical researcher, or lab scientist. You know molecular biology (PCR, DNA, methylation) but you are not a programmer. This handbook explains the computing concepts you need, in plain language with analogies to things you already know.

If you want to understand the science behind the pipeline, see [Materials and Methods](materials_and_methods.md). This handbook focuses on how to run it.

---

## Chapter 0 — Computing Basics for Lab Scientists

### 0.1 The Terminal (Command Line)

**What it is:** The Terminal is a text-based way to control your computer. Instead of clicking icons, you type commands. On Mac, open it from Applications → Utilities → Terminal (or press Cmd+Space and type "Terminal").

**The analogy:** Think of the Terminal as a very literal lab assistant. You tell it exactly what to do, step by step. If you type the wrong instruction, it does the wrong thing — but it never guesses.

**What you see when you open it:**

```
scott@MacBook-Pro ~ %
```

The `%` symbol is the **prompt** — it means "I am ready for your command." You type after it and press Enter.

**Your first commands:**

```bash
pwd
```
"Print working directory" — tells you which folder you are in. Like asking "which lab am I standing in?"

```bash
ls
```
"List" — shows files and folders in your current location. Like looking at the shelves.

```bash
ls -la
```
Options (flags) always start with a dash. `-l` means "long format", `-a` means "all files including hidden ones".

---

### 0.2 Folders, Paths, and Navigation

**Folders are called directories** in terminal language. Same thing.

**A path is the address of a file or folder.** Two kinds:

- **Absolute path:** Starts from the root of your hard drive. Always starts with `/`.
  - `/Users/scottailliet/Documents/data/` — like a full postal address
  - Your home folder: `/Users/scottailliet/` — abbreviated as `~`

- **Relative path:** Starts from wherever you are now.
  - `data/beta_files/` means "from here, go into data, then beta_files"
  - `.` means "this folder"
  - `..` means "the parent folder (one level up)"

**Moving between folders:**

```bash
cd Documents
cd ..
cd ~
cd /Users/scottailliet/cell-specific-epigenetic-dpcr-assay-design-pipeline
```

**Pro tip:** Type the first few letters of a folder name, then press **Tab** — the terminal auto-completes it. This saves time and prevents typos.

---

### 0.3 Virtual Environments (Instead of Conda)

This pipeline uses a Python **virtual environment** (`.venv`) instead of conda. A virtual environment is an isolated Python setup that keeps this project's packages separate from everything else on your Mac.

**The analogy:** Like a clean bench in the lab — this project's reagents stay on this bench and don't mix with other projects.

**Creating and using it:**

```bash
python3.11 -m venv .venv
source .venv/bin/activate
deactivate
```

**Key point:** You must activate the environment every time you open a new Terminal window and want to use the pipeline. The environment does not persist between Terminal sessions.

---

### 0.4 Installing Software

The `install_macos.sh` script handles everything automatically. It installs:
- **Homebrew** (Mac package manager) — if not already installed
- **System tools:** samtools, htslib (includes tabix), bedtools
- **Python packages:** primer3-py, openpyxl, pandas, numpy, reportlab, pysam
- **wgbstools** (from GitHub source — not on PyPI)

```bash
chmod +x install_macos.sh
./install_macos.sh
```

`chmod +x` means "make this file executable" — you only do this once. `./` means "run the file in this folder."

> **If you get "permission denied":** You forgot `chmod +x`. Run `chmod +x install_macos.sh` first, then `./install_macos.sh`.

> **If you get "command not found: #":** You pasted a `#` comment line from this handbook. Lines starting with `#` are explanations, not commands. Only paste lines that do NOT start with `#`.

---

### 0.5 Running Commands and Scripts

**Running the pipeline:**

```bash
python -m methyl_panel.pipeline --steps all \
    --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/
```

The `\` at the end of a line means "this command continues on the next line." You can also type it all on one line.

**Arguments (flags):**
- `--steps all` means "run all steps (0, 2, 3, 4, 41, 5, 6, 61, 7, 8, 9)"
- `--steps 1,2,3` means "run only steps 1, 2, and 3"
- `--min-tm 58` means "set the minimum Tm to 58°C"
- `--epcr-genome-dir data/hg19_bs/` points to the directory with `.npy` genome arrays for Step 4b ePCR screening (default: `data/hg19_bs/`). If not found, ePCR is skipped gracefully.
- `--help` shows all available arguments

```bash
python -m methyl_panel.pipeline --help
```

---

### 0.6 Files You Will Encounter

| Extension | What it is | How to open |
|-----------|-----------|-------------|
| `.py` | Python script (code) | Text editor or run with `python` |
| `.md` | Markdown document (like this handbook) | Text editor, or rendered in GitHub |
| `.xlsx` | Excel spreadsheet | Microsoft Excel, Google Sheets |
| `.bed.gz` | Compressed genomic coordinates file | `gunzip` to decompress |
| `.beta` | Methylation data file (Loyfer format) | Processed by wgbstools |
| `.csv` | Comma-separated values | Excel or text editor |
| `.json` | Data file (JavaScript Object Notation) | Text editor or Python |
| `.fa.gz` | Compressed genome FASTA | Used by samtools faidx |
| `.vcf.gz` | Compressed variant file (dbSNP) | Used by tabix |
| `.pdf` | PDF report (pipeline output) | Any PDF viewer |

---

### 0.7 Common Mistakes and How to Avoid Them

| Mistake | What happens | Fix |
|---------|-------------|-----|
| Not in the right folder | "File not found" errors | Run `pwd` to check; `cd` to the right place |
| Forgot to activate virtual environment | "Module not found" errors | Run `source .venv/bin/activate` |
| Typo in a command | "Command not found" | Use Tab-completion; copy-paste from this handbook |
| Wrong file path | "No such file or directory" | Check with `ls`; use absolute paths if unsure |
| Missing Tm values | "Tm values are REQUIRED" error | Always provide `--min-tm`, `--opt-tm`, `--max-tm` |
| Pasting `#` comment lines | `zsh: command not found: #` | Only paste lines without `#` |
| Forgot `chmod +x` | "permission denied" | Run `chmod +x` on the script first |
| Stuck at `quote>` prompt | Pasted a line with an apostrophe | Press Ctrl-C to cancel, re-run the command |

---

### 0.8 Quick Reference: The 5 Commands You Will Use Most

| Command | What it does | When to use it |
|---------|----------------|----------------|
| `source .venv/bin/activate` | Enter the pipeline environment | Every time you open a new Terminal |
| `./install_macos.sh` | Install all dependencies | First time only |
| `./download_data.sh` | Download all data (~15 GB) | First time only |
| `python -m methyl_panel.pipeline --steps all ...` | Run the pipeline | To design primers |
| `python -m methyl_panel.pipeline --help` | Show all options | When unsure |

---

## Chapter 1 — Installation

### 1.1 Prerequisites

- **macOS** (Intel or Apple Silicon)
- **Xcode Command Line Tools:** Install with `xcode-select --install` (a dialog will appear — click Install)
- **~20 GB free disk space** (15 GB for data, 5 GB for tools and genome)

### 1.2 Step-by-Step Installation

```bash
git clone https://github.com/ScottAilliet/cell-specific-epigenetic-dpcr-assay-design-pipeline.git
cd cell-specific-epigenetic-dpcr-assay-design-pipeline
chmod +x install_macos.sh
./install_macos.sh
```

> **Do not paste the `#` comment lines** — they are explanations, not commands. Only paste the lines without `#`. If you paste a `#` line, you will see `zsh: command not found: #` which is harmless but confusing.

The script will:
1. Install Homebrew (if not present)
2. Install system tools: samtools, htslib (includes tabix), bedtools
3. Create a Python virtual environment (`.venv`)
4. Install Python packages: primer3-py, openpyxl, pandas, numpy, reportlab, pysam, scipy
5. Install wgbstools from GitHub source (compiles C++ tools + pip install)
6. Verify all installations

This takes about 10–15 minutes depending on your internet speed.

### 1.3 Verify the Installation

```bash
source .venv/bin/activate

python -c "from methyl_panel.pipeline import main; print('Pipeline OK')"

python -c "import primer3; print(f'primer3-py {primer3.__version__}')"

samtools --version | head -1
```

If all commands print a version or "Pipeline OK", you're ready.

Also verify wgbstools:

```bash
wgbstools --version
```

If `wgbstools: command not found`, the install script should have created a symlink. If not:

```bash
ln -sf "$(pwd)/wgbs_tools/wgbstools" .venv/bin/wgbstools
wgbstools --version
```

### 1.2 Updating the Pipeline

When a new version is released (bug fixes, new features), update your local copy:

```bash
cd ~/cell-specific-epigenetic-dpcr-assay-design-pipeline
git pull
```

If `git pull` fails with "Your local changes would be overwritten", see the troubleshooting section (Chapter 7). After updating, you may need to re-run `./install_macos.sh` if dependencies changed.

To verify you have the latest version:

```bash
git log --oneline -1
```

Check the [change log](#appendix-c--change-log) at the end of this handbook to see what changed in each version.

---

## Chapter 2 — Data Download

### 2.1 What Data You Need

The pipeline requires:

| Data | Size | Purpose |
|------|------|---------|
| GSE186458 beta files (207 samples) | ~2.4 GB | WGBS methylation data for DMR discovery |
| GSE186458 blocks file | ~1 MB | Genomic block definitions (205 CpG blocks) |
| hg19 genome FASTA | ~3 GB | Reference genome for sequence fetching |
| dbSNP common variants VCF | ~10 GB | SNP screening (Step 6) |
| UCSC segdup + RepeatMasker tracks | ~50 MB | Repeat/segdup annotation (Step 6b) |

**Total: ~10–15 GB**

### 2.2 Download Everything

```bash
cd ~/cell-specific-epigenetic-dpcr-assay-design-pipeline
source .venv/bin/activate

chmod +x download_data.sh
./download_data.sh
```

This downloads:
1. Blocks file from GEO FTP
2. All 207 beta files from GEO (using the manifest in `data/sample_manifest.csv`)
3. hg19 genome via wgbstools
4. dbSNP common variants (MAF ≥ 1%) from NCBI
5. UCSC segmental duplications track
6. UCSC RepeatMasker track

**This takes 30–60 minutes** depending on your internet speed. The script skips files that are already downloaded, so you can re-run it if the download is interrupted.

> **Note:** The download script uses `wget -q` (quiet mode) for beta files, so there are no progress bars. You will see `[$COUNT/$TOTAL] OK:` lines as each file downloads. The dbSNP download (~10 GB) shows a progress bar.

### 2.3 The DMR Input — Two Options

The pipeline needs DMR blocks (differentially methylated regions) as input. There are two ways to provide them:

**Option A — Discover DMRs from scratch (recommended, fully self-contained):**

Use `--discover-dmrs` to run Step 0, which uses `wgbstools find_markers` to discover DMRs directly from the 207 beta files you downloaded. No external file needed — this is the recommended workflow for v3.2.0.

**Option B — Load a pre-computed DMR Excel file:**

Use `--dmr-xlsx` to load a pre-computed Excel file (Step 1). This is useful if you already have DMR blocks from a previous analysis. The file must be in v7.9 format (one sheet per cell type) or Block_Summary format.

### 2.4 If You Already Have Some Data

If you already have beta files or the blocks file from a previous project, you can copy them into the right locations:

```bash
cp /path/to/your/beta_files/*.beta data/beta_files/
cp /path/to/GSE186458_blocks.s205.bed.gz data/
```

Then run `./download_data.sh` — it will skip what you already have and download the rest.

---

## Chapter 3 — Running the Pipeline

### 3.1 The Pipeline Steps

| Step | Module | What it does | Needs full genome? |
|------|--------|-------------|-------------------|
| 0 | `phase0_dmr_discovery` | Discover DMRs from raw beta files (wgbstools find_markers) | No |
| 1 | `phase1_dmr_loader` | Load DMR blocks from Excel (legacy alternative to Step 0) | No |
| 2 | `phase2_bisulfite_convert` | Bisulfite-convert genomic sequences to 6 strands | Yes (samtools faidx) |
| 3 | `phase3_primer3_design` | Design primers with Primer3 | No (uses step 2 output) |
| 4 | `bisearch_epcr_core` | Hash-based ePCR specificity screening (16-mer, no external aligner) | Yes (.npy genome arrays) |
| 5 | `phase5_structure` | Screen for hairpin secondary structures | No |
| 6 | `phase6_snp` | Screen for common SNPs (dbSNP, MAF ≥ 1%) | Yes (dbSNP VCF) |
| 6b (61) | `phase8_repeats` | Annotate repeats and segmental duplications (UCSC tracks) | Yes (segdup + rmsk BED) |
| 7 | `phase7_dimer` | Predict primer-dimer formation (DimerDetective) | No |
| 8 | `output_xlsx` | Generate U-assays-style Excel output (~79 columns) | No |
| 9 | `output_pdf` | Generate U-assays-style PDF (one page per primer pair) | No |

**Step 0 vs Step 1:** You use either Step 0 (`--discover-dmrs`, recommended) or Step 1 (legacy, without `--discover-dmrs`), not both. Step 0 discovers DMRs from raw beta files. Step 1 loads them from a pre-computed Excel file. Both produce the same `dmr_blocks.json` output that steps 2–9 consume.

**Steps that work without ePCR genome arrays/dbSNP:** 0, 1, 2, 3, 5, 7, 8, 9 (need only beta files + blocks file + hg19 genome)

**Steps that need the full data download:** 4 (ePCR .npy genome arrays), 6 (dbSNP VCF), 6b (UCSC BED tracks)

### 3.2 Running with DMR Discovery (Recommended)

This is the recommended workflow for v3.2.0 — fully self-contained, no external DMR file needed:

```bash
source .venv/bin/activate

python -m methyl_panel.pipeline --steps all --discover-dmrs \
    --cell-type MONO \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/MONO/
```

**Required arguments for `--discover-dmrs`:**
- `--cell-type`: One of `MONO`, `BCELL`, `NK`, `GRAN`, `CD3T`, `CD8T`, `CD4T` (required)
- `--genome`: Path to the hg19 genome FASTA (can be `.gz`; default: `data/hg19/hg19.fa.gz`)
- `--min-tm`, `--opt-tm`, `--max-tm`: Tm range in °C — **required, no defaults**

**DMR discovery arguments (optional):**
- `--beta-dir`: Directory with .beta files (default: `data/beta_files/`)
- `--blocks-file`: Path to blocks BED file (default: `data/GSE186458_blocks.s205.bed.gz`)
- `--groups-csv`: Path to full atlas groups CSV (default: `data/sample_groups.csv`)
- `--threads`: Number of threads for find_markers (default: 2)
- `--top-markers`: Number of top markers per cell type (default: 1000)
- `--delta-means`: Min methylation difference between target and background (default: 0.3)
- `--unmeth-mean-thresh`: Target mean methylation must be below this (default: 0.15)
- `--meth-mean-thresh`: Background mean methylation must be above this (default: 0.65)
- `--min-bg-subgroup-meth`: Reject blocks where any background blood cell type subgroup has mean methylation below this (default: 0.70). Set to 0 to disable. For hypo markers, every background subgroup must be ≥ this value. For hyper markers (`--only-hyper`), every background subgroup must be ≤ (1 - this value).
- `--max-bg-samples`: Max background samples for per-CpG extraction (default: 30)
- `--use-full-atlas`: Use the full 207-sample atlas as find_markers background instead of blood-only (36 samples). Per-CpG extraction and per-subgroup filter still use blood-only background. **Not recommended for hypo markers** — atlas delta ranking selects regions most different from non-blood tissues, not from blood cells. Blood-only ranking is correct for blood separation. Remains available for experimentation.
- `--only-hyper`: Discover hypermethylated markers (target methylated, background unmethylated) instead of the default hypomethylated markers (target unmethylated, background methylated). The per-subgroup filter and DMS score are automatically inverted.

**What Step 0 does:**
1. Generates a groups CSV file — target samples (your cell type) vs background (other blood cell types only, 36 blood samples). Overlapping cell populations (e.g. T-CD3 for CD4T/CD8T) are excluded from the background.
2. Generates a beta list file — paths to all .beta files
3. Runs `wgbstools find_markers` — discovers DMRs with quality filters: delta_means ≥ 0.3, target mean < 0.15 (hypo) or > 0.65 (hyper), background mean > 0.65 (hypo) or < 0.15 (hyper), min 3 CpGs, top 1000. Uses `--only_hypo` (default) or `--only_hyper` (with `--only-hyper`). Ranks by blood delta — the correct metric for blood cell separation.
4. Parses the output BED file
5. Excludes sex chromosomes (chrX, chrY, chrM) — these have copy number or methylation differences between sexes that make dPCR quantification unreliable
6. Extracts per-CpG methylation from beta files for each DMR (blood-only background)
7. Computes DMS scores — for hypo markers: target near-zero, background near-one; for hyper markers: target near-one, background near-zero. Both include consistency and coverage components.
8. Filters by DMS score (DMS >= 0.6) — rejects blocks with poor per-CpG data quality (low coverage, high variability)
9. Computes per-subgroup background methylation — for hypo markers, verifies that EACH background blood cell type is sufficiently methylated (≥ 0.70); for hyper markers, verifies each is sufficiently unmethylated (≤ 0.30). Rejects blocks where any single subgroup fails, and also rejects blocks with missing subgroup data (cannot verify = reject).
10. Saves `dmr_blocks.json` — same format as Step 1, with a `direction` field (`"U"` for hypo, `"M"` for hyper), so steps 2–9 work unchanged

**Why delta values are lower with blood-only background:** The background is restricted to other blood cell types (36 samples), not the full 207-sample atlas. Blood cells are more epigenetically similar to each other than to non-blood tissues. MONO vs T cells gives delta ~0.88, while MONO vs liver gives delta ~0.95. CD4T vs CD8T gives delta ~0.85 because CD4 and CD8 T cells are closely related lineages. **This is correct for a WBC assay** — the deltas that matter are against cells that will be in the same tube. The per-subgroup filter ensures that even the most similar background cell type is sufficiently methylated.

**Runtime:** ~5 minutes per cell type with 2 threads on a Mac (find_markers ~1 min, per-CpG extraction ~4 min).

### 3.3 Running with a Pre-Computed Excel File

If you already have a DMR Excel file (v7.9 format):

```bash
python -m methyl_panel.pipeline --steps all \
    --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/
```

### 3.4 Running Without ePCR Genome Arrays and dbSNP

If you haven't downloaded the ePCR genome arrays or dbSNP VCF yet, you can still run steps 0, 2, 3, 5, 7, 8, 9:

```bash
python -m methyl_panel.pipeline --steps 0,2,3,5,7,8,9 --discover-dmrs \
    --cell-type MONO \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/MONO/
```

The ePCR_verdict and common_variant_score columns will be NULL in the output. You can run steps 4 and 6 later once the data is downloaded.

### 3.5 Running Specific Steps

You can run any combination of steps. This is useful for re-running just the QC steps after changing parameters:

```bash
# DMR discovery only:
python -m methyl_panel.pipeline --steps 0 --discover-dmrs \
    --cell-type MONO --output-dir results/MONO/

# Primer design only (after DMR discovery):
python -m methyl_panel.pipeline --steps 2,3 \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/MONO/

# QC steps only:
python -m methyl_panel.pipeline --steps 5,7 \
    --output-dir results/MONO/

# Output only:
python -m methyl_panel.pipeline --steps 8,9 \
    --output-dir results/MONO/
```

Each step reads from and writes to JSON files in the output directory, so steps can be run independently as long as the prerequisite JSON files exist.

### 3.6 Running a Single Cell Type

```bash
python -m methyl_panel.pipeline --steps all --discover-dmrs \
    --cell-type CD8T \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/CD8T/
```

Cell type IDs and their target samples in the atlas. Background is always the other blood cell types (36 total blood samples minus target samples):

| Cell type | ID | Atlas groups (target) | # Target | # Background |
|-----------|-----|------------------------|----------|--------------|
| Monocytes | `MONO` | Blood-Monocytes | 3 | 33 |
| B cells | `BCELL` | Blood-B, Blood-B-Mem | 5 | 31 |
| NK cells | `NK` | Blood-NK | 3 | 33 |
| Granulocytes | `GRAN` | Blood-Granulocytes | 3 | 33 |
| CD3 T cells (pan-T) | `CD3T` | All Blood-T-* groups | 22 | 14 |
| CD8 T cells | `CD8T` | Blood-T-CD8, Blood-T-Eff-CD8, Blood-T-EffMem-CD8, Blood-T-Naive-CD8 | 10 | 26 |
| CD4 T cells | `CD4T` | Blood-T-CD4, Blood-T-CenMem-CD4, Blood-T-EffMem-CD4, Blood-T-Naive-CD4 | 10 | 26 |

When using `--discover-dmrs`, the pipeline automatically merges these atlas groups into one target population. The background is restricted to other blood cell types only (36 blood samples total) — not the full atlas — because this is a WBC assay that needs to distinguish the target from other blood cells.

### 3.7 Running All 7 Cell Types

To run all 7 immune cell types with DMR discovery:

```bash
# Hypomethylated markers (blood-only background, strict 0.70 filter)
for CT in MONO BCELL NK GRAN CD3T CD8T CD4T; do
    python -m methyl_panel.pipeline --steps all --discover-dmrs \
        --cell-type $CT \
        --genome data/hg19/hg19.fa.gz \
        --min-tm 58 --opt-tm 60 --max-tm 62 \
        --output-dir results/$CT/
done
```

```bash
# Hypermethylated markers (target methylated, background unmethylated)
for CT in MONO BCELL NK GRAN CD3T CD8T CD4T; do
    python -m methyl_panel.pipeline --steps all --discover-dmrs \
        --cell-type $CT \
        --genome data/hg19/hg19.fa.gz \
        --min-tm 58 --opt-tm 60 --max-tm 62 \
        --only-hyper \
        --output-dir results_hyper/$CT/
done
```

Each cell type's results are saved in its own subdirectory (`results/MONO/`, `results/BCELL/`, etc. for hypo; `results_hyper/MONO/`, etc. for hyper). The find_markers step runs once per cell type (~1 min each), and per-CpG extraction takes ~4 min per cell type. Total runtime: ~35 minutes for all 7 (hypo), plus another ~35 minutes for hyper if you run both.

### 3.7 Choosing Tm Values

The Tm range is **required** — the pipeline refuses to run without it. This forces an explicit decision about the Tm window for each run.

| Platform | Recommended Tm range |
|----------|---------------------|
| Primer3Plus default | `--min-tm 55 --opt-tm 60 --max-tm 65` |
| Roche Digital LightCycler | `--min-tm 58 --opt-tm 60 --max-tm 62` |
| QIAcuity dPCR | `--min-tm 56 --opt-tm 58 --max-tm 60` |

To use Roche DLC salt conditions, load from a Primer3Plus settings file that has the Roche preset, or modify the salt parameters in `config.py`.

### 3.8 Using a Primer3Plus Settings File

If you have a Primer3Plus settings file (exported from the Primer3Plus web interface):

```bash
python -m methyl_panel.pipeline --steps all \
    --dmr-xlsx data/WBC_Panel_Top200_v7.9.xlsx \
    --genome data/hg19/hg19.fa.gz \
    --settings Primer3Plus_Settings.txt \
    --output-dir results/
```

The settings file overrides all Primer3 parameters, including Tm. If the settings file has Tm values, you don't need `--min-tm` etc.

### 3.9 Testing with a Small Subset

To test the pipeline quickly without processing all blocks:

```bash
python -m methyl_panel.pipeline --steps all --discover-dmrs \
    --cell-type MONO \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --max-blocks 50 \
    --output-dir results/test/
```

`--max-blocks 50` processes only the first 50 blocks that pass the CpG count filter. Note: the first blocks by delta_means may have only 3–4 CpGs (below the `--min-cpg 5` threshold), so use at least `--max-blocks 50` to ensure some blocks pass the filter and produce primers.

To skip re-running find_markers on subsequent runs (reuse the existing BED output):

```bash
python -m methyl_panel.pipeline --steps 0,2,3,5,7,8,9 --discover-dmrs \
    --cell-type MONO \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --max-blocks 50 \
    --skip-find-markers \
    --output-dir results/test/
```

### 3.10 Exporting DMR Regions to Excel

The pipeline's main output is primer assays (Section 4). If you also want the DMR regions themselves — the top N blocks per cell type with per-CpG methylation values — use the `export_dmr_excel.py` script.

**Export all 7 cell types (one sheet each + per-CpG detail):**

```bash
source .venv/bin/activate

python scripts/export_dmr_excel.py results/*/dmr_blocks.json --top 200 -o dmr_regions.xlsx
```

**Export a single cell type:**

```bash
python scripts/export_dmr_excel.py results/MONO/dmr_blocks.json --top 200 -o dmr_MONO.xlsx
```

**Arguments:**
- First argument: path(s) to `dmr_blocks.json` file(s). Supports glob patterns (`results/*/dmr_blocks.json`).
- `--top N`: Number of top DMR regions per cell type (default: 200). Regions are ranked by `delta_means` (methylation difference between target and background).
- `-o PATH`: Output Excel file path (default: `dmr_regions.xlsx`).

**Output Excel structure:**

| Sheet | Content |
|-------|---------|
| `All blocks summary` | All DMR blocks across all cell types (one row per block) |
| `MONO`, `BCELL`, `NK`, ... | One sheet per cell type, top N blocks (one row per block) |
| `Per-CpG detail` | Every CpG site across all blocks (one row per CpG) |

**Block-level columns** (in each cell type sheet and the summary):

| Column | Description |
|--------|-------------|
| Cell type | MONO, BCELL, NK, etc. |
| Rank | DMR rank by delta_means (1 = largest methylation difference) |
| Seq ID | Block identifier (e.g. MONO_0001) |
| Chromosome | Chromosome (e.g. chr16) |
| Start | Genomic start position (1-based, hg19) |
| End | Genomic end position (1-based, hg19) |
| Block length (bp) | Length of the DMR block |
| # CpGs | Number of CpG sites in the block |
| Gene | Associated gene name |
| Annotation | Genomic annotation (e.g. promoter-TSS, intron, intergenic) |
| Target mean methylation | Mean methylation across target samples (should be near 0 for hypo-DMRs) |
| Background mean methylation | Mean methylation across all background blood samples (should be near 1) |
| Delta means | Difference between background and target methylation |
| DMS score | 5-component score (0–1, higher = better primer candidate) |
| BG: Blood-B | Mean methylation of B cell background subgroup |
| BG: Blood-NK | Mean methylation of NK cell background subgroup |
| BG: Blood-T-CD8 | Mean methylation of CD8 T cell background subgroup |
| ... | One column per background blood cell type subgroup |
| Min subgroup methylation | Lowest methylation across all background subgroups (should be > 0.50) |
| Worst background subgroup | The blood cell type with the lowest methylation (most likely to cause false positives) |

**Per-CpG columns** (in the Per-CpG detail sheet):

| Column | Description |
|--------|-------------|
| Cell type | MONO, BCELL, NK, etc. |
| Seq ID | Parent block identifier |
| Rank | Parent block rank |
| Chromosome | Chromosome |
| Block start / Block end | Parent block coordinates |
| Gene | Associated gene |
| CpG label | CpG_1, CpG_2, etc. (position within block) |
| CpG position | Genomic position (hg19) |
| CpG index | Global CpG index in the genome |
| Target mean beta | Mean methylation at this CpG across target samples |
| Background mean beta | Mean methylation at this CpG across background samples |
| Delta beta | Absolute methylation difference at this CpG |

**Prerequisites:** You must have already run Step 0 (`--discover-dmrs`) for each cell type you want to export. The script reads the `dmr_blocks.json` files from the output directories.

---

## Chapter 4 — Understanding the Output

### 4.1 Output Files

After running the pipeline, the output directory contains:

| File | Description |
|------|-------------|
| `groups_MONO.csv` | Groups file for find_markers (target vs background) — Step 0 |
| `beta_list.txt` | List of beta file paths — Step 0 |
| `find_markers_output/` | wgbstools find_markers output (BED files, params) — Step 0 |
| `dmr_blocks.json` | DMR blocks (from Step 0 or Step 1) |
| `converted_sequences.json` | Bisulfite-converted sequences from Step 2 |
| `primers.json` | All primer pairs with QC results (updated by each step) |
| `primer_assays.xlsx` | U-assays-style Excel output (~79 columns) |
| `primer_assays.pdf` | U-assays-style PDF (one page per primer pair) |

The `export_dmr_excel.py` script (Section 3.10) can additionally produce `dmr_regions.xlsx` — a separate Excel with the DMR regions and per-CpG methylation values, which is useful if you want to inspect the DMRs themselves rather than the primer assays.

### 4.2 Excel Output Columns

The Excel output (`primer_assays.xlsx`) contains ~79 columns organized in four groups. Key columns are documented below.

**Identity (2 columns):**

| Column | Description |
|--------|-------------|
| `assay_id` | Unique assay identifier (e.g. MONO_0001_001) |
| `seq_id` | Source DMR block ID |

**DMR context (19 columns):** `cell_type`, `chrom`, `dmr_start`, `dmr_end`, `dmr_length_bp`, `gene`, `direction`, `dmr_rank`, `delta_means`, `delta_quants`, `delta_maxmin`, `ttest_pval`, `mwtest_pval`, `mvalue_ttest`, `target_mean_meth`, `background_mean_meth`, `dms_score`, `dmr_num_cpgs`, `annotation`

**Amplicon (8 columns):** `amplicon_start`, `amplicon_end`, `cpgs_in_amplicon`, `amp_target_mean_beta`, `amp_background_mean_beta`, `amp_delta`, `min_non_target_delta`, `worst_celltype`

**Primers + screening (50 columns):**

| Column | Description | Step |
|--------|-------------|------|
| `template_used` | Bisulfite strand (SM, AM, SU, AU) | 3 |
| `assay` | M (methylated) or U (unmethylated) | 3 |
| `left_primer` / `right_primer` | Primer sequences | 3 |
| `left_primer_cpg_display` / `right_primer_cpg_display` | CpG notation (y = CpG C) | 3 |
| `left_tm_C` / `right_tm_C` | Primer Tm (°C) | 3 |
| `product_size_bp` | Amplicon size (bp) | 3 |
| `c_total_tail`, `c_total`, `left_c_total`, `right_c_total`, `left_c_tail`, `right_c_tail` | CpG counts | 3 |
| `sense_meth_mismatch_score` etc. (4 columns) | Mismatch scores vs 4 bisulfite strands | 3 |
| `ePCR_verdict` | PASS or FAIL (template-aware) | 4 |
| `ePCR_product_count` | Number of PCR products found | 4 |
| `ePCR_products` | Product coordinates and details | 4 |
| `ePCR_error` | Error message (empty on PASS) | 4 |
| `ePCR_elapsed_sec` | Screening time (seconds) | 4 |
| `left_structure_mfe_kcal_mol` / `right_structure_mfe_kcal_mol` | Hairpin MFE | 5 |
| `left_self_dimer_dg`, `right_self_dimer_dg`, `heterodimer_dg`, `amplicon_mfe` | Extended structure metrics | 5 |
| `primer_dimer_prediction` | Dimer risk tier: high, medium, low | 7 |
| `primer_dimer_end_min_dg` | Minimum end stability ΔG (kcal/mol) | 7 |
| `heterodimer_full_dg`, `left_homodimer_dg`, `right_homodimer_dg` | Extended dimer metrics | 7 |
| `common_variant_score` | SNP score: 00=clean, 10=body, 01=3'end, 11=both | 6 |
| `snp_count` | Number of dbSNP variants overlapping primer | 6 |
| `segdup_overlap` | TRUE if primer overlaps a segmental duplication | 6b |
| `repeat_overlap` | TRUE if primer overlaps a repeat element | 6b |
| `repeat_family` | Repeat family name (e.g. SINE/Alu, LINE/L1) | 6b |
| `amplicon_snp_count`, `cpg_snp_count` | Amplicon-level SNP counts | 6/6b |
| `amplicon_segdup_overlap`, `amplicon_repeat_overlap`, `amplicon_repeat_family` | Amplicon-level repeat/segdup annotations | 6b |
| `mismatch_profile`, `mismatch_detail` | Detailed mismatch information | 4 |

### 4.3 PDF Output

Each page of the PDF contains:
- **Assay header:** assay_id, seq_id, template, cell type
- **6-strand sequence visualization:** S, SM, SU, A, AM, AU with CpG sites highlighted
- **Primer information table:** sequences, Tm, length, GC%, CpG counts
- **Assay characteristics:** product size, CpG counts, mismatch scores
- **QC filter results:** ePCR, structure MFE, dimer prediction, SNP score
- **DimerDetective footer:** thermodynamic parameters and risk tier explanation

### 4.4 Interpreting QC Results

**ePCR specificity (Step 4):**
- `ePCR_verdict = PASS`: Single PCR product found on the intended methylation state (methylated for SM/AM templates, unmethylated for SU/AU templates)
- `ePCR_verdict = FAIL`: No product, product on wrong methylation state, or multiple products — primer may lack specificity
- `ePCR_product_count`: Number of PCR products (1 = ideal, 0 = no amplification, >1 = off-target)
- `ePCR_products`: Formatted product coordinates (chrom:start-end, size, methylation state, mismatches)
- `NULL`: ePCR genome arrays not available — screening was skipped

**Secondary structure (Step 5):**
- MFE ≥ −1.5 kcal/mol: PASS (no stable hairpin)
- MFE < −1.5 kcal/mol: FAIL (stable hairpin may form, reducing PCR efficiency)
- More negative = more stable = worse

**Primer-dimer (Step 7, DimerDetective):**
- `low`: end_min_dg > −0.18 kcal/mol — no dimer risk
- `medium`: −2.48 < end_min_dg ≤ −0.18 — mixed zone, may form dimers
- `high`: end_min_dg ≤ −2.48 kcal/mol — high dimer risk
- Passes filter if end_min_dg ≥ −1.0 kcal/mol (conservative cutoff)

> **Bisulfite caveat:** DimerDetective was validated on 40–60% GC primers. Bisulfite MSP primers typically have 10–25% GC, outside the validation range. Treat dimer predictions as indicative, not definitive.

**Common SNPs (Step 6):**
- Score `00`: No common SNPs in primer — clean
- Score `10`: SNP in primer body (not at 3' end) — moderate risk
- Score `01`: SNP in 3' end (last 5 nt) — high risk, may disrupt primer binding
- Score `11`: SNP in both body and 3' end — highest risk
- `NULL`: dbSNP VCF not available — screening was skipped

**Repeats and segmental duplications (Step 6b):**
- `segdup_overlap = TRUE`: Primer overlaps a segmental duplication — may map to multiple genomic locations
- `repeat_overlap = TRUE`: Primer overlaps a repeat element (e.g. SINE/Alu, LINE/L1)
- `repeat_family`: Name of the overlapping repeat family (NULL if no overlap)
- `snp_count`: Total dbSNP variants overlapping the primer region
- `amplicon_snp_count`, `cpg_snp_count`: SNP counts at the amplicon and CpG-site level
- `NULL`: UCSC track files not available — screening was skipped

### 4.5 M-assays vs U-assays

The pipeline designs primers on 4 bisulfite-converted templates:
- **Methylated (M assays):** from SM (sense methylated) or AM (antisense methylated) — CpG C's are preserved
- **Unmethylated (U assays):** from SU (sense unmethylated) or AU (antisense unmethylated) — all C's converted to T

For hypomethylated DMRs (the target cell type has low methylation), most primers will be M-assays. This is expected: the methylated templates retain CpG sites that provide the sequence variation needed for Primer3 to design primers within the Tm window. The unmethylated templates become very AT-rich (all C→T), making it harder for Primer3 to find suitable primers.

---

## Chapter 5 — The 6 Bisulfite Strands

The pipeline converts each genomic region into 6 strands for primer design:

| Strand | Name | Description |
|--------|------|-------------|
| S | Sense | Original top strand (unconverted) |
| SM | Sense Methylated | Top strand, C→T except CpG C's preserved |
| SU | Sense Unmethylated | Top strand, all C→T |
| A | Antisense | Original bottom strand (reverse complement) |
| AM | Antisense Methylated | Bottom strand, C→T except CpG C's preserved |
| AU | Antisense Unmethylated | Bottom strand, all C→T |

**MSP primer design:**
- **Methylated primers** (M assays): designed from SM or AM — CpG C's are preserved, so the primer only matches methylated DNA
- **Unmethylated primers** (U assays): designed from SU or AU — all C's converted to T, so the primer only matches unmethylated DNA

The `template_used` column in the output tells you which strand each primer was designed from.

---

## Chapter 6 — Configuration

### 6.1 Primer3PlusConfig

All 161 Primer3Plus parameters are configurable. Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `primer_min_tm` | None (required) | Minimum Tm (°C) |
| `primer_opt_tm` | None (required) | Optimal Tm (°C) |
| `primer_max_tm` | None (required) | Maximum Tm (°C) |
| `primer_min_size` | 18 | Minimum primer length (bp) |
| `primer_opt_size` | 20 | Optimal primer length (bp) |
| `primer_max_size` | 23 | Maximum primer length (bp) |
| `primer_product_size_range` | [(60, 150)] | Amplicon size range (bp) |
| `primer_min_gc` | 30.0 | Minimum GC% |
| `primer_max_gc` | 80.0 | Maximum GC% |
| `primer_max_poly_x` | 5 | Max consecutive identical bases |
| `primer_gc_clamp` | 1 | Require GC clamp at 3' end |
| `primer_max_end_stability` | 9.0 | Max end stability (kcal/mol) |
| `primer_max_self_any` | 8.0 | Max self-complementarity |
| `primer_max_self_end` | 3.0 | Max self-end complementarity |
| `primer_pair_max_diff_tm` | 1.0 | Max Tm difference between primers |

### 6.2 Salt Presets

| Preset | Na+ (mM) | Mg2+ (mM) | dNTP (mM) | Primer (nM) | Correction |
|--------|---------|----------|----------|------------|------------|
| `primer3plus` | 50 | 3 | 1.2 | 250 | SantaLucia |
| `roche_dlc` | 100 | 4.5 | 1.2 | 250 | Owczarzy 2008 |

Default: `primer3plus`. To use Roche DLC conditions, set `salt_preset = "roche_dlc"` in config or load from a Primer3Plus settings file with Roche parameters.

### 6.3 PipelineConfig

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_cpg` | 5 | Minimum CpGs per DMR block |
| `preferred_min_cpg` | 7 | Preferred minimum CpGs |
| `min_cpg_per_primer` | 2 | Minimum CpGs per primer |
| `min_cpg_pair_total` | 4 | Minimum total CpGs per primer pair |
| `max_primers_per_block` | 10 | Max primer pairs per block |
| `delta_means` | 0.3 | DMR delta_means threshold |
| `only_hypo` | True | Only hypomethylated markers |

---


## Chapter 6b — Standalone ePCR Screening and Annotation

> **Who this is for:** You have already run the main pipeline (Chapter 3) and produced the master primer atlas. Now you need to screen all primers for bisulfite-genome specificity, compute DMS sub-scores, annotate with UCSC reference tracks, and extract per-CpG beta values. This chapter documents the standalone tools that perform these post-pipeline steps.

---

### 6b.1 Overview

The standalone ePCR screening workflow is a post-pipeline step that takes the master primer atlas (output of the batch design script) and:

1. **Screens** all primers for specificity against the bisulfite-converted genome using `bisearch_epcr.py`
2. **Selects** passing primers (keeps ALL that pass, no deduplication or caps)
3. **Computes DMS sub-scores** (A-E) for each DMR block using `compute_dms_subscores.py`
4. **Annotates** passing primers with UCSC hg19 reference tracks using `annotate_primers.py`
5. **Extracts per-CpG beta values** for all cell types using `extract_cpg_beta.py`
6. **Merges** everything into a final atlas using `merge_screening_results.py`

The consolidated runner `run_full_pipeline.sh` automates all six steps. It does NOT re-run the main pipeline or re-screen primers if screening results already exist.

**Key distinction from the pipeline's internal Step 4:** The pipeline's internal Step 4 (`methyl_panel/bisearch_epcr_core.py`) uses the same hash-based numpy indexing as the standalone tool but is called as a pipeline step with a 2-way verdict (PASS / FAIL). The standalone `bisearch_epcr.py` tool offers additional CLI options for batch screening with early stopping, cell-type grouping, and selection. See Section 6b.3 for details.

---

### 6b.2 The bisearch_epcr.py Tool

**File:** `bisearch_epcr.py` (repository root)

**What it does:** Screens bisulfite PCR primer pairs for specificity against the bisulfite-converted genome. It recreates the BiSearch ePCR algorithm (Tusnady et al. 2005 NAR 33:e9; Aranyi et al. 2006 BMC Bioinformatics 7:431) using 16-mer hash-based indexing with sorted numpy arrays and vectorised verification.

**No external aligner dependency:** This tool does not require bowtie2, samtools, or any external aligner. It builds its own in-memory index of the bisulfite-converted genome using numpy arrays and performs seed lookup via `np.searchsorted`. This makes it fully portable on macOS without any external bioinformatics tools.

**Algorithm per primer pair:**

1. **Seed extraction:** Take the 3'-most 16 bases of the forward primer (FP) and of the reverse-complemented reverse primer (RC(RP)). These are the seeds.
2. **Seed variant generation:** Generate all seed variants within the mismatch tolerance (max 3 mismatches, governed by the mismatch string `"1233333333333333"`). The mismatch string specifies the cumulative mismatch budget at each position: position 0 allows 1 mismatch, position 1 allows 2 cumulative, positions 2-15 allow 3 cumulative.
3. **Genome conversion:** For each strand (Watson forward, Crick reverse-complement) and each methylation state:
   - **CT (unmethylated):** All cytosines converted to thymine. Alphabet = {A, G, T}.
   - **CT_M (methylated):** CpG cytosines preserved, non-CpG cytosines converted to thymine. Alphabet = {A, C, G, T}.
4. **Seed lookup:** For each converted genome, look up all seed variants in the sorted index (`np.searchsorted`). Each hit is a potential primer binding site.
5. **Full-primer verification:** For each seed hit, verify the full primer (not just the seed) against the genomic sequence. Count mismatches across the entire primer length. Reject if mismatches exceed the tolerance.
6. **Product detection:** Find all (FP hit, RP hit) pairs on opposite strands within `MAX_PCR_LEN` (1000 bp) of each other. Each such pair is a predicted PCR product.
7. **Deduplication:** Products are deduplicated per (chrom, strand, fp_pos, rp_pos). A locus found in both methylation states (CT and CT_M) counts as one product, but the `states` field records all genome states where the product was found.
8. **Verdict assignment:** Based on the number of unique products and the template-aware logic (see Section 6b.3).

**Constants:**

| Constant | Value | Description |
|----------|-------|-------------|
| `SEED_LEN` | 16 | Length of the 3' seed region |
| `MAX_MM` | 3 | Maximum mismatches per primer |
| `MAX_PCR_LEN` | 1000 | Maximum PCR product size in bp |
| `MM_STRING` | `"1233333333333333"` | Cumulative mismatch budget per position |
| `STRAND_STATES` | `[('unmethylated', 'CT'), ('methylated', 'CT_M')]` | Genome states to search |

**Template codes:** S = Sense (Watson strand), A = Antisense (Crick strand), M = Methylated (CT_M), U = Unmethylated (CT). So SM = Sense-Methylated, SU = Sense-Unmethylated, AM = Antisense-Methylated, AU = Antisense-Unmethylated.

**CLI usage:**

```bash
# Screen a single primer pair:
python bisearch_epcr.py screen --fp ATCG... --rp GGTG... --genome-dir hg19_bs

# Screen a batch from CSV:
python bisearch_epcr.py screen --input primers.csv --output results.csv \
    --genome-dir hg19_bs --screen-per-ct 200

# Select passing primers from screening results:
python bisearch_epcr.py select --input results.csv --output selected.csv
```

**Batch screening with DMS-ranked early stopping:** When screening from a CSV, primers are processed in DMS-ranked order (highest delta_means first) within each cell type. The `--screen-per-ct` parameter (default 200) controls how many primers to screen per cell type. The `--max-pass-per-ct` parameter (default 0 = unlimited) controls how many passing primers to select per cell type before stopping early.

**Genome directory:** The `--genome-dir` argument points to a directory containing numpy `.npy` files, one per chromosome (e.g., `chr1.npy`, `chr2.npy`, ...). These are uint8 arrays of the original hg19 genome (conversion is done at runtime). Use `prepare_genomes.py` to generate these from a FASTA file.

---

### 6b.3 The 3-Way Template-Aware Verdict

The standalone ePCR tool uses a 3-way verdict system that is template-aware: it checks not only how many products a primer pair generates, but also whether the product is on the correct genome state (methylated or unmethylated) for the assay's intended template.

**Template-to-state mapping:** The `template_used` column in the input CSV encodes which genome state the assay was designed for:

| Template code | Intended genome state | Description |
|---------------|----------------------|-------------|
| SM, AM | methylated (CT_M) | Sense/Antisense, Methylated |
| SU, AU | unmethylated (CT) | Sense/Antisense, Unmethylated |

**Verdict logic:**

| Verdict | Condition | Meaning |
|---------|-----------|---------|
| **PASS** | Exactly 1 product AND (template unknown OR product on intended genome state only) | Primers amplify a unique locus on the correct methylation state |
| **FAIL** (wrong_genome) | Exactly 1 product, but on the wrong genome state | Primers amplify a unique locus, but on the wrong methylation state. For example, a UM assay (template SU) whose product is found only on the methylated (CT_M) genome. |
| **FAIL** (not_methylation_specific) | Exactly 1 product locus, but found on BOTH genome states | Primers amplify a unique locus, but cannot distinguish methylated from unmethylated DNA. The assay is not methylation-specific. |
| **FAIL** (generic) | 0 products or >1 products | No amplification (0 products) or off-target amplification (>1 products) |

**Why the 3-way verdict matters:** The pipeline's internal Phase 4 uses a 4-way verdict (SPECIFIC / CHECK / NON-SPECIFIC / NO_PRODUCT) that only counts products, without checking which genome state they are on. The standalone tool goes further: it verifies that the product is on the intended methylation state. This catches two failure modes that the 4-way verdict misses:

1. **wrong_genome:** A primer designed for the unmethylated template that only amplifies on the methylated genome. In the wet lab, this assay would only detect methylated DNA, not unmethylated, defeating the purpose of the assay.

2. **not_methylation_specific:** A primer that amplifies on both genome states. In the wet lab, this assay cannot distinguish methylated from unmethylated DNA, making it useless for methylation quantification.

**The `states` field:** Each `PCRProduct` object has a `states` set that records all genome states where the product was found. A product at the same locus on both CT and CT_M has `states = {'unmethylated', 'methylated'}`. The verdict logic checks whether the intended state (from `template_used`) is in this set, and whether the set contains only the intended state.

**When template is unknown:** If the `template_used` column is absent or empty, the tool falls back to a count-only check: PASS if exactly 1 product, FAIL otherwise. This preserves backward compatibility with primer sets that do not have template information.

---

### 6b.4 Selection Logic

The `select` subcommand selects passing primers from the screening results.

**Default behaviour: keep ALL passing primers.** There is no DMR deduplication (previously, only one primer per DMR was kept) and no cell-type cap (previously, a maximum of 20 per cell type was enforced). Both restrictions have been removed.

**Ranking within each cell type:**
1. `delta_means` descending (DMR importance: larger methylation difference = better marker)
2. `cpg_total` descending (primer quality: more CpG sites in the primer pair = better discrimination)

**What counts as "passing":** A primer passes if:
- `verdict == 'PASS'`
- No error
- Exactly 1 product
- If template is known: product is on the intended genome state only (not wrong_genome, not not_methylation_specific)

---

### 6b.5 The run_full_pipeline.sh Script

**File:** `run_full_pipeline.sh` (repository root)

This is the consolidated pipeline runner that automates all post-pipeline steps in a single command. It does NOT re-run the main pipeline or re-screen primers if screening results already exist.

**Prerequisites:** Before running, edit the USER CONFIG section at the top of the script:

```bash
EPCR_DIR=~/epcr_screening           # Directory with results_v2.csv, atlas xlsx, hg19_bs/
REPO_DIR=~/cell-specific-epigenetic-dpcr-assay-design-pipeline  # Cloned repository
BETA_DIR=~/path/to/beta/files        # Directory containing wgbstools .beta files
GROUPS_CSV=~/path/to/groups.csv      # Sample-to-cell-type mapping CSV
CPG_INDEX_FILE=~/path/to/CpG.bed.gz  # wgbstools CpG index file
```

**The 7 steps:**

| Step | Action | Script | Time |
|------|--------|--------|------|
| 1 | Pull latest code from GitHub | `git pull` | Seconds |
| 2 | Copy scripts to working directory | `cp` commands | Seconds |
| 3 | Re-run selection on existing screening results | `bisearch_epcr.py select` | Seconds |
| 4 | Compute DMS sub-scores | `compute_dms_subscores.py` | 30-60 min |
| 5 | Annotate primers with UCSC tracks | `annotate_primers.py` | 10-20 min |
| 6 | Extract per-CpG beta values | `extract_cpg_beta.py` | 10-20 min |
| 7 | Merge everything into final atlas | `merge_screening_results.py` | Seconds |

**Step 3 does NOT re-screen.** It only re-runs the `select` subcommand on the existing `results_v2.csv` with the new "keep all passing primers" logic. This is instant. The ePCR screening itself was already done and takes ~2.3 hours; it is not repeated.

**Usage:**

```bash
# Edit the USER CONFIG section, then:
bash run_full_pipeline.sh
```

---

### 6b.6 compute_dms_subscores.py

**File:** `scripts/compute_dms_subscores.py`

**What it does:** Computes the DMS A-E sub-scores for every unique DMR block in the atlas. This is the same calculation as the per-CpG DMS scoring in Phase 0, but run as a standalone post-pipeline step so that the sub-scores can be merged into the final atlas.

**Inputs:**
- Master primer atlas xlsx (with DMR block coordinates)
- Beta files directory (wgbstools `.beta` binary files, one per sample)
- Groups CSV (maps sample names to cell types; auto-detects columns containing "sample"/"file"/"gsm" and "group"/"type"/"cell")
- CpG index file (`CpG.bed.gz` from wgbstools, bgzipped + tabixed: columns chr, loc, site)

**Process:**
1. Reads the atlas xlsx to extract unique DMR blocks (cell_type, rank, chrom, dmr_start, dmr_end)
2. Loads the CpG index file and converts genomic coordinates to CpG indices via `genomic_to_cpg_range()`
3. For each DMR block, reads the beta files for all target and background samples
4. Computes A-E sub-scores using the `cleanliness_score()` function (identical to Phase 0)
5. Outputs `dms_subscores.csv`

**Beta file format:** Binary, 2 x uint8 per CpG (#meth, #covered), 1-based CpG index. Byte offset = (start_cpg - 1) * 2. CpGs with coverage < `MIN_COV` (5) are treated as missing data (NaN).

**Output columns:** cell_type, rank, chrom, dmr_start, dmr_end, start_cpg, end_cpg, n_cpgs, dms_score, score_A_tg_near0, score_B_bg_near1, score_C_tg_consist, score_D_bg_consist, score_E_coverage

**Runtime:** 30-60 minutes for ~47K DMR blocks (reads beta files for each block).

---

### 6b.7 annotate_primers.py

**File:** `scripts/annotate_primers.py`

**What it does:** Annotates passing primers with SNP, repeat, and segmental duplication overlap data from UCSC hg19 reference tracks. This fills 5 annotation columns that were previously empty in the atlas.

**UCSC tracks downloaded (one-time, cached locally):**

| Track | UCSC table | What it provides |
|-------|-----------|------------------|
| Common SNPs | snp151Common | SNPs with MAF >= 1% |
| All SNPs | snp151 | All known SNPs |
| Repeats | rmsk | RepeatMasker repetitive elements |
| Segmental duplications | genomicSuperDups | Segmental duplication regions |

**Primer binding position finding:** The atlas xlsx does not store explicit primer genomic positions. The script finds them by searching the bisulfite-converted genome within the DMR block + 500bp flanking:
- Forward primer (FP): searched as-is against CT/CT_M converted forward strand
- Reverse primer (RP): reverse-complemented first, then searched against CT/CT_M converted forward strand
- Also tries GA/GA_M conversions for antisense primers
- The hg19_bs/ directory contains numpy `.npy` files of the original hg19 genome (conversion done at runtime)

**Annotation columns filled:**

| Column | Type | Description |
|--------|------|-------------|
| `common_variant_score` | int | Count of common SNPs (MAF >= 1%) overlapping primer binding sites |
| `snp_count` | int | Total known SNPs overlapping primer binding sites |
| `segdup_overlap` | bool | True if primer overlaps a segmental duplication |
| `repeat_overlap` | bool | True if primer overlaps a repetitive element |
| `repeat_family` | str | Name of the repetitive element family (e.g., "SINE", "LINE", "LTR") if any |

**Only processes passing primers** (verdict = PASS from the screening results CSV).

**Output:** `primer_annotations.csv` with assay_id, seq_id, template_used, chrom, fp_start, fp_end, rp_start, rp_end, and the 5 annotation columns above.

**Runtime:** 10-20 minutes (dominated by UCSC track downloads on first run, ~100MB total).

---

### 6b.8 extract_cpg_beta.py

**File:** `scripts/extract_cpg_beta.py`

**What it does:** For each DMR block that has at least one passing assay, extracts per-CpG methylation beta values for every cell type and every sample in the atlas. This produces the raw data needed to verify DMS scores and to visualise methylation patterns across cell types.

**Inputs:**
- Screening results CSV (to identify selected DMR blocks: verdict = PASS)
- Beta files directory
- Groups CSV
- CpG index file

**Process:**
1. Reads the results CSV to identify unique DMR blocks with passing assays
2. Loads the CpG index file to map genomic coordinates to CpG indices
3. For each DMR block, reads beta files for ALL cell types, ALL samples
4. `read_beta_region()` returns (fracs, meth, cov) arrays: fracs is NaN where coverage < MIN_COV (5)
5. Outputs a long-format CSV

**Output columns:** dmr_id, cell_type_target, chrom, dmr_start, dmr_end, cpg_index, cpg_pos, cell_type, sample_name, beta_value, n_meth, n_cov, is_target

- `is_target`: True if this cell type is the target for this DMR block, False if it is background
- `beta_value`: methylation fraction (n_meth / n_cov), NaN if coverage < 5
- One row per (DMR block, CpG, sample) combination

**Runtime:** 10-20 minutes for selected DMR blocks.

---

### 6b.9 merge_screening_results.py

**File:** `scripts/merge_screening_results.py`

**What it does:** Merges ePCR screening results, DMS sub-scores, and primer annotations into the master primer atlas Excel file, producing the final screened atlas.

**Inputs:**
- Master primer atlas xlsx (original atlas with sheets: All_Primers, Best_Per_DMR, Summary, DMR_Blocks)
- `results_v2.csv` (full screening output)
- `results_v2_selected.csv` (selected panel, auto-generated by select subcommand)
- `dms_subscores.csv` (optional, from compute_dms_subscores.py)
- `primer_annotations.csv` (optional, from annotate_primers.py)

**Merge key:** `assay_id` + `left_primer` + `right_primer` + `template_used` (each row in the atlas is a unique primer-template combination)

**Columns removed from all sheets:**

| Column | Reason |
|--------|--------|
| `direction` | Unreliable; removed per user decision |
| `cleanliness_score` | Renamed to `dms_score`; duplicate removed |
| `sense_meth_mismatch_score` | Superseded by ePCR screening |
| `sense_unmeth_mismatch_score` | Superseded by ePCR screening |
| `anti_meth_mismatch_score` | Superseded by ePCR screening |
| `anti_unmeth_mismatch_score` | Superseded by ePCR screening |
| `bowtie_passes_filter` | Bowtie-specific; replaced by ePCR_verdict |
| `bowtie_intended_genome` | Bowtie-specific; replaced by ePCR_verdict |
| `mismatch_profile` | Bowtie-specific; replaced by ePCR screening |
| `mismatch_detail` | Bowtie-specific; replaced by ePCR screening |

**Columns retained:** The `assay` column (contains "M" or "U") is kept per user decision. It is distinct from the removed `direction` column.

**Columns added:**

| Column | Source | Description |
|--------|--------|-------------|
| `ePCR_verdict` | Screening results | PASS or FAIL |
| `ePCR_product_count` | Screening results | Number of PCR products found |
| `ePCR_products` | Screening results | Product coordinates and sizes |
| `ePCR_error` | Screening results | Error message (e.g., "wrong_genome: ...", "not_methylation_specific: ...") |
| `ePCR_elapsed_sec` | Screening results | Screening time in seconds |
| `dms_score` | DMS sub-scores | Composite DMS score (mean(A,B,C,D) x E) |
| `score_A_tg_near0` | DMS sub-scores | Sub-score A |
| `score_B_bg_near1` | DMS sub-scores | Sub-score B |
| `score_C_tg_consist` | DMS sub-scores | Sub-score C |
| `score_D_bg_consist` | DMS sub-scores | Sub-score D |
| `score_E_coverage` | DMS sub-scores | Sub-score E |
| `common_variant_score` | Primer annotations | Count of common SNPs overlapping primer sites |
| `snp_count` | Primer annotations | Total known SNPs overlapping primer sites |
| `segdup_overlap` | Primer annotations | Boolean: overlaps segmental duplication |
| `repeat_overlap` | Primer annotations | Boolean: overlaps repetitive element |
| `repeat_family` | Primer annotations | Repetitive element family name |

**Output sheets:**
- **All_Primers:** Original columns (minus removed) + ePCR columns + DMS sub-scores + annotation columns
- **Selected_Panel:** All passing primers with all columns
- **Best_Per_DMR, Summary, DMR_Blocks:** Preserved from original (removed columns dropped)
- **Data_Dictionary:** Generated by `data_dictionary.py`, with detailed explanation of every column in every sheet

**Output file:** `master_atlas_screened_v3.xlsx`

---

### 6b.10 data_dictionary.py

**File:** `scripts/data_dictionary.py`

**What it does:** Generates the Data_Dictionary sheet that is embedded in the final atlas xlsx. It documents every column in every sheet, with type, description, and source.

**Key documentation:**
- `dms_score`: Full A-E formula documentation (mean(A,B,C,D) x E, with each sub-score defined)
- `score_A_tg_near0` through `score_E_coverage`: Individual sub-score descriptions
- `ePCR_verdict`: PASS or FAIL
- `ePCR_error`: Documents all three FAIL subtypes (wrong_genome, not_methylation_specific, generic)
- `selected_for_panel`: "ALL passing primers are kept (no DMR deduplication, no cell-type cap)"
- `assay`: Retained (contains "M" or "U")
- Annotation columns: Detailed UCSC track descriptions

**Output:** 85 rows, with sequential column numbers per sheet.

---

### 6b.11 Output Files Summary

After running the full pipeline (`run_full_pipeline.sh`), the following files are produced in the working directory:

| File | Format | Description |
|------|--------|-------------|
| `master_atlas_screened_v3.xlsx` | Excel | Final merged atlas with all columns |
| `cpg_beta_values.csv` | CSV | Per-CpG beta values for all cell types, all selected DMR blocks |
| `dms_subscores.csv` | CSV | DMS A-E sub-scores for all DMR blocks |
| `primer_annotations.csv` | CSV | UCSC track annotations for passing primers |
| `results_v2_selected.csv` | CSV | Selected panel (all passing primers) |
| `results_v2.csv` | CSV | Full ePCR screening results (pre-existing, not re-generated) |

---

### 6b.12 Helper Scripts

**`prepare_genomes.py`** (repository root): Converts a hg19 FASTA file into per-chromosome numpy `.npy` arrays (uint8) for use by `bisearch_epcr.py` and `annotate_primers.py`. The genome conversion (bisulfite CT/CT_M) is done at runtime by the screening tools, so only the original genome needs to be stored.

**`validate_epcr.py`** (repository root): Validates ePCR screening results by spot-checking specific primer pairs against the genome. Useful for verifying that the hash-based indexing produces correct results.

**`requirements_epcr.txt`** (repository root): Lists the additional Python dependencies for the ePCR screening workflow (numpy, pandas, openpyxl, pyarrow).

**`docs/bisearch_epcr_README.md`**: Quick-start documentation for the standalone ePCR tool.

**`docs/config_patch.txt`**: Configuration patch notes for the ePCR workflow.

---

## Chapter 7 — Troubleshooting

### 7.1 Installation Problems

**Problem:** `brew: command not found`
- **Cause:** Homebrew not installed
- **Solution:** The install script should handle this. If not, install manually from https://brew.sh

**Problem:** `python3.11: command not found`
- **Cause:** Python 3.11 not installed
- **Solution:** `brew install python@3.11`

**Problem:** `ModuleNotFoundError: No module named 'methyl_panel'`
- **Cause:** Virtual environment not activated, or not in the repo folder
- **Solution:** Run `source .venv/bin/activate` from the repo folder

**Problem:** `ModuleNotFoundError: No module named 'primer3'`
- **Cause:** primer3-py not installed in the virtual environment
- **Solution:** `pip install primer3-py==2.3.0`

**Problem:** `brew install tabix` fails — "No available formula with the name 'tabix'"`
- **Cause:** Homebrew removed the standalone tabix formula; it's now bundled with htslib
- **Solution:** Use `brew install htslib` instead. The install script has been fixed to use htslib.

**Problem:** `zsh: permission denied: ./install_macos.sh` or `./download_data.sh`
- **Cause:** The script file is not executable. This happens on first clone, and also after `git checkout --` or `git pull` resets file permissions.
- **Solution:** Run `chmod +x install_macos.sh` (or `chmod +x download_data.sh`), then run the script again.

**Problem:** `git pull` fails with "Your local changes to the following files would be overwritten by merge"
- **Cause:** You modified a file locally (e.g. `chmod +x` changed the file, or you edited it) and the remote version also changed. The error message names the conflicting file(s).
- **Solution:** Discard your local copy of the conflicting file(s) and take the remote version. For each file listed in the error message:
  ```bash
  git checkout -- <filename>
  ```
  Then pull and restore permissions:
  ```bash
  git pull
  chmod +x install_macos.sh download_data.sh
  ```
  Example: if the error says `install_macos.sh` would be overwritten:
  ```bash
  git checkout -- install_macos.sh
  git pull
  chmod +x install_macos.sh download_data.sh
  ```
  `git checkout -- <filename>` replaces your local file with the last committed version. Then `git pull` gets the latest changes. Then `chmod +x` because `git checkout --` resets file permissions.

**Problem:** `pip install -e .` fails with `ModuleOrPackageNotFoundError: No file/folder found for package wgbs-tools`
- **Cause:** wgbstools is not a standard Python package — it has no `__init__.py` and its `pyproject.toml` has the `packages` line commented out. Poetry cannot find a package to install.
- **Solution:** This is expected. wgbstools is installed via a symlink, not pip. Run:
  ```bash
  cd wgbs_tools && python3 setup.py && cd ..
  ln -sf "$(pwd)/wgbs_tools/wgbstools" .venv/bin/wgbstools
  wgbstools --version
  ```
  The `wgbstools` executable is a symlink to `src/python/wgbs_tools.py`. When Python runs it, it resolves the symlink and adds `src/python/` to `sys.path`, which lets it find `find_markers.py` and other modules. The install script has been fixed to use this approach.

**Problem:** `pipeline.py: error: unrecognized arguments: --discover-dmrs`
- **Cause:** Your local code is out of date. The `--discover-dmrs` flag was added in v2.2.0. You need to `git pull` to get the latest code.
- **Solution:**
  ```bash
  git pull
  git log --oneline -1
  ```
  You should see a recent commit (not `6c47ed1` or earlier). If `git pull` says "already up to date" but `--discover-dmrs` is still not recognized, you may have a **nested directory** problem — the repo was cloned inside a folder with the same name. Check:
  ```bash
  pwd
  git remote -v
  grep -c "discover-dmrs" methyl_panel/pipeline.py
  ```
  If the `grep` returns `0`, you are in the wrong copy of the repo. Look for another `cell-specific-epigenetic-dpcr-assay-design-pipeline/` directory (possibly inside the current one) that has the updated code. The correct directory is the one where `grep -c "discover-dmrs" methyl_panel/pipeline.py` returns a number greater than `0`.

### 7.2 Data Problems

**Problem:** `FileNotFoundError: data/WBC_Panel_Top200_v7.9.xlsx`
- **Cause:** You're using Step 1 (load from Excel) but the DMR Excel file is not in the `data/` folder
- **Solution:** Use `--discover-dmrs` instead to discover DMRs from scratch from the beta files. Or, if you have a pre-computed Excel file, copy it to `data/`.

**Problem:** `wgbstools: command not found` (when using --discover-dmrs)
- **Cause:** wgbstools not installed or not on PATH
- **Solution:** Re-run `./install_macos.sh`. If that doesn't fix it, manually install:
  ```bash
  cd wgbs_tools && python3 setup.py && cd ..
  ln -sf "$(pwd)/wgbs_tools/wgbstools" .venv/bin/wgbstools
  ```
  Note: `pip install -e .` does NOT work for wgbstools — the package has no `__init__.py` and `pyproject.toml` has the packages line commented out. The symlink approach is the correct installation method.

**Problem:** `wgbstools find_markers failed with return code 1` with `ModuleNotFoundError: No module named 'scipy'`
- **Cause:** scipy is not installed in the venv. This can happen if the wgbstools `pip install -e .` step failed and rolled back other package installations.
- **Solution:** Install scipy directly in the venv:
  ```bash
  source .venv/bin/activate
  pip install scipy
  python -c "import scipy; print(scipy.__version__)"
  ```
  Then re-run the pipeline. The pipeline also checks for scipy before running find_markers and will tell you the exact command to run if it's missing.

**Problem:** `wgbstools find_markers failed with return code 1` (other causes)
- **Cause:** Missing beta files, blocks file, or groups file. Or insufficient memory.
- **Solution:** Ensure `data/beta_files/` has all 207 .beta files (run `./download_data.sh`). Ensure `data/GSE186458_blocks.s205.bed.gz` exists. Check the error message printed by find_markers.

**Problem:** `No markers found` or `0 blocks` from find_markers
- **Cause:** The target cell type has too few samples, or the delta_means threshold is too high, or beta files are missing for the target samples
- **Solution:** Check that the target cell type's beta files are present in `data/beta_files/`. The groups file is automatically filtered to only include samples with available beta files.

**Problem:** `No primers found` with `--max-blocks 10`
- **Cause:** The first 10 blocks by delta_means may have only 3–4 CpGs, below the `--min-cpg 5` threshold
- **Solution:** Use `--max-blocks 50` or more. The first blocks ranked by delta_means tend to have fewer CpGs; blocks with 5+ CpGs appear later in the ranking.

**Problem:** `samtools faidx failed` or `No sequence returned`
- **Cause:** hg19 genome FASTA not found, not indexed, or the chromosome is not in the genome
- **Solution:** Run `./download_data.sh` to download and index the genome. If only some chromosomes are available, blocks on missing chromosomes are automatically skipped.

**Problem:** hg19 genome download fails or `wgbstools init_genome` fails
- **Cause:** wgbstools `init_genome` may fail on some macOS configurations. The download script now downloads hg19 directly from UCSC instead.
- **Solution:** Re-run `./download_data.sh` after `git pull`. If it still fails, download manually:
  ```bash
  wget -O data/hg19/hg19.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
  samtools faidx data/hg19/hg19.fa.gz
  ```

**Problem:** dbSNP download fails with 404 Not Found
- **Cause:** NCBI reorganized the dbSNP FTP paths. The file is now `00-common_all.vcf.gz` (with `00-` prefix), not `common_all.vcf.gz`.
- **Solution:** Re-run `./download_data.sh` after `git pull`. The URL has been corrected.

**Problem:** `File truncated at line 1` when indexing genome
- **Cause:** The FASTA file was compressed with gzip instead of bgzip
- **Solution:** `gunzip -c file.fa.gz | bgzip -c > file.fa.bgz && mv file.fa.bgz file.fa.gz && samtools faidx file.fa.gz`

**Problem:** Beta file download fails for some samples
- **Cause:** NCBI FTP server may be slow or unavailable, or filenames are truncated
- **Solution:** Re-run `./download_data.sh` — it skips files that are already downloaded and uses directory listing to find truncated filenames

**Problem:** `zsh: command not found: #`
- **Cause:** You pasted a `#` comment line from this handbook
- **Solution:** Only paste lines that do NOT start with `#`. The `#` lines are explanations.

**Problem:** Terminal shows `quote>` and seems stuck
- **Cause:** You pasted a line containing an apostrophe (like "you're"), and zsh is waiting for a closing quote
- **Solution:** Press Ctrl-C to cancel, then re-run the command without the comment line

### 7.3 Pipeline Problems

**Problem:** `Tm values are REQUIRED`
- **Cause:** You didn't provide `--min-tm`, `--opt-tm`, `--max-tm`
- **Solution:** Always provide all three Tm values, or use `--settings` with a file that contains them

**Problem:** `No primers found for block X`
- **Cause:** The DMR block doesn't have enough CpGs, or the sequence is too short, or no primers passed Primer3's filters
- **Solution:** Try widening the Tm range, or increasing `--flank` (default 100), or lowering `--min-cpg`

**Problem:** `Skipping X (chrN): samtools faidx failed`
- **Cause:** The chromosome is not in your genome FASTA file
- **Solution:** This is expected if you don't have the full hg19 genome. The pipeline skips these blocks and continues with the ones it can fetch. Download the full genome with `./download_data.sh` to process all blocks.

**Problem:** ePCR screening shows `NULL` for all primers
- **Cause:** ePCR genome arrays not found in the specified directory
- **Solution:** Provide `--epcr-genome-dir` pointing to the directory with `.npy` files, or run `python methyl_panel/prepare_genomes.py --output-dir data/hg19_bs --genome data/hg19/hg19.fa` to generate them

**Problem:** SNP screening shows `NULL` for all primers
- **Cause:** dbSNP VCF not provided
- **Solution:** Provide `--dbsnp data/dbsnp_common.vcf.gz` (downloaded by `download_data.sh`)

**Problem:** Pipeline is slow
- **Cause:** Processing all blocks × 4 templates is computationally intensive
- **Solution:** Use `--cell-type` to process one cell type at a time, or `--max-blocks 10` for testing

**Problem:** All primers are M-assays (no U-assays)
- **Cause:** This is expected for hypomethylated DMRs. The unmethylated templates (SU/AU) are very AT-rich after bisulfite conversion, making it difficult for Primer3 to design primers within the Tm window.
- **Solution:** This is normal behavior, not an error. If you need U-assays, try widening the Tm range or increasing `--flank`.

---

## Chapter 8 — Repository Structure

```
cell-specific-epigenetic-dpcr-assay-design-pipeline/
├── methyl_panel/                    # Python package (main pipeline)
│   ├── __init__.py                  # Package init
│   ├── config.py                    # Primer3PlusConfig (161 params) + PipelineConfig
│   ├── pipeline.py                  # CLI entry point (steps 0-9, 6b=61)
│   ├── phase0_dmr_discovery.py      # Step 0: DMR discovery via wgbstools find_markers
│   ├── phase1_dmr_loader.py         # Step 1: Load DMR blocks from Excel (legacy)
│   ├── phase2_bisulfite_convert.py  # Step 2: Bisulfite convert to 6 strands
│   ├── phase3_primer3_design.py     # Step 3: Primer3 primer design
│   ├── bisearch_epcr_core.py        # Step 4: ePCR screening (hash-based numpy indexing)
│   ├── prepare_genomes.py           # Convert FASTA to .npy for ePCR
│   ├── phase5_structure.py          # Step 5: Hairpin MFE screening
│   ├── phase6_snp.py                # Step 6: dbSNP common variant screening
│   ├── phase7_dimer.py              # Step 7: DimerDetective primer-dimer prediction
│   ├── phase8_repeats.py            # Step 6b: UCSC repeat/segdup annotation
│   ├── output_xlsx.py               # Step 8: U-assays-style Excel output
│   └── output_pdf.py                # Step 9: U-assays-style PDF output
├── scripts/                         # Post-pipeline screening & annotation scripts
│   ├── compute_dms_subscores.py     # DMS A-E sub-score computation from beta files
│   ├── annotate_primers.py          # UCSC hg19 track primer annotation
│   ├── extract_cpg_beta.py          # Per-CpG beta value extraction
│   ├── merge_screening_results.py   # Merge ePCR + DMS + annotations into final atlas
│   ├── data_dictionary.py           # Automated column documentation for final atlas
│   ├── export_dmr_excel.py          # Export DMR regions + per-CpG methylation to Excel
│   └── generate_master_excel.py     # Generate master primer assay Excel
├── data/                            # Data files
│   ├── sample_manifest.csv          # 207 sample download manifest
│   ├── sample_groups.csv            # Full atlas cell type groups (82 types)
│   └── dmr_regions_v2.2.8_260714.xlsx  # DMR input for pipeline runs
├── docs/
│   ├── handbook.md                  # This document
│   ├── bisearch_epcr_README.md      # Quick-start for standalone ePCR tool
│   ├── data_setup.md                # Data download and setup instructions
│   └── figures/                     # Pipeline figures (PNG + SVG)
├── tests/                           # Unit tests
├── bisearch_epcr.py                 # Standalone ePCR screening tool (Chapter 6b)
├── validate_epcr.py                 # Validate ePCR screening results
├── run_full_pipeline.sh             # Consolidated post-pipeline runner
├── annotate_primers.sh              # Post-pipeline annotation shell script
├── download_data.sh                 # Data download script (~10 GB)
├── install_macos.sh                 # macOS installation script
├── setup_local.sh                   # Local setup script
├── requirements.txt                 # Python dependencies
├── README.md                        # Project overview
└── .gitignore
```

---

## Chapter 9 — Data Source

Loyfer N et al. (2023). A human methylation atlas for normal and cancer cell types. *Nature* 617, 690–697.

Dataset: [GSE186458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186458) — 207 samples, 82 cell types, hg19 (GRCh37).

---

## Appendix A — Glossary

| Term | Definition |
|------|-----------|
| **DMR** | Differentially Methylated Region — a genomic block where methylation differs between cell types |
| **MSP** | Methylation-Specific PCR — PCR using primers designed on bisulfite-converted DNA |
| **CpG** | A cytosine followed by a guanine (5'→3') — the site of DNA methylation in mammals |
| **MFE** | Minimum Free Energy — the stability of a secondary structure (more negative = more stable) |
| **ΔG** | Gibbs free energy — thermodynamic stability (negative = stable, positive = unstable) |
| **MAF** | Minor Allele Frequency — how common a SNP is in the population |
| **Tm** | Melting temperature — the temperature at which half of DNA duplexes dissociate |
| **SM/AM** | Sense/Antisense Methylated — bisulfite-converted strand with CpG C's preserved |
| **SU/AU** | Sense/Antisense Unmethylated — bisulfite-converted strand with all C→T |
| **DimerDetective** | A primer-dimer prediction method based on 3' end thermodynamic stability |
| **U-assays** | The output format (XLSX + PDF) matching the reference assay tracking format |
| **wgbstools** | Tools for analyzing whole-genome bisulfite sequencing data |
| **Beta file** | Binary methylation data file (2 bytes per CpG: methylated count + total count) |
| **DMS score** | Differential Methylation Suitability — per-CpG quality score (0-1) quantifying how suitable a DMR is for bisulfite PCR primer design. Computed from five sub-scores (A-E): A = 1 - max(target CpG mean methylation), B = min(background CpG mean methylation), C = max(0, 1 - 2*mean(target CpG std)), D = max(0, 1 - 2*mean(background CpG std)), E = mean(fraction CpGs with adequate coverage). Final score = mean(A,B,C,D) x E. Threshold >= 0.6 = good candidate. Formerly called "cleanliness score"; renamed in v3.1.0. Source code function is still named `cleanliness_score()` for backward compatibility. |
| **bgzip** | Block-compressed gzip — required by samtools faidx (not regular gzip) |

---

## Appendix B — Quick Start (Copy-Paste)

```bash
git clone https://github.com/ScottAilliet/cell-specific-epigenetic-dpcr-assay-design-pipeline.git
cd cell-specific-epigenetic-dpcr-assay-design-pipeline
chmod +x install_macos.sh
./install_macos.sh
```

```bash
chmod +x download_data.sh
./download_data.sh
```

```bash
source .venv/bin/activate
```

```bash
# Discover DMRs from scratch and design primers for monocytes
python -m methyl_panel.pipeline --steps 0,2,3,5,7,8,9 \
    --discover-dmrs --cell-type MONO \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/
```

```bash
# Same for CD8 T cells
python -m methyl_panel.pipeline --steps 0,2,3,5,7,8,9 \
    --discover-dmrs --cell-type CD8T \
    --genome data/hg19/hg19.fa.gz \
    --min-tm 58 --opt-tm 60 --max-tm 62 \
    --output-dir results/CD8T/
```

```bash
# All 7 immune cell types — hypomethylated markers (blood-only, strict 0.70 filter)
for CT in MONO BCELL NK GRAN CD3T CD8T CD4T; do
    python -m methyl_panel.pipeline --steps 0,2,3,5,7,8,9 \
        --discover-dmrs --cell-type $CT \
        --genome data/hg19/hg19.fa.gz \
        --min-tm 58 --opt-tm 60 --max-tm 62 \
        --output-dir results/$CT/
done
```

```bash
# Hypermethylated markers (target methylated, background unmethylated)
for CT in MONO BCELL NK GRAN CD3T CD8T CD4T; do
    python -m methyl_panel.pipeline --steps 0,2,3,5,7,8,9 \
        --discover-dmrs --cell-type $CT \
        --genome data/hg19/hg19.fa.gz \
        --min-tm 58 --opt-tm 60 --max-tm 62 \
        --only-hyper \
        --output-dir results_hyper/$CT/
done
```

```bash
open results/MONO/primer_assays.xlsx
open results/MONO/primer_assays.pdf
```

---

## Appendix C — Change Log

### v3.1.0 (2026-08-25)

**Terminology:**
- Renamed "cleanliness score" to "DMS score" (Differential Methylation Suitability) throughout the handbook and output files. The source code function is still named `cleanliness_score()` for backward compatibility.
- Fixed DMS formula documentation to match the actual code: sub-scores A and B use max()/min() (not mean()), and sub-scores C and D use a factor-of-2 scaling (not raw 1-std).

**New standalone ePCR tool (`bisearch_epcr.py`, repository root):**
- Hash-based numpy indexing with sorted arrays (no bowtie2 or samtools dependency)
- 16-mer seed extraction, mismatch string "1233333333333333", max 3 mismatches
- Searches both strands (Watson/Crick) and both methylation states (CT/CT_M)
- Products deduplicated per (chrom, strand, fp_pos, rp_pos); a locus on both states counts once
- 3-way template-aware verdict: PASS / FAIL wrong_genome / FAIL not_methylation_specific
- Template-to-state mapping: SM/AM -> methylated (CT_M), SU/AU -> unmethylated (CT)
- Selection: keeps ALL passing primers (no DMR deduplication, no cell-type cap)
- CLI: `screen` subcommand (single pair or batch from CSV), `select` subcommand

**New scripts:**
- `scripts/compute_dms_subscores.py`: Computes DMS A-E sub-scores for all DMR blocks from beta files
- `scripts/annotate_primers.py`: Downloads 4 UCSC hg19 tracks (snp151Common, snp151, rmsk, genomicSuperDups) and annotates passing primers with SNP, repeat, and segmental duplication overlaps
- `scripts/extract_cpg_beta.py`: Extracts per-CpG beta values for all cell types and samples for every selected DMR block; outputs long-format CSV
- `scripts/merge_screening_results.py`: Merges ePCR results, DMS sub-scores, and primer annotations into `master_atlas_screened_v3.xlsx`; removes direction, cleanliness_score, and bowtie2-specific columns; generates Data_Dictionary sheet
- `scripts/data_dictionary.py`: Automated column documentation for the final atlas (85 rows)
- `run_full_pipeline.sh` (repo root): Consolidated 7-step post-pipeline runner (git pull, copy scripts, re-select, compute DMS, annotate, extract CpG beta, merge). Does NOT re-run ePCR screening.
- `prepare_genomes.py` (repo root): Converts hg19 FASTA to per-chromosome numpy .npy arrays for ePCR
- `validate_epcr.py` (repo root): Validates ePCR screening results by spot-checking primer pairs
- `methyl_panel/bisearch_epcr_core.py`: Core ePCR engine (hash-based numpy indexing)
- `methyl_panel/phase4_bisearch_epcr.py`: Pipeline integration for ePCR screening

**Column changes in `master_atlas_screened_v3.xlsx`:**
- Removed: `direction` (unreliable), `cleanliness_score` (renamed to `dms_score`), and 8 bowtie2-specific columns (sense_meth_mismatch_score, sense_unmeth_mismatch_score, anti_meth_mismatch_score, anti_unmeth_mismatch_score, bowtie_passes_filter, bowtie_intended_genome, mismatch_profile, mismatch_detail)
- Retained: `assay` column (contains "M" or "U"), distinct from the removed `direction` column
- Added: `ePCR_verdict`, `ePCR_product_count`, `ePCR_products`, `ePCR_error`, `ePCR_elapsed_sec`, `dms_score`, `score_A_tg_near0` through `score_E_coverage`, `common_variant_score`, `snp_count`, `segdup_overlap`, `repeat_overlap`, `repeat_family`

**Handbook updates:**
- Added Chapter 6b: Standalone ePCR Screening and Annotation (12 sections covering the tool, verdict logic, selection logic, all scripts, output files, and helper scripts)
- Renamed "cleanliness score" to "DMS score" in current documentation (Sections 3.2, 4.2)
- Added note in v2.2.7 changelog that "cleanliness score" was renamed to "DMS score" in v3.1.0
- Updated Chapter 8 (Repository Structure) to include new files
- Added DMS score entry to Glossary (Appendix A)
- Updated version header to v3.1.0

### v2.2.7 (2026-07-14)

- `phase0_dmr_discovery.py`: Fixed per-subgroup filter bypass. Blocks with empty `bg_subgroup_meth` (no subgroup data due to low WGBS coverage) were silently passing the filter. Blocks with partial `bg_subgroup_meth` (some subgroups missing) were checked only on present subgroups. Now both are rejected — missing subgroup data means we cannot verify that subgroup, so the safe choice is to reject. Rejection report now categorizes as "no subgroup data", "partial subgroup data", or threshold violations.
- `phase0_dmr_discovery.py`: Added cleanliness score filter (CS >= 0.6) at Step 0e1, before the per-subgroup filter. Blocks with poor per-CpG data quality (low coverage, high variability, CS=0) are now rejected. Previously these were computed but never filtered.
  - *Note: "cleanliness score" was renamed to "DMS score" (Differential Methylation Suitability) in v3.1.0. The formula and threshold are unchanged.*
- `phase0_dmr_discovery.py`: Added sex chromosome exclusion (chrX, chrY, chrM) at Step 0d2, before per-CpG extraction. chrX has copy number differences between sexes and X-inactivation methylation differences; chrY is male-only; chrM has different copy number biology. All make dPCR quantification unreliable.
- `export_dmr_excel.py`: Fixed hyper sheet naming. Sheets are now named `{CT}_hypo` and `{CT}_hyper` instead of `{CT}` and `{CT}1`. Added Direction column (`U`/`M`) to block-level sheets, summary sheet, and Per-CpG detail sheet.

### v2.2.6 (2026-07-14)

- `pipeline.py`, `phase0_dmr_discovery.py`, `config.py`: Default `--top-markers` changed from 300 to 1000. With 300 candidates and the strict 0.70 per-subgroup filter, only 9 MONO and 6 NK markers survived — too few for a panel. With 1000, the filter has 3x more correctly-ranked candidates to work with.
- Deprecated `--use-full-atlas` for hypomethylated markers. Testing showed the full atlas background finds more raw candidates but ranks them by atlas delta (difference from liver/brain/colon), not blood delta. The top 300 by atlas delta are the regions most different from non-blood tissues — the worst candidates for blood subgroup separation. Blood-only find_markers ranks by blood delta, which is the correct metric for what the per-subgroup filter checks. The flag remains available but is not in recommended commands. Hypermethylated markers (`--only-hyper`) remain recommended — they search a different set of regions.

### v2.2.5 (2026-07-14)

- `phase0_dmr_discovery.py`: Added hybrid atlas approach. New `generate_groups_file_full_atlas()` creates a groups file with all 207 atlas samples as background (minus BG_EXCLUDE groups). `discover_dmrs()` now generates two groups files: a blood-only file for per-CpG extraction and per-subgroup filter, and a full-atlas file for find_markers. Activated by `use_full_atlas=True` parameter and `--use-full-atlas` CLI flag. Finds more candidate DMRs for closely related T-cell subtypes (CD4T, CD8T) while maintaining blood-only quality control.
- `phase0_dmr_discovery.py`: Added hypermethylated marker support. `run_find_markers()` now accepts `only_hyper=True` which passes `--only_hyper` to wgbstools find_markers (swaps target/background internally). The per-subgroup filter is inverted for hyper markers: rejects if any subgroup > (1 - min_bg_subgroup_meth) instead of < min_bg_subgroup_meth. `_cleanliness_score()` updated with `only_hyper` parameter: component A = min(target mean) for hyper (target near 1), B = 1 - max(background mean) for hyper (background near 0). `direction` field ("M" for hyper, "U" for hypo) added to DMRBlock and JSON output. Activated by `--only-hyper` CLI flag.
- `phase0_dmr_discovery.py`, `pipeline.py`: Raised `--min-bg-subgroup-meth` default from 0.50 to 0.70. Every background blood cell type must now be at least 70% methylated (hypo) or at most 30% methylated (hyper) at a DMR region. Ensures strong binary discrimination at the DMR region level rather than relying on primer design to compensate for marginal separation.
- `export_dmr_excel.py`: Updated for hypermethylated markers. `format_block_row()` detects hyper from `direction` field and shows "Max subgroup methylation" column instead of "Min subgroup methylation". The "Worst background subgroup" column shows the subgroup with the highest methylation for hyper markers (instead of lowest for hypo). Summary sheet includes both Min and Max columns.
- `pipeline.py`: Added `--use-full-atlas` and `--only-hyper` CLI arguments, passed through to `discover_dmrs()`.

### v2.2.4 (2026-07-14)

- `phase0_dmr_discovery.py`: Fixed subgroup-aware background subsampling. `extract_percpg_methylation_with_indices()` used `background_betas[:max_bg_samples]` which silently dropped entire subgroups (Blood-B-Mem was entirely dropped for MONO/NK/GRAN because both samples were at positions 31-33 of 33). Added `_subsample_background()` function that groups background samples by atlas subgroup, guarantees at least 1 sample per subgroup, then distributes remaining slots proportionally using the largest-remainder method. Handles edge cases: more groups than max (takes 1 from largest groups), no group map (falls back to simple truncation), empty input.
- `export_dmr_excel.py`: Fixed per-cell-type sheets showing all 14 subgroup columns. Each sheet now only includes columns for subgroups present in that cell type's `bg_subgroup_meth` data. The "All blocks summary" sheet retains all subgroup columns for cross-cell-type comparison. Previously, CD4T had 5 empty columns (4 CD4 target groups + excluded Blood-T-CD3), CD3T had 9 empty columns — making the data look misaligned.

### v2.2.3 (2026-07-14)

- `phase0_dmr_discovery.py`: Added per-subgroup background check. After per-CpG extraction, the pipeline computes mean methylation for each background blood cell type subgroup independently (e.g. Blood-T-CD8, Blood-B, Blood-NK). Blocks where any subgroup has mean methylation below `--min-bg-subgroup-meth` (default 0.50) are rejected. This catches the key problem with CD4T/CD8T: a region where the average background looks fine (0.88) but CD8 T cells are partially unmethylated (0.55) — which would cause false positives in real blood. Added `_match_beta_files_with_groups()` to map background beta files to their original atlas group names. Added `bg_subgroup_meth` field to DMRBlock and JSON output.
- `phase0_dmr_discovery.py`: Added `BG_EXCLUDE` — groups excluded from background in addition to targets. Blood-T-CD3 (pan-T cells containing both CD4+ and CD8+) is excluded from CD4T and CD8T backgrounds. Without this, the per-subgroup filter rejected 68/73 CD4T markers because T-CD3 contains CD4+ T cells that are unmethylated at CD4-specific regions. `generate_groups_file()` now filters out excluded groups before assigning target/background labels.
- `phase0_dmr_discovery.py`, `pipeline.py`: Added `--unmeth-mean-thresh` (default 0.15) and `--meth-mean-thresh` (default 0.65) — these were previously disabled (find_markers defaults of 1.0 and 0.0), allowing mediocre markers with target=0.3/bg=0.6 to pass. `delta_means` stays at 0.3 (was briefly 0.4 but reverted — too few markers for CD4T/CD8T).
- `phase1_dmr_loader.py`: Added `bg_subgroup_meth: Dict[str, float]` field to DMRBlock dataclass.
- `export_dmr_excel.py`: Added per-subgroup methylation columns (one per blood cell type, e.g. "BG: Blood-T-CD8"), "Min subgroup methylation", and "Worst background subgroup" columns to block-level sheets. First pass collects all unique subgroup names across all input files to build consistent columns.

### v2.2.2 (2026-07-13)

- `export_dmr_excel.py`: New script that exports DMR regions with per-CpG methylation to Excel. Reads `dmr_blocks.json` files from pipeline output directories and produces a multi-sheet workbook: one sheet per cell type (block-level summary), an "All blocks summary" sheet, and a "Per-CpG detail" sheet with one row per CpG site. Use `--top N` to control how many regions per cell type (default 200).
- `pipeline.py`, `phase0_dmr_discovery.py`, `config.py`: Default `--top-markers` changed from 200 to 300. find_markers now returns the top 300 DMR candidates per cell type (was 200). Override with `--top-markers N` for any other number.
- `phase0_dmr_discovery.py`: Added scipy pre-check in `run_find_markers()`. Before invoking wgbstools, the pipeline runs `import scipy` in the venv Python. If scipy is missing, it raises a `RuntimeError` with the exact install command (`<venv_python> -m pip install scipy`) instead of letting wgbstools crash with `ModuleNotFoundError: No module named 'scipy'`.
- `install_macos.sh`: Added scipy re-verification after the wgbstools install step. If scipy was rolled back by a failed `pip install -e .` (as happened on macOS), the script automatically reinstalls it.

### v2.2.1 (2026-07-10)

- `phase0_dmr_discovery.py`: Background restricted to blood cell types only. Previously, `generate_groups_file()` used all 207 atlas samples as background (including liver, brain, colon, etc.). Now only the 36 blood samples (14 `Blood-*` groups) are used. This is biologically correct for a WBC assay — DMRs must distinguish the target from other blood cells, not from non-blood tissues. Added `BLOOD_GROUPS` constant listing all 14 blood cell type group names.
- Removed multiplex Step 10 (`phase10_multiplex.py`) and all associated CLI arguments (`--multiplex-dirs`, `--cross-dimer-cutoff`, `--tm-tolerance`, `--min-amplicon-diff`).

### v2.2.0 (2026-07-10)

**New feature — DMR discovery from scratch:**
- `phase0_dmr_discovery.py`: New module implementing Step 0. Generates per-cell-type groups files (target vs background), runs `wgbstools find_markers` on raw beta files, parses BED output, extracts per-CpG methylation from beta files, computes 5-component cleanliness scores, and saves `dmr_blocks.json` in the same format as Step 1.
- `pipeline.py`: Added `step0_discover_dmrs()` and new CLI arguments (`--discover-dmrs`, `--beta-dir`, `--blocks-file`, `--groups-csv`, `--wgbstools-path`, `--threads`, `--top-markers`, `--max-bg-samples`, `--skip-find-markers`). When `--discover-dmrs` is set, Step 0 replaces Step 1.
- `install_macos.sh`: Fixed wgbstools installation — uses symlink to `.venv/bin/wgbstools` (not `pip install -e .`, which fails because the package has no `__init__.py`). Added `scipy` dependency. Removed unnecessary `init_genome` call.
- Cell type target mapping: 7 immune cell types mapped to atlas groups (MONO→Blood-Monocytes, BCELL→Blood-B+Blood-B-Mem, NK→Blood-NK, GRAN→Blood-Granulocytes, CD3T→9 T cell subtypes, CD4T→4 CD4 subtypes, CD8T→4 CD8 subtypes).
- Backward compatible: Step 1 (load from Excel) still works via `--dmr-xlsx`.

**Verified by sandbox testing:** Step 0 runs end-to-end on chr22 with 10 beta files (3 Monocytes + 7 background). find_markers finds 1,873 markers, outputs top 200. Per-CpG methylation extracted (target ~0.01, background ~0.93). 197/200 blocks have cleanliness score ≥ 0.6. Full pipeline (steps 0,2,3,5,7,8,9, --max-blocks 50) produces 30 primer pairs, XLSX (30 rows, 27 cols), PDF (76 KB).

### v2.1.2 (2026-07-09)

- `download_data.sh`: Fixed dbSNP URL — NCBI renamed `common_all.vcf.gz` to `00-common_all.vcf.gz`. Now downloads pre-built `.tbi` tabix index from NCBI instead of building it locally.
- `download_data.sh`: Fixed hg19 genome download — replaced wgbstools `init_genome` (which fails silently on some macOS configs) with direct download from UCSC (`hg19.fa.gz`, ~900 MB). Includes automatic bgzip recompression if needed.
- Handbook: Added troubleshooting for dbSNP 404 and hg19/wgbstools init_genome failure.

### v2.1.1 (2026-07-09)

- `download_data.sh`: Fixed wget consuming stdin from the manifest file, which killed the download loop after the first iteration. The manifest is now read on file descriptor 3 and wget stdin is redirected from /dev/null.
- Handbook: Added troubleshooting for `git pull` conflict ("Your local changes would be overwritten") and `chmod +x` needed after `git checkout --`.

### v2.1.0 (2026-07-09)

**Code fixes:**
- `phase1_dmr_loader.py`: Rewritten to auto-detect and load v7.9 Excel format (one sheet per cell type). Added cell type ID mapping (Mono→MONO, B→BCELL, etc.). Fallback CpG position distribution when `block_cpg_coords` is missing.
- `config.py`: Added missing `PipelineConfig` fields (`max_primers_per_block`, `min_cpg_per_primer`, `min_cpg_pair_total`) that caused step 3 to crash.
- `pipeline.py`: Added graceful skip for blocks on chromosomes not in the genome. Updated default DMR and genome paths. Updated `--cell-type` help text.
- `download_data.sh`: Removed `set -e` (was killing script on first download failure). Use FTP directory listing to handle NCBI filename truncation (74/207 files affected). Added progress counts.
- `install_macos.sh`: Changed `tabix` to `htslib` (Homebrew removed standalone tabix formula).

**Handbook fixes:**
- Corrected DMR file path to `WBC_Panel_Top200_v7.9.xlsx`
- Added cell type ID ↔ sheet name mapping table
- Added "Running Without Bowtie v1 and dbSNP" section
- Added "Running All 7 Cell Types" section
- Added "Testing with a Small Subset" section
- Added troubleshooting entries for: permission denied, `#` comment lines, `quote>` stuck prompt, bgzip vs gzip, NCBI filename truncation, missing chromosomes, M-assays only
- Added "What Changed from v1" section
- Added "M-assays vs U-assays" explanation
- Preserved old handbook as `handbook_v1.md`

**Verified by stress testing:** All 7 cell types (MONO, BCELL, NK, GRAN, CD3T, CD8T, CD4T) run end-to-end through steps 1,2,3,5,7,8,9 with 20 blocks each. Total: 123 primer pairs, 7 XLSX files (27 columns each), 7 PDF files.
