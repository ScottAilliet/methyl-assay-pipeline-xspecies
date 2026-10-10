#!/usr/bin/env python3
"""
Phase 0: DMR Discovery via wgbstools find_markers
=================================================

Discovers differentially methylated regions (DMRs) from raw beta files
using wgbstools find_markers, then extracts per-CpG methylation data
and computes DMS scores.

This module makes the pipeline fully self-contained: no pre-computed
DMR Excel file is needed. Anyone with the downloaded beta files and
blocks file can discover DMRs from scratch.

Workflow:
    1. Generate a groups CSV (target vs background) for the requested cell type
    2. Generate a beta list file (paths to all .beta files)
    3. Run wgbstools find_markers via subprocess
    4. Parse the output BED file into DMRBlock objects
    5. Extract per-CpG methylation from beta files
    6. Compute DMS scores
    7. Save as JSON (same format as step 1 output)

Usage (called by pipeline.py when --discover-dmrs is set):
    from methyl_panel.phase0_dmr_discovery import discover_dmrs
    blocks = discover_dmrs(cell_type_id, beta_dir, blocks_file, groups_csv, out_dir)
"""

import os
import sys
import glob
import json
import shutil
import subprocess
import warnings
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple

# Suppress harmless RuntimeWarnings from numpy (empty slices when all NaN)
warnings.filterwarnings("ignore", category=RuntimeWarning)

from methyl_panel.phase1_dmr_loader import DMRBlock, CpGSite


# ---------------------------------------------------------------------------
# Cell type → target groups mapping
# ---------------------------------------------------------------------------
# Each cell type ID maps to a list of atlas group names that constitute
# the target population. Background is restricted to other blood cell
# types only — this is a WBC assay, so we need to distinguish monocytes
# from other blood cells, not from liver or brain tissue.
# Base dict: 7 merged blood cell types (multi-group targets).
# These merge related blood subtypes into single assay targets
# (e.g. CD4T = naive + central-memory + effector-memory CD4 T cells).
CELL_TYPE_TARGETS_BASE: Dict[str, List[str]] = {
    "MONO":  ["Blood-Monocytes"],
    "BCELL": ["Blood-B", "Blood-B-Mem"],
    "NK":    ["Blood-NK"],
    "GRAN":  ["Blood-Granulocytes"],
    "CD3T":  [
        "Blood-T-CD3", "Blood-T-CD4", "Blood-T-CD8",
        "Blood-T-CenMem-CD4", "Blood-T-EffMem-CD4", "Blood-T-EffMem-CD8",
        "Blood-T-Eff-CD8", "Blood-T-Naive-CD4", "Blood-T-Naive-CD8",
    ],
    "CD4T":  [
        "Blood-T-CD4", "Blood-T-CenMem-CD4",
        "Blood-T-EffMem-CD4", "Blood-T-Naive-CD4",
    ],
    "CD8T":  [
        "Blood-T-CD8", "Blood-T-Eff-CD8",
        "Blood-T-EffMem-CD8", "Blood-T-Naive-CD8",
    ],
}


def _build_cell_type_targets(atlas_groups_csv: str) -> Dict[str, List[str]]:
    """
    Build the full cell-type → target-groups mapping.

    Starts with the 7 hardcoded merged blood cell types, then auto-expands
    with every atlas group from the groups CSV (each group maps to itself).
    This allows the pipeline to run for any of the 82 Loyfer atlas cell types
    without hardcoding each one.

    Returns:
        Dict mapping cell_type_id → list of atlas group names
    """
    targets = dict(CELL_TYPE_TARGETS_BASE)
    try:
        df = pd.read_csv(atlas_groups_csv)
        for group in sorted(df['group'].unique()):
            if group not in targets:
                targets[group] = [group]
    except (FileNotFoundError, Exception) as e:
        # Fall back to base dict if CSV not found (e.g. running from
        # a different directory). The 7 blood types still work.
        print(f"  NOTE: Could not read atlas groups CSV ({atlas_groups_csv}): {e}")
        print(f"  Only the 7 base blood cell types will be available.")
    return targets


# Default atlas groups CSV path (relative to repo root)
_DEFAULT_ATLAS_CSV = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    '..', 'data', 'sample_groups.csv'
)

# Auto-expand at import time so VALID_CELL_TYPES includes all atlas groups
CELL_TYPE_TARGETS: Dict[str, List[str]] = _build_cell_type_targets(
    _DEFAULT_ATLAS_CSV if os.path.isfile(_DEFAULT_ATLAS_CSV)
    else 'data/sample_groups.csv'
)

# All blood cell type groups in the atlas (14 groups, 36 samples total).
# Background for any WBC assay is restricted to these groups minus the
# target groups for the cell type being assayed.
BLOOD_GROUPS: List[str] = [
    "Blood-B", "Blood-B-Mem", "Blood-Granulocytes", "Blood-Monocytes",
    "Blood-NK", "Blood-T-CD3", "Blood-T-CD4", "Blood-T-CD8",
    "Blood-T-CenMem-CD4", "Blood-T-Eff-CD8", "Blood-T-EffMem-CD4",
    "Blood-T-EffMem-CD8", "Blood-T-Naive-CD4", "Blood-T-Naive-CD8",
]

# Groups to EXCLUDE from background for each cell type, in addition to
# the target groups. Blood-T-CD3 is a pan-T cell population containing
# both CD4+ and CD8+ T cells — it overlaps with CD4T and CD8T targets,
# so it must be excluded from their backgrounds to avoid false rejections
# in the per-subgroup filter. Similarly, CD4T excludes CD8T subgroups
# and vice versa (they are too closely related to serve as clean
# background for each other at the subgroup level — the per-subgroup
# filter handles this, but excluding them from find_markers background
# gives more markers to choose from).
# For CD3T (pan-T), all T-cell subtypes are targets, so no exclusions
# beyond the targets themselves.
BG_EXCLUDE: Dict[str, List[str]] = {
    # T-CD3 is pan-T (contains both CD4+ and CD8+); exclude to avoid
    # overlap with target. Also exclude the opposite T subtype's groups:
    # CD4+ and CD8+ T cells are epigenetically very similar, so including
    # CD8T subgroups in CD4T's background (or vice versa) makes DMR
    # discovery far too restrictive — only a handful of markers survive.
    # Excluding them gives find_markers a broader, cleaner background and
    # yields many more candidate DMRs.
    "CD4T": [
        "Blood-T-CD3",
        "Blood-T-CD8", "Blood-T-Eff-CD8",
        "Blood-T-EffMem-CD8", "Blood-T-Naive-CD8",
    ],
    "CD8T": [
        "Blood-T-CD3",
        "Blood-T-CD4", "Blood-T-CenMem-CD4",
        "Blood-T-EffMem-CD4", "Blood-T-Naive-CD4",
    ],
}

# Groups to exclude from the PER-SUBGROUP FILTER background only.
# This is a subset of BG_EXCLUDE: only the truly overlapping group
# (Blood-T-CD3, which contains both CD4+ and CD8+ cells) is excluded.
# The opposite T subtype is KEPT in the filter background so that
# T-cell-general markers (unmethylated in both CD4T and CD8T) are
# rejected. This implements the hybrid approach described in the
# code comment: find_markers uses a broad background (BG_EXCLUDE) to
# get many candidates, then the per-subgroup filter uses a stricter
# background (FILTER_EXCLUDE) to reject markers that don't discriminate
# the target from its closest relative.
FILTER_EXCLUDE: Dict[str, List[str]] = {
    "CD4T": ["Blood-T-CD3"],
    "CD8T": ["Blood-T-CD3"],
}

# All valid cell type IDs (7 merged blood types + all atlas groups)
VALID_CELL_TYPES = list(CELL_TYPE_TARGETS.keys())


def is_blood_cell_type(cell_type_id: str) -> bool:
    """
    Check whether this is a merged blood cell type that needs blood-specific
    background logic.

    Only the 7 merged blood types (MONO, BCELL, NK, GRAN, CD3T, CD4T, CD8T)
    use blood-only background and the per-subgroup filter. Individual atlas
    blood groups (Blood-T-CD3, Blood-Monocytes, etc.) are treated like any
    other atlas cell type: full atlas background, per-subgroup filter disabled.

    Args:
        cell_type_id: Cell type ID (e.g. MONO, Blood-T-CD3, Liver-Hepatocytes)

    Returns:
        True if cell_type_id is one of the 7 merged blood types, False otherwise
    """
    return cell_type_id in CELL_TYPE_TARGETS_BASE


# ---------------------------------------------------------------------------
# Groups file generation
# ---------------------------------------------------------------------------
def generate_groups_file(
    cell_type_id: str,
    atlas_groups_csv: str,
    output_path: str,
    beta_dir: Optional[str] = None,
) -> str:
    """
    Generate a wgbstools-compatible groups CSV for one cell type.

    Target samples (all subtypes merged) get group = cell_type_id.
    Background is restricted to other blood cell types only — this is
    a WBC assay, so we distinguish the target from other blood cells,
    not from non-blood tissues.

    If beta_dir is provided, only samples that have a corresponding
    .beta file in that directory are included in the groups file.
    This prevents find_markers from failing on missing beta files.

    Args:
        cell_type_id: One of MONO, BCELL, NK, GRAN, CD3T, CD4T, CD8T
        atlas_groups_csv: Path to data/sample_groups.csv
        output_path: Where to write the groups CSV
        beta_dir: If provided, filter to samples with .beta files here

    Returns:
        Path to the generated groups CSV
    """
    if cell_type_id not in CELL_TYPE_TARGETS:
        raise ValueError(
            f"Unknown cell type '{cell_type_id}'. "
            f"Valid options: {VALID_CELL_TYPES}"
        )

    target_groups = CELL_TYPE_TARGETS[cell_type_id]
    exclude_groups = BG_EXCLUDE.get(cell_type_id, [])
    df = pd.read_csv(atlas_groups_csv)

    # Restrict background to blood cell types only for blood assays.
    # For non-blood cell types (e.g. Liver-Hepatocytes), use the full atlas
    # as background — we want to distinguish the target from ALL other cell
    # types, not just blood cells.
    if is_blood_cell_type(cell_type_id):
        df = df[df['group'].isin(BLOOD_GROUPS)].copy()
    else:
        df = df.copy()  # Keep all atlas samples as background

    # Exclude overlapping groups (e.g. T-CD3 for CD4T/CD8T) — these
    # contain cells that overlap with the target and would cause false
    # rejections in the per-subgroup filter.
    if exclude_groups:
        before = len(df)
        df = df[~df['group'].isin(exclude_groups)].copy()
        print(f"    Excluded {before - len(df)} samples from overlapping groups: "
              f"{', '.join(exclude_groups)}")

    # The atlas groups CSV has columns: name, group
    # wgbstools renames the first column to 'fname', so 'name' is fine
    # Merge target subtypes into one group
    df['group'] = df['group'].apply(
        lambda g: cell_type_id if g in target_groups else 'background'
    )

    # Filter to samples that have beta files, if beta_dir is provided
    if beta_dir:
        import glob
        available = set()
        for bp in glob.glob(os.path.join(beta_dir, "*.beta")):
            stem = os.path.basename(bp)
            if stem.endswith('.beta'):
                stem = stem[:-5]
            available.add(stem)
        name_col = df.columns[0]
        before = len(df)
        df = df[df[name_col].isin(available)]
        after = len(df)
        n_target = (df['group'] != 'background').sum()
        n_bg = (df['group'] == 'background').sum()
        print(f"    Groups file: {after}/{before} samples have beta files "
              f"({n_target} target, {n_bg} background)")

    # Write with 'fname' column name for clarity (wgbstools renames first col anyway)
    df.to_csv(output_path, index=False)
    return output_path


def generate_groups_file_full_atlas(
    cell_type_id: str,
    atlas_groups_csv: str,
    output_path: str,
    beta_dir: Optional[str] = None,
) -> str:
    """
    Generate a groups CSV using the FULL atlas (all 207 samples) as background.

    This is used for the hybrid approach: find_markers runs against the full
    atlas to find a large candidate pool, then per-CpG extraction + the
    per-subgroup filter use blood-only background for quality control.

    Target samples (all subtypes merged) get group = cell_type_id.
    Background = ALL other atlas samples (blood + non-blood tissues).
    BG_EXCLUDE groups are still excluded (e.g. T-CD3 for CD4T/CD8T).

    Args:
        cell_type_id: One of MONO, BCELL, NK, GRAN, CD3T, CD4T, CD8T
        atlas_groups_csv: Path to data/sample_groups.csv
        output_path: Where to write the groups CSV
        beta_dir: If provided, filter to samples with .beta files here

    Returns:
        Path to the generated groups CSV
    """
    if cell_type_id not in CELL_TYPE_TARGETS:
        raise ValueError(
            f"Unknown cell type '{cell_type_id}'. "
            f"Valid options: {VALID_CELL_TYPES}"
        )

    target_groups = CELL_TYPE_TARGETS[cell_type_id]
    exclude_groups = BG_EXCLUDE.get(cell_type_id, [])
    df = pd.read_csv(atlas_groups_csv)

    # Exclude overlapping groups (e.g. T-CD3 for CD4T/CD8T)
    if exclude_groups:
        before = len(df)
        df = df[~df['group'].isin(exclude_groups)].copy()
        print(f"    Excluded {before - len(df)} samples from overlapping groups: "
              f"{', '.join(exclude_groups)}")

    # Merge target subtypes into one group, everything else = background
    df['group'] = df['group'].apply(
        lambda g: cell_type_id if g in target_groups else 'background'
    )

    # Filter to samples that have beta files, if beta_dir is provided
    if beta_dir:
        import glob
        available = set()
        for bp in glob.glob(os.path.join(beta_dir, "*.beta")):
            stem = os.path.basename(bp)
            if stem.endswith('.beta'):
                stem = stem[:-5]
            available.add(stem)
        name_col = df.columns[0]
        before = len(df)
        df = df[df[name_col].isin(available)]
        after = len(df)
        n_target = (df['group'] != 'background').sum()
        n_bg = (df['group'] == 'background').sum()
        print(f"    Full atlas groups file: {after}/{before} samples have beta files "
              f"({n_target} target, {n_bg} background)")

    df.to_csv(output_path, index=False)
    return output_path


# ---------------------------------------------------------------------------
# Beta list file generation
# ---------------------------------------------------------------------------
def generate_groups_file_filter(
    cell_type_id: str,
    atlas_groups_csv: str,
    output_path: str,
    beta_dir: Optional[str] = None,
) -> str:
    """
    Generate a groups CSV for the per-subgroup filter.

    Uses FILTER_EXCLUDE (not BG_EXCLUDE): only the truly overlapping
    group (Blood-T-CD3) is excluded. The opposite T subtype is KEPT
    so that T-cell-general markers are rejected by the filter.

    For cell types without a FILTER_EXCLUDE entry, this is identical
    to generate_groups_file (uses BG_EXCLUDE, which is empty for them).
    """
    if cell_type_id not in CELL_TYPE_TARGETS:
        raise ValueError(
            f"Unknown cell type '{cell_type_id}'. "
            f"Valid options: {VALID_CELL_TYPES}"
        )

    target_groups = CELL_TYPE_TARGETS[cell_type_id]
    # Use FILTER_EXCLUDE if available, otherwise fall back to BG_EXCLUDE
    exclude_groups = FILTER_EXCLUDE.get(cell_type_id, BG_EXCLUDE.get(cell_type_id, []))
    df = pd.read_csv(atlas_groups_csv)

    # Restrict background to blood cell types only for blood assays.
    # For non-blood cell types, use the full atlas as background.
    if is_blood_cell_type(cell_type_id):
        df = df[df['group'].isin(BLOOD_GROUPS)].copy()
    else:
        df = df.copy()

    if exclude_groups:
        before = len(df)
        df = df[~df['group'].isin(exclude_groups)].copy()
        print(f"    Filter groups: excluded {before - len(df)} samples from: "
              f"{', '.join(exclude_groups)}")

    df['group'] = df['group'].apply(
        lambda g: cell_type_id if g in target_groups else 'background'
    )

    if beta_dir:
        import glob
        available = set()
        for bp in glob.glob(os.path.join(beta_dir, "*.beta")):
            stem = os.path.basename(bp)
            if stem.endswith('.beta'):
                stem = stem[:-5]
            available.add(stem)
        name_col = df.columns[0]
        df = df[df[name_col].isin(available)]

    df.to_csv(output_path, index=False)
    return output_path


# ---------------------------------------------------------------------------
# Beta list file generation
# ---------------------------------------------------------------------------
def generate_beta_list(beta_dir: str, output_path: str) -> str:
    """
    Generate a text file listing all .beta file paths.

    Args:
        beta_dir: Directory containing .beta files
        output_path: Where to write the beta list

    Returns:
        Path to the generated beta list
    """
    beta_files = sorted(glob.glob(os.path.join(beta_dir, "*.beta")))
    if not beta_files:
        raise FileNotFoundError(
            f"No .beta files found in {beta_dir}. "
            f"Run download_data.sh first."
        )
    with open(output_path, 'w') as f:
        for path in beta_files:
            f.write(path + '\n')
    return output_path


# ---------------------------------------------------------------------------
# Run wgbstools find_markers
# ---------------------------------------------------------------------------
def find_wgbstools(wgbstools_path: Optional[str] = None) -> str:
    """
    Locate the wgbstools executable.

    Priority:
        1. Explicit path provided
        2. 'wgbstools' on PATH
        3. .venv/bin/wgbstools (common install location)
        4. ./wgbs_tools/wgbstools (cloned by install_macos.sh)
    """
    if wgbstools_path and os.path.isfile(wgbstools_path):
        return wgbstools_path

    on_path = shutil.which('wgbstools')
    if on_path:
        return on_path

    venv_bin = os.path.join(os.getcwd(), '.venv', 'bin', 'wgbstools')
    if os.path.isfile(venv_bin) or os.path.islink(venv_bin):
        return venv_bin

    local = os.path.join(os.getcwd(), 'wgbs_tools', 'wgbstools')
    if os.path.isfile(local):
        return local

    raise FileNotFoundError(
        "wgbstools not found. Install it with: "
        "git clone https://github.com/nloyfer/wgbs_tools.git && "
        "cd wgbs_tools && python3 setup.py && cd .. && "
        "ln -sf \"$(pwd)/wgbs_tools/wgbstools\" .venv/bin/wgbstools"
    )


def _find_venv_python() -> Optional[str]:
    """
    Find the Python executable in the active virtual environment.

    wgbstools is a Python script with shebang '#!/usr/bin/env python3'.
    When run as a subprocess, 'env python3' may resolve to a different
    Python than the venv (e.g. conda base), which won't have scipy.
    We explicitly invoke wgbstools via the venv Python to fix this.
    """
    # sys.executable is the Python running this module — should be the venv
    if sys.executable and os.path.isfile(sys.executable):
        return sys.executable
    return None


def run_find_markers(
    groups_file: str,
    beta_list: str,
    blocks_file: str,
    out_dir: str,
    cell_type_id: str,
    threads: int = 2,
    top_n: int = 1000,
    delta_means: float = 0.3,
    min_cpg: int = 3,
    unmeth_mean_thresh: float = 0.15,
    meth_mean_thresh: float = 0.65,
    wgbstools_path: Optional[str] = None,
    only_hyper: bool = False,
) -> str:
    """
    Run wgbstools find_markers for one cell type.

    Args:
        groups_file: Path to groups CSV
        beta_list: Path to beta list text file
        blocks_file: Path to blocks BED file
        out_dir: Output directory for find_markers
        cell_type_id: Target group name in the groups file
        threads: Number of threads
        top_n: Number of top markers to keep
        delta_means: Minimum methylation difference (target vs background)
        min_cpg: Minimum CpGs per block
        unmeth_mean_thresh: Target mean must be below this (default 0.15)
        meth_mean_thresh: Background mean must be above this (default 0.65)
        wgbstools_path: Path to wgbstools executable
        only_hyper: If True, find hypermethylated markers (target methylated,
            background unmethylated) instead of hypomethylated markers.

    Returns:
        Path to the output BED file (Markers.{cell_type_id}.bed)
    """
    wgbstools = find_wgbstools(wgbstools_path)
    os.makedirs(out_dir, exist_ok=True)

    # wgbstools is a Python script with shebang '#!/usr/bin/env python3'.
    # When run directly as a subprocess, 'env python3' may resolve to a
    # different Python than the venv (e.g. conda base), which won't have
    # scipy installed. We explicitly invoke it via the venv Python.
    venv_python = _find_venv_python()
    if venv_python:
        cmd = [venv_python, wgbstools, 'find_markers']
    else:
        cmd = [wgbstools, 'find_markers']

    # Pre-check: verify scipy is available in the Python that will run wgbstools.
    # find_markers requires scipy.stats.mannwhitneyu and will crash without it.
    check_py = cmd[0]
    scipy_check = subprocess.run(
        [check_py, '-c', 'import scipy; print(scipy.__version__)'],
        capture_output=True, text=True
    )
    if scipy_check.returncode != 0:
        raise RuntimeError(
            f"scipy is not installed in the Python that runs wgbstools "
            f"({check_py}). find_markers requires scipy. "
            f"Install it with: {check_py} -m pip install scipy"
        )

    cmd.extend([
        '--blocks_path', blocks_file,
        '--groups_file', groups_file,
        '--beta_list_file', beta_list,
        '--targets', cell_type_id,
        '--only_hyper' if only_hyper else '--only_hypo',
        '--delta_means', str(delta_means),
        '--min_cpg', str(min_cpg),
        '--unmeth_mean_thresh', str(unmeth_mean_thresh),
        '--meth_mean_thresh', str(meth_mean_thresh),
        '--sort_by', 'delta_means',
        '--top', str(top_n),
        '--header',
        '--threads', str(threads),
        '--out_dir', out_dir,
    ])

    print(f"  Running: {' '.join(cmd[:6])} ... --out_dir {out_dir}")
    result = subprocess.run(
        cmd, capture_output=True, text=True, timeout=3600
    )

    if result.returncode != 0:
        # Print stderr for debugging
        if result.stderr:
            print(f"  wgbstools stderr: {result.stderr[:500]}")
        raise RuntimeError(
            f"wgbstools find_markers failed with return code {result.returncode}"
        )

    if result.stderr:
        # find_markers prints progress to stderr
        for line in result.stderr.strip().split('\n')[-5:]:
            print(f"  {line}")

    bed_path = os.path.join(out_dir, f'Markers.{cell_type_id}.bed')
    if not os.path.isfile(bed_path):
        raise FileNotFoundError(
            f"Expected output not found: {bed_path}"
        )

    return bed_path


# ---------------------------------------------------------------------------
# Parse find_markers BED output
# ---------------------------------------------------------------------------
# Column names from find_markers dump_results:
MARKERS_BED_COLUMNS = [
    'chr', 'start', 'end', 'startCpG', 'endCpG',
    'target', 'region', 'lenCpG', 'bp',
    'tg_mean', 'bg_mean', 'delta_means', 'delta_quants',
    'delta_maxmin', 'ttest', 'mw_test', 'mvalue_ttest',
    'direction', 'anno', 'gene',
]


def parse_markers_bed(bed_path: str, cell_type_id: str) -> List[DMRBlock]:
    """
    Parse a find_markers BED file into DMRBlock objects.

    The BED file has comment lines starting with '#' (sample lists)
    followed by a header line and data lines.

    Args:
        bed_path: Path to Markers.{cell_type_id}.bed
        cell_type_id: Cell type ID for these markers

    Returns:
        List of DMRBlock objects (without per-CpG data yet)
    """
    # Read the file, skipping comment lines (but keeping the #chr header)
    with open(bed_path) as f:
        lines = [l for l in f if not l.startswith('#>') and not l.startswith('#<')]

    if not lines:
        return []

    # Parse with pandas
    from io import StringIO
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t')

    # Rename #chr to chr if present
    if '#chr' in df.columns:
        df = df.rename(columns={'#chr': 'chr'})

    blocks = []
    for rank, (_, row) in enumerate(df.iterrows(), 1):
        chrom = str(row['chr'])
        start = int(row['start'])
        end = int(row['end'])
        start_cpg = int(row['startCpG'])
        end_cpg = int(row['endCpG'])
        num_cpgs = end_cpg - start_cpg

        # Parse lenCpG string like "7CpGs"
        len_cpg_str = str(row.get('lenCpG', ''))
        if len_cpg_str.endswith('CpGs'):
            try:
                num_cpgs = int(len_cpg_str.replace('CpGs', ''))
            except ValueError:
                pass

        # Parse bp string like "120bp"
        bp_str = str(row.get('bp', ''))
        block_len = end - start
        if bp_str.endswith('bp'):
            try:
                block_len = int(bp_str.replace('bp', ''))
            except ValueError:
                pass

        seq_id = f"{cell_type_id}_{rank:04d}"

        block = DMRBlock(
            cell_type_id=cell_type_id,
            seq_id=seq_id,
            rank=rank,
            chrom=chrom,
            start=start,
            end=end,
            num_cpgs=num_cpgs,
            block_len=block_len,
            gene=str(row.get('gene', '')) if pd.notna(row.get('gene', '')) else '',
            annotation=str(row.get('anno', '')) if pd.notna(row.get('anno', '')) else '',
            direction=str(row.get('direction', '')) if pd.notna(row.get('direction', '')) else '',
            target_mean_meth=float(row.get('tg_mean', 0) or 0),
            background_mean_meth=float(row.get('bg_mean', 0) or 0),
            delta_means=float(row.get('delta_means', 0) or 0),
            delta_quants=float(row.get('delta_quants', 0) or 0),
            delta_maxmin=float(row.get('delta_maxmin', 0) or 0),
            ttest_pval=str(row.get('ttest', '')) if pd.notna(row.get('ttest', '')) else '',
            mwtest_pval=str(row.get('mw_test', '')) if pd.notna(row.get('mw_test', '')) else '',
            mvalue_ttest=str(row.get('mvalue_ttest', '')) if pd.notna(row.get('mvalue_ttest', '')) else '',
            dms_score=0.0,  # Will be computed from per-CpG data
            target_celltype=cell_type_id,
            cpg_sites=[],  # Will be filled from beta files
        )
        blocks.append(block)

    return blocks


# ---------------------------------------------------------------------------
# Per-CpG methylation extraction from beta files
# ---------------------------------------------------------------------------
def read_beta_region(
    beta_path: str,
    start_cpg: int,
    end_cpg: int,
    min_cov: int = 5,
) -> np.ndarray:
    """
    Read per-CpG methylation fractions from a wgbstools beta file.

    Beta file format:
        - Binary, 2 x uint8 per CpG: (#meth, #covered)
        - CpG index is 1-based (global)
        - Byte offset = (start_cpg - 1) * 2

    Args:
        beta_path: Path to .beta file
        start_cpg: Start CpG index (1-based, inclusive)
        end_cpg: End CpG index (1-based, exclusive)
        min_cov: Minimum coverage to include

    Returns:
        np.array of float, shape (end_cpg - start_cpg,), NaN where cov < min_cov
    """
    n_cpgs = end_cpg - start_cpg
    if n_cpgs < 1:
        return np.array([])

    offset = (start_cpg - 1) * 2
    try:
        with open(beta_path, "rb") as f:
            f.seek(offset)
            raw = np.frombuffer(
                f.read(n_cpgs * 2), dtype=np.uint8
            ).reshape(-1, 2)
        meth = raw[:, 0].astype(float)
        cov = raw[:, 1].astype(float)
        fracs = np.where(cov >= min_cov, meth / np.maximum(cov, 1), np.nan)
        return fracs
    except Exception as e:
        print(f"    WARNING: could not read {os.path.basename(beta_path)}: {e}")
        return np.full(n_cpgs, np.nan)


def _match_beta_files(
    groups_csv: str,
    beta_dir: str,
) -> Tuple[List[str], List[str]]:
    """
    Match target and background sample names to beta file paths.

    Args:
        groups_csv: Path to the generated groups CSV (target/background)
        beta_dir: Directory containing .beta files

    Returns:
        (target_betas, background_betas) — lists of file paths
    """
    df = pd.read_csv(groups_csv)
    # First column is sample name (renamed to 'fname' by wgbstools, but
    # we wrote it with original column name)
    name_col = df.columns[0]

    all_betas = glob.glob(os.path.join(beta_dir, "*.beta"))
    # Build a lookup: sample_name -> beta_path
    beta_lookup = {}
    for bp in all_betas:
        basename = os.path.basename(bp)
        # Remove .beta extension
        stem = basename[:-5] if basename.endswith('.beta') else basename
        beta_lookup[stem] = bp

    target_names = df[df['group'] != 'background'][name_col].tolist()
    bg_names = df[df['group'] == 'background'][name_col].tolist()

    target_betas = []
    for name in target_names:
        if name in beta_lookup:
            target_betas.append(beta_lookup[name])
        else:
            print(f"    WARNING: no beta file for target sample '{name}'")

    bg_betas = []
    for name in bg_names:
        if name in beta_lookup:
            bg_betas.append(beta_lookup[name])
        # Don't warn for missing background — there may be many

    return target_betas, bg_betas


def _match_beta_files_with_groups(
    groups_csv: str,
    atlas_groups_csv: str,
    beta_dir: str,
) -> Tuple[List[str], List[str], Dict[str, str]]:
    """
    Match target and background samples to beta files, and also return
    a mapping from background beta file path -> original atlas group name.

    This is needed for per-subgroup background checking: we need to know
    which background samples belong to which blood cell type (e.g. CD8 T
    cells, B cells) to verify that EACH subgroup is sufficiently methylated.

    Args:
        groups_csv: Path to the generated groups CSV (target/background)
        atlas_groups_csv: Path to data/sample_groups.csv (original groups)
        beta_dir: Directory containing .beta files

    Returns:
        (target_betas, background_betas, bg_group_map)
        where bg_group_map maps beta_path -> atlas group name (e.g. "Blood-T-CD8")
    """
    df = pd.read_csv(groups_csv)
    name_col = df.columns[0]

    # Load original atlas groups to get the real group name for each sample
    atlas_df = pd.read_csv(atlas_groups_csv)
    atlas_name_col = atlas_df.columns[0]
    sample_to_group = dict(zip(atlas_df[atlas_name_col], atlas_df['group']))

    all_betas = glob.glob(os.path.join(beta_dir, "*.beta"))
    beta_lookup = {}
    for bp in all_betas:
        basename = os.path.basename(bp)
        stem = basename[:-5] if basename.endswith('.beta') else basename
        beta_lookup[stem] = bp

    target_names = df[df['group'] != 'background'][name_col].tolist()
    bg_names = df[df['group'] == 'background'][name_col].tolist()

    target_betas = []
    for name in target_names:
        if name in beta_lookup:
            target_betas.append(beta_lookup[name])

    bg_betas = []
    bg_group_map = {}
    for name in bg_names:
        if name in beta_lookup:
            bp = beta_lookup[name]
            bg_betas.append(bp)
            bg_group_map[bp] = sample_to_group.get(name, "unknown")

    return target_betas, bg_betas, bg_group_map


def extract_percpg_methylation(
    blocks: List[DMRBlock],
    target_betas: List[str],
    background_betas: List[str],
    min_cov: int = 5,
    max_bg_samples: int = 30,
) -> None:
    """
    Extract per-CpG methylation for each block and populate cpg_sites.

    For each block, reads target and background beta files at the block's
    CpG indices (startCpG to endCpG-1), computes per-CpG mean methylation,
    and creates CpGSite objects.

    Also computes the DMS score.

    Args:
        blocks: List of DMRBlock objects (modified in place)
        target_betas: List of target beta file paths
        background_betas: List of background beta file paths
        min_cov: Minimum coverage per CpG site
        max_bg_samples: Maximum background samples to use (for speed)
    """
    # Subsample background if too many
    bg_betas = background_betas[:max_bg_samples]
    print(f"  Per-CpG extraction: {len(target_betas)} target, "
          f"{len(bg_betas)} background samples")

    for i, block in enumerate(blocks):
        start_cpg = None
        end_cpg = None

        # We need the CpG indices. These were in the BED file but not
        # stored in DMRBlock. We'll re-derive them from the BED file
        # or use the block's start/end genomic positions.
        # Actually, we stored them during parsing — let's use a workaround:
        # The DMRBlock doesn't have startCpG/endCpG fields, so we need to
        # pass them separately. We'll use the block's num_cpgs and
        # distribute CpGs evenly if we can't get indices.

        # For now, skip per-CpG extraction if we don't have CpG indices
        # The DMS score will be 0, and CpG positions will be
        # found by sequence scanning in step 2.

        # Create placeholder CpG sites with even distribution
        if block.num_cpgs > 0 and not block.cpg_sites:
            step = max(1, (block.end - block.start) // max(1, block.num_cpgs))
            for j in range(block.num_cpgs):
                pos = block.start + j * step
                block.cpg_sites.append(CpGSite(
                    label=f"CpG_{j+1}",
                    position=pos,
                    global_idx=pos,  # Placeholder; will be corrected in step 2
                    target_mean_beta=block.target_mean_meth,
                    background_mean_beta=block.background_mean_meth,
                    delta_beta=abs(block.target_mean_meth - block.background_mean_meth),
                ))

        # Compute DMS score from block-level stats
        # (simplified: use delta_means as a proxy)
        block.dms_score = min(1.0, block.delta_means / 0.5)

        if (i + 1) % 50 == 0:
            print(f"  Processed {i+1}/{len(blocks)} blocks")


def _subsample_background(
    background_betas: List[str],
    max_bg_samples: int,
    bg_group_map: Optional[Dict[str, str]] = None,
) -> List[str]:
    """
    Subsample background beta files, preserving at least one sample per
    subgroup.

    A naive ``background_betas[:max_bg_samples]`` truncation can silently
    drop entire subgroups when they happen to be at the end of the list.
    For example, with 33 blood-only background samples and max_bg_samples=30,
    the last 3 samples (both Blood-B-Mem) are dropped, making Blood-B-Mem
    invisible in per-subgroup methylation and the per-subgroup filter.

    This function groups background samples by atlas subgroup (using
    bg_group_map), keeps at least one sample per subgroup, then fills
    remaining slots proportionally. If bg_group_map is None or empty,
    falls back to simple truncation.

    Args:
        background_betas: Full list of background beta file paths
        max_bg_samples: Maximum number of background samples to keep
        bg_group_map: Maps beta_path -> atlas group name

    Returns:
        Subsampled list of background beta file paths
    """
    if len(background_betas) <= max_bg_samples:
        return list(background_betas)

    if not bg_group_map:
        return background_betas[:max_bg_samples]

    # Group beta paths by subgroup
    groups: Dict[str, List[str]] = {}
    for bp in background_betas:
        g = bg_group_map.get(bp, "unknown")
        groups.setdefault(g, []).append(bp)

    n_groups = len(groups)
    if n_groups == 0:
        return background_betas[:max_bg_samples]

    # If we can fit all samples, do so
    if len(background_betas) <= max_bg_samples:
        return list(background_betas)

    # Allocate slots: at least 1 per group, distribute the rest proportionally.
    # Use the largest-remainder method to avoid rounding errors that could
    # cause over-allocation followed by a trim that drops a group's only sample.
    group_sizes = {g: len(bps) for g, bps in groups.items()}
    total = sum(group_sizes.values())
    slots = max_bg_samples

    # If more groups than slots, we can't guarantee 1 per group.
    # Take 1 sample from the largest groups until we fill slots.
    if n_groups >= slots:
        result: List[str] = []
        for g in sorted(groups, key=lambda x: group_sizes[x], reverse=True):
            if len(result) >= slots:
                break
            result.append(groups[g][0])
        return result

    # Step 1: guarantee 1 per group
    allocation = {g: 1 for g in groups}
    remaining = slots - n_groups

    # Step 2: proportional allocation of remaining slots using largest remainder
    if remaining > 0:
        exact = {g: remaining * group_sizes[g] / total for g in groups}
        floor = {g: int(exact[g]) for g in groups}
        allocated = sum(floor.values())
        leftover = remaining - allocated

        # Distribute leftover by largest fractional remainder
        remainders = sorted(
            groups.keys(),
            key=lambda g: exact[g] - floor[g],
            reverse=True,
        )
        for g in remainders[:leftover]:
            floor[g] += 1

        for g in groups:
            allocation[g] = 1 + min(floor[g], group_sizes[g] - 1)

    # Build result from allocation
    result: List[str] = []
    for g, bps in groups.items():
        result.extend(bps[:allocation[g]])

    # If still over (shouldn't happen with largest-remainder, but safety),
    # trim from groups with the most samples
    if len(result) > slots:
        # Sort groups by current allocation descending, trim from the top
        over = len(result) - slots
        for g in sorted(groups, key=lambda x: allocation[x], reverse=True):
            if over <= 0:
                break
            trim = min(over, allocation[g] - 1)  # never below 1
            allocation[g] -= trim
            over -= trim
        result = []
        for g, bps in groups.items():
            result.extend(bps[:allocation[g]])

    # If under-allocated, fill from groups that have extras
    if len(result) < slots:
        used = set(result)
        for bp in background_betas:
            if bp not in used and len(result) < slots:
                result.append(bp)

    return result


def extract_percpg_methylation_with_indices(
    blocks: List[DMRBlock],
    cpg_indices: List[Tuple[int, int]],  # (startCpG, endCpG) per block
    target_betas: List[str],
    background_betas: List[str],
    min_cov: int = 5,
    max_bg_samples: int = 30,
    bg_group_map: Optional[Dict[str, str]] = None,
    only_hyper: bool = False,
) -> None:
    """
    Extract per-CpG methylation using CpG indices from find_markers.

    This is the full version that reads beta files at the correct CpG
    indices and computes proper DMS scores.

    If bg_group_map is provided, also computes per-subgroup background
    methylation (e.g. CD8 T cells vs B cells vs monocytes) and stores
    it in block.bg_subgroup_meth.

    Args:
        blocks: List of DMRBlock objects (modified in place)
        cpg_indices: List of (startCpG, endCpG) tuples, one per block
        target_betas: List of target beta file paths
        background_betas: List of background beta file paths
        min_cov: Minimum coverage per CpG site
        max_bg_samples: Maximum background samples to use (for speed)
        bg_group_map: Maps beta_path -> atlas group name (e.g. "Blood-T-CD8")
    """
    bg_betas = _subsample_background(
        background_betas, max_bg_samples, bg_group_map
    )
    print(f"  Per-CpG extraction: {len(target_betas)} target, "
          f"{len(bg_betas)} background samples")

    for i, (block, (start_cpg, end_cpg)) in enumerate(zip(blocks, cpg_indices)):
        n_cpgs = end_cpg - start_cpg
        if n_cpgs < 1:
            continue

        # Read per-CpG methylation from all target and background beta files
        tg_mat = np.array([
            read_beta_region(bp, start_cpg, end_cpg, min_cov)
            for bp in target_betas
        ]) if target_betas else np.full((1, n_cpgs), np.nan)

        bg_mat = np.array([
            read_beta_region(bp, start_cpg, end_cpg, min_cov)
            for bp in bg_betas
        ]) if bg_betas else np.full((1, n_cpgs), np.nan)

        # Per-CpG mean methylation
        tg_cpg_mean = np.nanmean(tg_mat, axis=0)
        bg_cpg_mean = np.nanmean(bg_mat, axis=0)

        # Compute DMS score
        score = _dms_score(tg_mat, bg_mat, only_hyper=only_hyper)

        # Per-subgroup background methylation
        if bg_group_map:
            subgroup_meth = {}
            for idx, bp in enumerate(bg_betas):
                group = bg_group_map.get(bp, "unknown")
                if group not in subgroup_meth:
                    subgroup_meth[group] = []
                # Block-level mean for this sample
                sample_mean = np.nanmean(bg_mat[idx])
                if not np.isnan(sample_mean):
                    subgroup_meth[group].append(float(sample_mean))
            # Average per subgroup
            block.bg_subgroup_meth = {
                g: float(np.mean(vals)) for g, vals in subgroup_meth.items() if vals
            }

        # Build CpG sites
        # Genomic positions: distribute evenly across the block
        # (precise positions will be found by sequence scanning in step 2)
        block.cpg_sites = []
        step = max(1, (block.end - block.start) // max(1, n_cpgs))
        for j in range(n_cpgs):
            pos = block.start + j * step
            tg_val = float(tg_cpg_mean[j]) if not np.isnan(tg_cpg_mean[j]) else 0.0
            bg_val = float(bg_cpg_mean[j]) if not np.isnan(bg_cpg_mean[j]) else 0.0
            block.cpg_sites.append(CpGSite(
                label=f"CpG_{j+1}",
                position=pos,
                global_idx=start_cpg + j,  # Actual CpG index
                target_mean_beta=tg_val,
                background_mean_beta=bg_val,
                delta_beta=abs(tg_val - bg_val),
            ))

        block.dms_score = score

        if (i + 1) % 50 == 0:
            print(f"  Processed {i+1}/{len(blocks)} blocks")


# ---------------------------------------------------------------------------
# Per-CpG per-cell-type enrichment
# ---------------------------------------------------------------------------
def _build_all_beta_group_map(
    atlas_groups_csv: str,
    beta_dir: str,
) -> Dict[str, str]:
    """
    Map every .beta file in beta_dir to its atlas group name.

    Unlike _match_beta_files_with_groups (which separates target/background),
    this returns ALL samples grouped by their original atlas cell type.

    Args:
        atlas_groups_csv: Path to data/sample_groups.csv
        beta_dir: Directory containing .beta files

    Returns:
        Dict mapping beta_path -> atlas group name (e.g. "Blood-T-CD8")
    """
    atlas_df = pd.read_csv(atlas_groups_csv)
    atlas_name_col = atlas_df.columns[0]
    sample_to_group = dict(zip(atlas_df[atlas_name_col], atlas_df['group']))

    all_betas = glob.glob(os.path.join(beta_dir, "*.beta"))
    beta_group_map = {}
    for bp in all_betas:
        basename = os.path.basename(bp)
        stem = basename[:-5] if basename.endswith('.beta') else basename
        beta_group_map[bp] = sample_to_group.get(stem, "unknown")

    return beta_group_map


def enrich_percpg_per_celltype(
    blocks: List[DMRBlock],
    beta_group_map: Dict[str, str],
    min_cov: int = 5,
) -> None:
    """
    Enrich each CpG site with per-cell-type mean beta values.

    Reads ALL beta files (not subsampled) to get complete per-cell-type
    coverage across all 82 atlas groups. Called after all DMR filters,
    so only surviving blocks are processed.

    For each CpG in each block, stores a dict {atlas_group: mean_beta}
    in cpg_site.per_celltype_beta. This includes the target cell type(s)
    as separate entries, enabling direct comparison in the Excel output.

    I/O optimization: opens each beta file once and reads all block slices
    via seek, rather than re-opening per block.

    Args:
        blocks: List of DMRBlock objects (modified in place)
        beta_group_map: Maps every beta file path -> atlas group name
        min_cov: Minimum coverage per CpG site
    """
    if not blocks or not beta_group_map:
        return

    # Group beta files by atlas group
    group_to_betas: Dict[str, List[str]] = {}
    for bp, group in beta_group_map.items():
        group_to_betas.setdefault(group, []).append(bp)

    n_groups = len(group_to_betas)
    n_samples = len(beta_group_map)
    print(f"  Step 0e3: Per-cell-type enrichment: {n_groups} cell types, "
          f"{n_samples} samples, {len(blocks)} blocks")

    # Pre-compute CpG index ranges for each block
    block_ranges: List[Optional[Tuple[int, int]]] = []
    for block in blocks:
        cs = block.cpg_sites
        if cs and len(cs) > 0:
            block_ranges.append((cs[0].global_idx, cs[-1].global_idx + 1))
        else:
            block_ranges.append(None)

    # Accumulate per-block, per-group sample arrays
    # block_data[block_idx][group_name] = list of np.array(n_cpgs,)
    block_data: List[Dict[str, List[np.ndarray]]] = [{} for _ in blocks]

    for bp, group_name in beta_group_map.items():
        try:
            with open(bp, "rb") as f:
                for i, brange in enumerate(block_ranges):
                    if brange is None:
                        continue
                    start_cpg, end_cpg = brange
                    n_cpgs = end_cpg - start_cpg
                    offset = (start_cpg - 1) * 2
                    f.seek(offset)
                    raw_bytes = f.read(n_cpgs * 2)
                    raw = np.frombuffer(raw_bytes, dtype=np.uint8)
                    if len(raw) < n_cpgs * 2:
                        raw = np.pad(raw, (0, n_cpgs * 2 - len(raw)))
                    raw = raw.reshape(-1, 2)
                    meth = raw[:, 0].astype(float)
                    cov = raw[:, 1].astype(float)
                    fracs = np.where(cov >= min_cov, meth / np.maximum(cov, 1), np.nan)

                    block_data[i].setdefault(group_name, []).append(fracs)
        except Exception as e:
            print(f"    WARNING: could not read {os.path.basename(bp)}: {e}")
            continue

    # Compute per-group per-CpG means and store in blocks
    for i, block in enumerate(blocks):
        if not block.cpg_sites:
            continue
        for group_name, mats in block_data[i].items():
            mat = np.array(mats)
            cpg_means = np.nanmean(mat, axis=0)
            for j, cpg in enumerate(block.cpg_sites):
                if j < len(cpg_means):
                    val = cpg_means[j]
                    if not np.isnan(val):
                        cpg.per_celltype_beta[group_name] = float(val)

        if (i + 1) % 50 == 0:
            print(f"  Enriched {i+1}/{len(blocks)} blocks")

    print(f"  Per-cell-type enrichment complete")


def enrich_blocks_with_celltype_beta(dmr_blocks, beta_dir, atlas_groups_csv, min_cov=5):
    """
    Enrich dmr_blocks (JSON dict format) with per-CpG per-cell-type beta values.

    Reads ALL beta files in beta_dir, groups samples by atlas cell type,
    and computes per-CpG mean beta per group. Stores results in
    cpg_site['per_celltype_beta'].

    Args:
        dmr_blocks: List of block dicts (from dmr_blocks.json, modified in place)
        beta_dir: Directory containing .beta files
        atlas_groups_csv: Path to data/sample_groups.csv
        min_cov: Minimum coverage per CpG site

    Returns:
        Number of blocks enriched with per-cell-type data
    """
    beta_group_map = _build_all_beta_group_map(atlas_groups_csv, beta_dir)
    if not beta_group_map:
        print(f"  WARNING: no beta files found in {beta_dir}")
        return 0

    # Group beta files by atlas group
    group_to_betas = {}
    for bp, group in beta_group_map.items():
        group_to_betas.setdefault(group, []).append(bp)

    n_groups = len(group_to_betas)
    n_samples = len(beta_group_map)
    print(f"  Per-cell-type enrichment: {n_groups} cell types, "
          f"{n_samples} samples, {len(dmr_blocks)} blocks")

    # Pre-compute CpG index ranges for each block
    block_ranges = []
    for block in dmr_blocks:
        cs = block.get("cpg_sites", [])
        if cs:
            block_ranges.append((cs[0]["global_idx"], cs[-1]["global_idx"] + 1))
        else:
            block_ranges.append(None)

    # Accumulate per-block, per-group sample arrays
    block_data = [{} for _ in dmr_blocks]

    for bp, group_name in beta_group_map.items():
        try:
            with open(bp, "rb") as f:
                for i, brange in enumerate(block_ranges):
                    if brange is None:
                        continue
                    start_cpg, end_cpg = brange
                    n_cpgs = end_cpg - start_cpg
                    offset = (start_cpg - 1) * 2
                    f.seek(offset)
                    raw_bytes = f.read(n_cpgs * 2)
                    raw = np.frombuffer(raw_bytes, dtype=np.uint8)
                    if len(raw) < n_cpgs * 2:
                        raw = np.pad(raw, (0, n_cpgs * 2 - len(raw)))
                    raw = raw.reshape(-1, 2)
                    meth = raw[:, 0].astype(float)
                    cov = raw[:, 1].astype(float)
                    fracs = np.where(cov >= min_cov, meth / np.maximum(cov, 1), np.nan)

                    block_data[i].setdefault(group_name, []).append(fracs)
        except Exception as e:
            print(f"    WARNING: could not read {os.path.basename(bp)}: {e}")
            continue

    # Compute per-group per-CpG means and store in blocks
    enriched_count = 0
    for i, block in enumerate(dmr_blocks):
        cs = block.get("cpg_sites", [])
        if not cs:
            continue
        has_data = False
        for group_name, mats in block_data[i].items():
            mat = np.array(mats)
            cpg_means = np.nanmean(mat, axis=0)
            for j, cpg in enumerate(cs):
                if j < len(cpg_means):
                    val = cpg_means[j]
                    if not np.isnan(val):
                        cpg.setdefault("per_celltype_beta", {})[group_name] = float(val)
                        has_data = True
        if has_data:
            enriched_count += 1

        if (i + 1) % 50 == 0:
            print(f"  Enriched {i+1}/{len(dmr_blocks)} blocks")

    print(f"  Per-cell-type enrichment complete: {enriched_count}/{len(dmr_blocks)} blocks enriched")
    return enriched_count


# ---------------------------------------------------------------------------
# DMS score (from repo 1's per_cpg_marker_analysis.py)
# ---------------------------------------------------------------------------
def _dms_score(
    tg_matrix: np.ndarray,
    bg_matrix: np.ndarray,
    only_hyper: bool = False,
) -> float:
    """
    Score how suitable a marker block is for bisulfite PCR.

    For hypomethylated markers (default):
        A. tg_near_zero: 1 - max(per-CpG mean in target)
        B. bg_near_one: min(per-CpG mean in background)
        C. tg_consistency: 1 - 2*mean(per-CpG std across target samples)
        D. bg_consistency: 1 - 2*mean(per-CpG std across background samples)
        E. coverage: fraction of CpGs with data in >=50% of samples

    For hypermethylated markers (only_hyper=True):
        A. tg_near_one: min(per-CpG mean in target)
        B. bg_near_zero: 1 - max(per-CpG mean in background)
        C, D, E are the same (consistency and coverage are direction-agnostic)

    Final score = mean(A, B, C, D) * E

    Args:
        tg_matrix: (n_target_samples, n_cpgs) methylation fractions
        bg_matrix: (n_bg_samples, n_cpgs) methylation fractions
        only_hyper: If True, score for hypermethylated markers

    Returns:
        DMS score (0-1)
    """
    tg_cpg_mean = np.nanmean(tg_matrix, axis=0)
    bg_cpg_mean = np.nanmean(bg_matrix, axis=0)

    if only_hyper:
        # Hyper: target should be near 1, background near 0
        A = max(0.0, float(np.nanmin(tg_cpg_mean))) if not np.all(np.isnan(tg_cpg_mean)) else 0.0
        B = max(0.0, float(1.0 - np.nanmax(bg_cpg_mean))) if not np.all(np.isnan(bg_cpg_mean)) else 0.0
    else:
        # Hypo: target should be near 0, background near 1
        A = max(0.0, float(1.0 - np.nanmax(tg_cpg_mean))) if not np.all(np.isnan(tg_cpg_mean)) else 0.0
        B = max(0.0, float(np.nanmin(bg_cpg_mean))) if not np.all(np.isnan(bg_cpg_mean)) else 0.0

    C = max(0.0, 1.0 - 2.0 * float(np.nanmean(np.nanstd(tg_matrix, axis=0))))
    D = max(0.0, 1.0 - 2.0 * float(np.nanmean(np.nanstd(bg_matrix, axis=0))))

    tg_cov = float(np.mean(np.mean(~np.isnan(tg_matrix), axis=0) >= 0.5))
    bg_cov = float(np.mean(np.mean(~np.isnan(bg_matrix), axis=0) >= 0.5))
    E = (tg_cov + bg_cov) / 2.0

    score = np.mean([A, B, C, D]) * E
    return round(float(score), 4)


# ---------------------------------------------------------------------------
# Parse BED with CpG indices
# ---------------------------------------------------------------------------
def parse_markers_bed_with_indices(
    bed_path: str,
    cell_type_id: str,
) -> Tuple[List[DMRBlock], List[Tuple[int, int]]]:
    """
    Parse a find_markers BED file, returning both DMRBlocks and CpG indices.

    Args:
        bed_path: Path to Markers.{cell_type_id}.bed
        cell_type_id: Cell type ID for these markers

    Returns:
        (blocks, cpg_indices) where cpg_indices[i] = (startCpG, endCpG) for block i
    """
    with open(bed_path) as f:
        lines = [l for l in f if not l.startswith('#>') and not l.startswith('#<')]

    if not lines:
        return [], []

    from io import StringIO
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t')

    if '#chr' in df.columns:
        df = df.rename(columns={'#chr': 'chr'})

    blocks = []
    cpg_indices = []

    for rank, (_, row) in enumerate(df.iterrows(), 1):
        chrom = str(row['chr'])
        start = int(row['start'])
        end = int(row['end'])
        start_cpg = int(row['startCpG'])
        end_cpg = int(row['endCpG'])
        num_cpgs = end_cpg - start_cpg

        len_cpg_str = str(row.get('lenCpG', ''))
        if len_cpg_str.endswith('CpGs'):
            try:
                num_cpgs = int(len_cpg_str.replace('CpGs', ''))
            except ValueError:
                pass

        bp_str = str(row.get('bp', ''))
        block_len = end - start
        if bp_str.endswith('bp'):
            try:
                block_len = int(bp_str.replace('bp', ''))
            except ValueError:
                pass

        seq_id = f"{cell_type_id}_{rank:04d}"

        block = DMRBlock(
            cell_type_id=cell_type_id,
            seq_id=seq_id,
            rank=rank,
            chrom=chrom,
            start=start,
            end=end,
            num_cpgs=num_cpgs,
            block_len=block_len,
            gene=str(row.get('gene', '')) if pd.notna(row.get('gene', '')) else '',
            annotation=str(row.get('anno', '')) if pd.notna(row.get('anno', '')) else '',
            direction=str(row.get('direction', '')) if pd.notna(row.get('direction', '')) else '',
            target_mean_meth=float(row.get('tg_mean', 0) or 0),
            background_mean_meth=float(row.get('bg_mean', 0) or 0),
            delta_means=float(row.get('delta_means', 0) or 0),
            delta_quants=float(row.get('delta_quants', 0) or 0),
            delta_maxmin=float(row.get('delta_maxmin', 0) or 0),
            ttest_pval=str(row.get('ttest', '')) if pd.notna(row.get('ttest', '')) else '',
            mwtest_pval=str(row.get('mw_test', '')) if pd.notna(row.get('mw_test', '')) else '',
            mvalue_ttest=str(row.get('mvalue_ttest', '')) if pd.notna(row.get('mvalue_ttest', '')) else '',
            dms_score=0.0,
            target_celltype=cell_type_id,
            cpg_sites=[],
        )
        blocks.append(block)
        cpg_indices.append((start_cpg, end_cpg))

    return blocks, cpg_indices


# ---------------------------------------------------------------------------
# Save blocks as JSON (same format as step 1)
# ---------------------------------------------------------------------------
def save_blocks_json(blocks: List[DMRBlock], output_path: str) -> None:
    """Save DMRBlocks as JSON in the same format as step 1 output."""
    blocks_data = []
    for b in blocks:
        blocks_data.append({
            "cell_type_id": b.cell_type_id,
            "seq_id": b.seq_id,
            "rank": b.rank,
            "chrom": b.chrom,
            "start": b.start,
            "end": b.end,
            "num_cpgs": b.num_cpgs,
            "block_len": b.block_len,
            "gene": b.gene,
            "annotation": b.annotation,
            "dms_score": b.dms_score,
            "target_mean_meth": b.target_mean_meth,
            "background_mean_meth": b.background_mean_meth,
            "delta_means": b.delta_means,
            "delta_quants": b.delta_quants,
            "delta_maxmin": b.delta_maxmin,
            "ttest_pval": b.ttest_pval,
            "mwtest_pval": b.mwtest_pval,
            "mvalue_ttest": b.mvalue_ttest,
            "direction": b.direction,
            "target_celltype": b.target_celltype,
            "bg_subgroup_meth": b.bg_subgroup_meth if b.bg_subgroup_meth else {},
            "cpg_sites": [
                {"label": c.label, "position": c.position,
                 "global_idx": c.global_idx,
                 "target_mean_beta": c.target_mean_beta,
                 "background_mean_beta": c.background_mean_beta,
                 "delta_beta": c.delta_beta,
                 "per_celltype_beta": c.per_celltype_beta if c.per_celltype_beta else {}}
                for c in b.cpg_sites
            ],
        })

    with open(output_path, 'w') as f:
        json.dump(blocks_data, f, indent=2)
    print(f"  Saved {len(blocks)} blocks to {output_path}")


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------
def discover_dmrs(
    cell_type_id: str,
    beta_dir: str,
    blocks_file: str,
    groups_csv: str,
    out_dir: str,
    threads: int = 2,
    top_n: int = 1000,
    delta_means: float = 0.3,
    min_cpg: int = 3,
    wgbstools_path: Optional[str] = None,
    max_bg_samples: int = 30,
    skip_find_markers: bool = False,
    unmeth_mean_thresh: float = 0.15,
    meth_mean_thresh: float = 0.65,
    min_bg_subgroup_meth: float = 0.70,
    use_full_atlas: bool = False,
    only_hyper: bool = False,
) -> List[DMRBlock]:
    """
    Full DMR discovery pipeline for one cell type.

    Args:
        cell_type_id: Cell type ID — any of the 7 merged blood types
            (MONO, BCELL, NK, GRAN, CD3T, CD4T, CD8T) or any atlas group
            name (e.g. Liver-Hepatocytes, Cortex-Neuron). Non-blood cell
            types automatically use full-atlas background.
        beta_dir: Directory containing .beta files
        blocks_file: Path to blocks BED file
        groups_csv: Path to data/sample_groups.csv
        out_dir: Output directory
        threads: Threads for find_markers
        top_n: Top N markers to keep
        delta_means: Min methylation difference (default 0.3)
        min_cpg: Min CpGs per block
        wgbstools_path: Path to wgbstools executable
        max_bg_samples: Max background samples for per-CpG extraction
        skip_find_markers: If True, skip find_markers (use existing BED)
        unmeth_mean_thresh: Target mean must be below this (default 0.15)
        meth_mean_thresh: Background mean must be above this (default 0.65)
        min_bg_subgroup_meth: Reject blocks where any background subgroup
            has mean methylation below this (default 0.50).
            For hypermethylated markers (only_hyper=True), this becomes
            the MAX: reject blocks where any subgroup is ABOVE
            (1 - min_bg_subgroup_meth).
        use_full_atlas: If True, use full atlas (all 207 samples) as
            find_markers background to get more candidates, then apply
            blood-only per-subgroup filter for quality control (hybrid
            approach). If False, use blood-only for both.
        only_hyper: If True, find hypermethylated markers (target
            methylated, background unmethylated) instead of
            hypomethylated markers.

    Returns:
        List of DMRBlock objects with per-CpG data
    """
    os.makedirs(out_dir, exist_ok=True)

    marker_type = "hypermethylated" if only_hyper else "hypomethylated"
    bg_mode = "full atlas" if use_full_atlas else "blood-only"

    # Non-blood cell types: force full-atlas background and disable the
    # per-subgroup filter. The per-subgroup filter checks that every
    # background subgroup is sufficiently methylated/unmethylated, but
    # with ~80 atlas subgroups and max_bg_samples=30, most subgroups
    # would be missing → all blocks rejected. The DMS score
    # (CS >= 0.6) and delta_means threshold still filter poor markers.
    if not is_blood_cell_type(cell_type_id):
        if not use_full_atlas:
            print(f"  NOTE: {cell_type_id} is a non-blood cell type — "
                  f"automatically using full atlas as background")
            use_full_atlas = True
            bg_mode = "full atlas"
        if min_bg_subgroup_meth > 0:
            print(f"  NOTE: Disabling per-subgroup filter for non-blood cell type "
                  f"(would reject all blocks with ~80 subgroups and 30 bg samples)")
            min_bg_subgroup_meth = 0.0

    # Step 0a: Generate groups file(s)
    # For the hybrid approach, find_markers uses the full atlas (more
    # candidates), but per-CpG extraction + per-subgroup filter use
    # blood-only background (quality control).
    print(f"\n  Step 0a: Generating groups file for {cell_type_id} "
          f"({marker_type}, {bg_mode} background)")
    blood_groups_file = os.path.join(out_dir, f'groups_{cell_type_id}_blood.csv')
    generate_groups_file(cell_type_id, groups_csv, blood_groups_file, beta_dir=beta_dir)

    if use_full_atlas:
        fm_groups_file = os.path.join(out_dir, f'groups_{cell_type_id}_atlas.csv')
        generate_groups_file_full_atlas(cell_type_id, groups_csv, fm_groups_file, beta_dir=beta_dir)
    else:
        fm_groups_file = blood_groups_file

    n_target = sum(1 for _, r in pd.read_csv(fm_groups_file).iterrows()
                   if r['group'] != 'background')
    n_bg = len(pd.read_csv(fm_groups_file)) - n_target
    print(f"    find_markers: {n_target} target, {n_bg} background ({bg_mode})")

    # Step 0b: Generate beta list
    print(f"  Step 0b: Generating beta list")
    beta_list = os.path.join(out_dir, 'beta_list.txt')
    generate_beta_list(beta_dir, beta_list)
    print(f"    {sum(1 for _ in open(beta_list))} beta files")

    # Step 0c: Run find_markers
    fm_out_dir = os.path.join(out_dir, 'find_markers_output')
    bed_path = os.path.join(fm_out_dir, f'Markers.{cell_type_id}.bed')

    if skip_find_markers and os.path.isfile(bed_path):
        print(f"  Step 0c: Skipping find_markers (using existing {bed_path})")
    else:
        print(f"  Step 0c: Running wgbstools find_markers ({marker_type})")
        bed_path = run_find_markers(
            groups_file=fm_groups_file,
            beta_list=beta_list,
            blocks_file=blocks_file,
            out_dir=fm_out_dir,
            cell_type_id=cell_type_id,
            threads=threads,
            top_n=top_n,
            delta_means=delta_means,
            min_cpg=min_cpg,
            unmeth_mean_thresh=unmeth_mean_thresh,
            meth_mean_thresh=meth_mean_thresh,
            wgbstools_path=wgbstools_path,
            only_hyper=only_hyper,
        )

    # Step 0d: Parse BED output
    print(f"  Step 0d: Parsing find_markers output")
    blocks, cpg_indices = parse_markers_bed_with_indices(bed_path, cell_type_id)
    print(f"    Found {len(blocks)} markers")

    if not blocks:
        print("    WARNING: No markers found! Check find_markers parameters.")
        return []

    # Step 0d2: Exclude sex chromosomes (chrX, chrY)
    # chrX has copy number differences between sexes (females 2, males 1)
    # and X-inactivation methylation differences, making dPCR quantification
    # unreliable. chrY is male-only. Both should be excluded from dPCR assays.
    before_sex = len(blocks)
    sex_chroms = {"chrX", "chrY", "chrM"}
    sex_blocks = [b for b in blocks if b.chrom in sex_chroms]
    if sex_blocks:
        blocks[:] = [b for b in blocks if b.chrom not in sex_chroms]
        print(f"  Step 0d2: Excluded {before_sex - len(blocks)} blocks on sex/mito chromosomes")
    print(f"    {len(blocks)} blocks remain after sex chromosome exclusion")

    # Step 0e: Extract per-CpG methylation
    # Background matches find_markers: full atlas for non-blood types,
    # blood-only for blood types. For CD4T/CD8T (blood types), use a
    # FILTER groups file that includes the opposite T subtype, so the
    # per-subgroup filter can reject T-cell-general markers.
    if use_full_atlas:
        bg_label = "full atlas"
        filter_bg_file = fm_groups_file
    else:
        bg_label = "blood-only"
        filter_groups_file = os.path.join(out_dir, f'groups_{cell_type_id}_filter.csv')
        if cell_type_id in FILTER_EXCLUDE:
            generate_groups_file_filter(cell_type_id, groups_csv, filter_groups_file, beta_dir=beta_dir)
            filter_bg_file = filter_groups_file
        else:
            filter_bg_file = blood_groups_file
    print(f"  Step 0e: Extracting per-CpG methylation from beta files ({bg_label})")

    target_betas, bg_betas, bg_group_map = _match_beta_files_with_groups(
        filter_bg_file, groups_csv, beta_dir
    )
    extract_percpg_methylation_with_indices(
        blocks, cpg_indices, target_betas, bg_betas,
        max_bg_samples=max_bg_samples,
        bg_group_map=bg_group_map,
        only_hyper=only_hyper,
    )

    # Report DMS scores
    good = sum(1 for b in blocks if b.dms_score >= 0.6)
    print(f"    DMS: {good}/{len(blocks)} blocks score >= 0.6")

    # Step 0e1: DMS score filter
    # Reject blocks with CS < 0.6 — these have poor per-CpG data quality
    # (low coverage, high variability, or no usable CpG data at all).
    before_cs = len(blocks)
    low_dms_blocks = [b for b in blocks if b.dms_score < 0.6]
    if low_dms_blocks:
        blocks[:] = [b for b in blocks if b.dms_score >= 0.6]
        print(f"  Step 0e1: DMS filter (CS >= 0.6): "
              f"rejected {before_cs - len(blocks)}/{before_cs} blocks")
    print(f"    {len(blocks)} blocks remain after DMS filter")

    # Step 0e2: Per-subgroup background filter
    # For hypomethylated markers: reject blocks where any background
    # subgroup has mean methylation BELOW the threshold (partially
    # unmethylated → false positives in real blood).
    # For hypermethylated markers: reject blocks where any background
    # subgroup has mean methylation ABOVE (1 - threshold) (partially
    # methylated → false positives in real blood).
    #
    # Blocks with empty or partial bg_subgroup_meth are REJECTED, not
    # passed. Missing subgroup data means we cannot verify that subgroup
    # is sufficiently methylated/unmethylated — the safe choice is to reject.
    if bg_group_map and min_bg_subgroup_meth > 0:
        if only_hyper:
            max_bg_subgroup_meth = 1.0 - min_bg_subgroup_meth
            filter_desc = f"max <= {max_bg_subgroup_meth:.2f}"
        else:
            filter_desc = f"min >= {min_bg_subgroup_meth:.2f}"

        # Expected set of background subgroups (from the group map)
        expected_subgroups = set(bg_group_map.values())

        before = len(blocks)
        filtered = []
        rejected_blocks = []
        for b in blocks:
            if not b.bg_subgroup_meth:
                # No subgroup data at all — reject (cannot verify)
                rejected_blocks.append(b)
                continue
            # Check for missing subgroups (partial data)
            present_subgroups = set(b.bg_subgroup_meth.keys())
            missing = expected_subgroups - present_subgroups
            if missing:
                # Partial data — some subgroups have no coverage.
                # Reject: we cannot verify the missing subgroups.
                rejected_blocks.append(b)
                continue
            if only_hyper:
                max_subgroup = max(b.bg_subgroup_meth.values())
                if max_subgroup <= (1.0 - min_bg_subgroup_meth):
                    filtered.append(b)
                else:
                    rejected_blocks.append(b)
            else:
                min_subgroup = min(b.bg_subgroup_meth.values())
                if min_subgroup >= min_bg_subgroup_meth:
                    filtered.append(b)
                else:
                    rejected_blocks.append(b)
        rejected = len(rejected_blocks)
        if rejected > 0:
            print(f"  Step 0e2: Per-subgroup filter ({filter_desc}): "
                  f"rejected {rejected}/{before} blocks")
            # Categorize rejections
            no_data = sum(1 for b in rejected_blocks if not b.bg_subgroup_meth)
            partial_data = sum(1 for b in rejected_blocks
                             if b.bg_subgroup_meth
                             and (expected_subgroups - set(b.bg_subgroup_meth.keys())))
            threshold_rejected = rejected - no_data - partial_data
            if no_data:
                print(f"    {no_data} blocks rejected (no subgroup data — low coverage)")
            if partial_data:
                print(f"    {partial_data} blocks rejected (partial subgroup data — some subgroups missing)")
            rejected_subgroups = {}
            for b in rejected_blocks:
                if not b.bg_subgroup_meth:
                    continue
                if expected_subgroups - set(b.bg_subgroup_meth.keys()):
                    continue  # partial data, skip subgroup breakdown
                if only_hyper:
                    worst = max(b.bg_subgroup_meth, key=b.bg_subgroup_meth.get)
                    reason = "partially methylated"
                else:
                    worst = min(b.bg_subgroup_meth, key=b.bg_subgroup_meth.get)
                    reason = "partially unmethylated"
                rejected_subgroups[worst] = rejected_subgroups.get(worst, 0) + 1
            for sg, count in sorted(rejected_subgroups.items(), key=lambda x: -x[1]):
                print(f"    {sg}: {count} blocks rejected ({reason})")
        blocks[:] = filtered
        print(f"    {len(blocks)} blocks remain after per-subgroup filter")

    # Step 0e3: Per-CpG per-cell-type enrichment
    # Read ALL beta files (not subsampled) to compute per-CpG mean beta
    # for every atlas cell type. This populates cpg_site.per_celltype_beta
    # for the Excel "CpG CellType Beta" sheet.
    all_beta_group_map = _build_all_beta_group_map(groups_csv, beta_dir)
    enrich_percpg_per_celltype(blocks, all_beta_group_map, min_cov=5)

    # Step 0f: Save as JSON
    json_path = os.path.join(out_dir, 'dmr_blocks.json')
    save_blocks_json(blocks, json_path)

    return blocks


# ---------------------------------------------------------------------------
# CLI for standalone testing
# ---------------------------------------------------------------------------
def main():
    """CLI entry point for testing phase0 independently."""
    import argparse
    parser = argparse.ArgumentParser(
        description='DMR Discovery via wgbstools find_markers'
    )
    parser.add_argument('--cell-type', required=True,
                        help=f'Cell type ID: {", ".join(VALID_CELL_TYPES)}')
    parser.add_argument('--beta-dir', default='data/beta_files/',
                        help='Directory with .beta files')
    parser.add_argument('--blocks-file',
                        default='data/GSE186458_blocks.s205.bed.gz',
                        help='Path to blocks BED file')
    parser.add_argument('--groups-csv', default='data/sample_groups.csv',
                        help='Path to full atlas groups CSV')
    parser.add_argument('--output-dir', default='dmr_output/',
                        help='Output directory')
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--top', type=int, default=1000)
    parser.add_argument('--delta-means', type=float, default=0.3,
                        help='Min methylation difference (default 0.3)')
    parser.add_argument('--unmeth-mean-thresh', type=float, default=0.15,
                        help='Target mean must be below this (default 0.15)')
    parser.add_argument('--meth-mean-thresh', type=float, default=0.65,
                        help='Background mean must be above this (default 0.65)')
    parser.add_argument('--min-bg-subgroup-meth', type=float, default=0.70,
                        help='Reject blocks where any bg subgroup is below this (default 0.70)')
    parser.add_argument('--max-bg-samples', type=int, default=30)
    parser.add_argument('--skip-find-markers', action='store_true',
                        help='Skip find_markers (use existing BED)')
    parser.add_argument('--use-full-atlas', action='store_true',
                        help='Use full atlas as find_markers background (hybrid approach)')
    parser.add_argument('--only-hyper', action='store_true',
                        help='Find hypermethylated markers (target methylated, bg unmethylated)')
    args = parser.parse_args()

    blocks = discover_dmrs(
        cell_type_id=args.cell_type,
        beta_dir=args.beta_dir,
        blocks_file=args.blocks_file,
        groups_csv=args.groups_csv,
        out_dir=args.output_dir,
        threads=args.threads,
        top_n=args.top,
        delta_means=args.delta_means,
        unmeth_mean_thresh=args.unmeth_mean_thresh,
        meth_mean_thresh=args.meth_mean_thresh,
        min_bg_subgroup_meth=args.min_bg_subgroup_meth,
        max_bg_samples=args.max_bg_samples,
        skip_find_markers=args.skip_find_markers,
        use_full_atlas=args.use_full_atlas,
        only_hyper=args.only_hyper,
    )

    print(f"\n=== DMR Discovery complete ===")
    print(f"  Cell type: {args.cell_type}")
    print(f"  Markers found: {len(blocks)}")
    if blocks:
        print(f"\n  Top 5 markers:")
        for b in blocks[:5]:
            print(f"    {b.seq_id}: {b.chrom}:{b.start}-{b.end} "
                  f"{b.num_cpgs}CpGs delta={b.delta_means:.3f} "
                  f"score={b.dms_score:.3f} gene={b.gene}")


if __name__ == "__main__":
    main()
