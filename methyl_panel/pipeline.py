#!/usr/bin/env python3
"""
Cell-Specific Epigenetic dPCR Assay Design Pipeline CLI

Modular step-by-step pipeline for designing bisulfite-PCR primers
for cell-type-specific methylation markers.

Usage:
    # DMR discovery from raw beta files (fully self-contained):
    python -m methyl_panel.pipeline --steps all --discover-dmrs \\
        --cell-type MONO --genome data/hg19/hg19.fa.gz \\
        --min-tm 58 --opt-tm 60 --max-tm 62 --output-dir results/

    # DMR discovery for any atlas cell type (e.g. Liver-Hepatocytes):
    python -m methyl_panel.pipeline --steps all --discover-dmrs \\
        --cell-type Liver-Hepatocytes --genome data/hg19/hg19.fa.gz \\
        --direction both --output-dir results/Liver-Hepatocytes/

    # Load DMR blocks from a pre-computed Excel file:
    python -m methyl_panel.pipeline --steps 1,2,3 --dmr-xlsx DMR.xlsx --genome hg19.fa.gz
    python -m methyl_panel.pipeline --steps all --output-dir results/
    python -m methyl_panel.pipeline --steps 3,4,5,6,7 --primers primers.json

Steps:
    0. Discover DMRs from raw beta files (wgbstools find_markers)  [--discover-dmrs]
    1. Load DMR blocks from Excel
    2. Bisulfite convert genomic sequences
    3. Design primers with Primer3
    4. 16-mer hash ePCR specificity screening
    5. Secondary structure screening
    6. Common SNP screening
    7. Primer-dimer prediction (DimerDetective)
    8. Generate XLSX output (U-assays format)
    9. Generate PDF output (U-assays format)
"""

import os
import sys
import json
import argparse
import time

import numpy as np
from typing import List, Optional, Dict

from methyl_panel.config import Primer3PlusConfig, PipelineConfig


def _validate_seq_id_uniqueness(blocks_data, source_path=""):
    """Check that all seq_ids in blocks_data are unique.

    Prevents use of stale dmr_blocks.json files generated before the v2.2.21
    seq_id fix. When direction='both', the old code assigned the same seq_id
    to hypo and hyper DMRs of the same rank (e.g. CellType_0001 for both),
    causing wrong DMR metadata lookups in Excel/PDF output. The fix appends
    _U (hypo) or _M (hyper) suffixes. This guard catches any file that still
    has duplicates and refuses to proceed.

    Args:
        blocks_data: List of block dicts with 'seq_id' keys.
        source_path: Path the data was loaded from (for error message).

    Raises:
        RuntimeError: If duplicate seq_ids are found.
    """
    seq_ids = [b["seq_id"] for b in blocks_data if "seq_id" in b]
    seen = set()
    duplicates = set()
    for sid in seq_ids:
        if sid in seen:
            duplicates.add(sid)
        seen.add(sid)

    if duplicates:
        location = f" from {source_path}" if source_path else ""
        raise RuntimeError(
            f"Duplicate seq_ids found{location}: {sorted(duplicates)[:5]}"
            f"{'...' if len(duplicates) > 5 else ''}. "
            f"This dmr_blocks.json was generated before the v2.2.21 seq_id fix. "
            f"Re-run from step 0 with --discover-dmrs to regenerate with "
            f"unique _U/_M suffixes, or delete the stale file."
        )


def step0_discover_dmrs(args):
    """Step 0: Discover DMRs from raw beta files using wgbstools find_markers.

    When --direction is 'both', runs discovery twice (hypo then hyper),
    merges results, and tags each block with its direction.
    """
    from methyl_panel.phase0_dmr_discovery import discover_dmrs, VALID_CELL_TYPES

    # Determine which directions to run
    direction = getattr(args, 'direction', 'both')
    if direction == 'both':
        directions = ['hypo', 'hyper']
    else:
        directions = [direction]

    all_blocks_data = []

    for d in directions:
        print(f"\n=== Step 0: DMR Discovery (wgbstools find_markers) — direction={d} ===")

        if not args.cell_type:
            raise ValueError(
                "--cell-type is required when using --discover-dmrs. "
                f"Valid options: {', '.join(VALID_CELL_TYPES)}"
            )

        blocks = discover_dmrs(
            cell_type_id=args.cell_type,
            beta_dir=args.beta_dir,
            blocks_file=args.blocks_file,
            groups_csv=args.groups_csv,
            out_dir=args.output_dir,
            threads=args.threads,
            top_n=args.top_markers,
            delta_means=args.delta_means,
            unmeth_mean_thresh=args.unmeth_mean_thresh,
            meth_mean_thresh=args.meth_mean_thresh,
            min_bg_subgroup_meth=args.min_bg_subgroup_meth,
            wgbstools_path=args.wgbstools_path,
            max_bg_samples=args.max_bg_samples,
            skip_find_markers=args.skip_find_markers,
            use_full_atlas=args.use_full_atlas,
            only_hyper=(d == 'hyper'),
        )

        # Add direction suffix to seq_id to avoid duplicates when direction='both'.
        # Without this, hypo rank-1 and hyper rank-1 both get seq_id "CellType_0001",
        # causing wrong DMR metadata lookups downstream.
        if len(directions) > 1:
            direction_suffix = '_U' if d == 'hypo' else '_M'
            for b in blocks:
                b.seq_id = f"{b.seq_id}{direction_suffix}"

        # Convert DMRBlock objects to the same dict format that step 1 produces
        # so that steps 2-9 can consume them identically
        for b in blocks:
            all_blocks_data.append({
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
                "direction": b.direction if b.direction else d,
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

        print(f"  DMR discovery complete ({d}): {len(blocks)} blocks for {args.cell_type}")

    # Save merged blocks to JSON for downstream steps
    out_path = os.path.join(args.output_dir, "dmr_blocks.json")
    with open(out_path, 'w') as f:
        json.dump(all_blocks_data, f, indent=2)
    print(f"\n  Total blocks (all directions): {len(all_blocks_data)}")
    print(f"  Saved to {out_path}")
    return all_blocks_data


def step1_load_dmrs(args):
    """Step 1: Load DMR blocks from Excel."""
    from methyl_panel.phase1_dmr_loader import load_dmr_blocks, filter_blocks_by_cpg_count

    print("\n=== Step 1: Loading DMR blocks ===")
    cell_type = args.cell_type if args.cell_type else None
    blocks = load_dmr_blocks(args.dmr_xlsx, cell_type)
    blocks = filter_blocks_by_cpg_count(blocks, args.min_cpg, args.preferred_min_cpg)

    print(f"Loaded {len(blocks)} blocks")
    if cell_type:
        print(f"  Filtered to cell type: {cell_type}")

    # Save to JSON for downstream steps
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
                 "delta_beta": c.delta_beta}
                for c in b.cpg_sites
            ],
        })

    out_path = os.path.join(args.output_dir, "dmr_blocks.json")
    with open(out_path, 'w') as f:
        json.dump(blocks_data, f, indent=2)
    print(f"  Saved to {out_path}")
    return blocks_data


def step2_bisulfite_convert(args, blocks_data=None):
    """Step 2: Bisulfite convert genomic sequences."""
    from methyl_panel.phase2_bisulfite_convert import convert_sequence, fetch_genomic_sequence

    print("\n=== Step 2: Bisulfite converting sequences ===")
    if blocks_data is None:
        blocks_path = os.path.join(args.output_dir, "dmr_blocks.json")
        with open(blocks_path) as f:
            blocks_data = json.load(f)
        _validate_seq_id_uniqueness(blocks_data, blocks_path)

    converted = []
    skipped = 0
    for i, block in enumerate(blocks_data):
        if args.max_blocks and (i - skipped) >= args.max_blocks:
            break

        # Fetch genomic sequence with flanking
        flank = args.flank
        try:
            seq = fetch_genomic_sequence(
                args.genome, block["chrom"],
                block["start"], block["end"], flank
            )
        except (RuntimeError, Exception) as e:
            print(f"  Skipping {block['seq_id']} ({block['chrom']}): {e}")
            skipped += 1
            continue

        if not seq or 'N' * 10 in seq:
            print(f"  Skipping {block['seq_id']} ({block['chrom']}): sequence not available")
            skipped += 1
            continue

        strands = convert_sequence(seq)

        converted.append({
            "seq_id": block["seq_id"],
            "cell_type_id": block["cell_type_id"],
            "direction": block.get("direction", "hypo"),
            "chrom": block["chrom"],
            "start": block["start"],
            "end": block["end"],
            "genomic_seq": seq,
            "strands": {
                "S": strands.sense,
                "SM": strands.sense_methylated,
                "SU": strands.sense_unmethylated,
                "A": strands.antisense,
                "AM": strands.antisense_methylated,
                "AU": strands.antisense_unmethylated,
            },
            "cpg_positions": strands.cpg_positions,
        })

        if (i + 1) % 50 == 0:
            print(f"  Converted {i+1}/{len(blocks_data)} blocks")

    out_path = os.path.join(args.output_dir, "converted_sequences.json")
    with open(out_path, 'w') as f:
        json.dump(converted, f)
    print(f"  Converted {len(converted)} blocks, saved to {out_path}")
    return converted


def step3_design_primers(args, converted_data=None):
    """Step 3: Design primers with Primer3."""
    from methyl_panel.phase3_primer3_design import design_primers_for_block

    print("\n=== Step 3: Designing primers ===")
    if converted_data is None:
        path = os.path.join(args.output_dir, "converted_sequences.json")
        with open(path) as f:
            converted_data = json.load(f)

    # Load config
    if args.settings:
        config = Primer3PlusConfig.from_primer3plus_file(args.settings)
    else:
        config = Primer3PlusConfig()

    # Apply salt preset (e.g. roche_dlc for Roche DLC conditions)
    config.salt_preset = args.salt_preset
    config.__post_init__()  # Re-apply preset to update salt concentrations

    # Set Tm (required)
    if args.min_tm is not None:
        config.primer_min_tm = args.min_tm
    if args.opt_tm is not None:
        config.primer_opt_tm = args.opt_tm
    if args.max_tm is not None:
        config.primer_max_tm = args.max_tm

    # Set product size range
    config.primer_product_size_range = [(args.product_size_min, args.product_size_max)]

    config.validate_tm()
    pipeline_config = PipelineConfig()
    pipeline_config.min_cpg_pair_total = args.min_cpg_pair_total
    pipeline_config.min_cpg_per_primer = args.min_cpg_per_primer
    pipeline_config.cpg_placement_mode = args.cpg_placement_mode

    all_primers = []
    for i, block in enumerate(converted_data):
        if args.max_blocks and i >= args.max_blocks:
            break

        primers = design_primers_for_block(
            block["genomic_seq"], config, pipeline_config,
            f"{block['seq_id']}_001", block["seq_id"], block["cell_type_id"],
            direction=block.get("direction", "hypo"),
            chrom=block.get("chrom"),
            dmr_start=block.get("start"),
            dmr_end=block.get("end"),
        )

        all_primers.extend(primers)

        if (i + 1) % 20 == 0:
            print(f"  Designed primers for {i+1}/{len(converted_data)} blocks ({len(all_primers)} total)")

    # Save to JSON
    primers_data = []
    for p in all_primers:
        primers_data.append({
            "assay_id": p.assay_id, "seq_id": p.seq_id,
            "cell_type_id": p.cell_type_id, "template_used": p.template_used,
            "direction": p.direction,
            "left_primer": p.left_primer, "right_primer": p.right_primer,
            "left_primer_display": p.left_primer_display,
            "right_primer_display": p.right_primer_display,
            "left_tm": p.left_tm, "right_tm": p.right_tm,
            "product_size": p.product_size,
            "c_total_tail": p.c_total_tail, "c_total": p.c_total,
            "left_c_total": p.left_c_total, "right_c_total": p.right_c_total,
            "left_c_tail": p.left_c_tail, "right_c_tail": p.right_c_tail,
            "sense_meth_mismatch_score": p.sense_meth_mismatch_score,
            "sense_unmeth_mismatch_score": p.sense_unmeth_mismatch_score,
            "anti_meth_mismatch_score": p.anti_meth_mismatch_score,
            "anti_unmeth_mismatch_score": p.anti_unmeth_mismatch_score,
            "left_gc_percent": p.left_gc_percent, "right_gc_percent": p.right_gc_percent,
            "penalty": p.penalty,
            "left_start": p.left_start, "left_len": p.left_len,
            "right_start": p.right_start, "right_len": p.right_len,
            "chrom": p.chrom,
            "dmr_start": p.dmr_start,
            "dmr_end": p.dmr_end,
        })

    out_path = os.path.join(args.output_dir, "primers.json")
    with open(out_path, 'w') as f:
        json.dump(primers_data, f, indent=2)
    print(f"  Designed {len(all_primers)} primer pairs, saved to {out_path}")
    return primers_data


def step4_epcr(args, primers_data=None):
    """Step 4: 16-mer hash-based ePCR specificity screening.

    Uses the 16-mer hash-based BiSearch ePCR reimplementation from
    phase4_bisulfite_ePCR.py. Screens each primer pair against the bisulfite-
    converted genome and produces a template-aware verdict.

    Requires per-chromosome .npy genome files (built by prepare_genomes.py).
    If the genome directory is not found, ePCR fields are set to None and
    the step is skipped gracefully.
    """
    from methyl_panel.phase4_bisulfite_ePCR import screen_pair

    print("\n=== Step 4: 16-mer hash ePCR specificity screening ===")
    if primers_data is None:
        path = os.path.join(args.output_dir, "primers.json")
        with open(path) as f:
            primers_data = json.load(f)

    genome_dir = getattr(args, 'epcr_genome_dir', None)

    if not genome_dir or not os.path.isdir(genome_dir):
        print(f"  ePCR genome directory not found: {genome_dir}")
        print("  Run prepare_genomes.py to build converted genome files.")
        print("  Skipping ePCR screening.")
        for p in primers_data:
            p["ePCR_verdict"] = None
            p["ePCR_product_count"] = None
            p["ePCR_products"] = None
            p["ePCR_error"] = None
            p["ePCR_elapsed_sec"] = None
    else:
        # Template-to-state mapping for template-aware verdict
        def _template_to_state(template):
            if not template:
                return None
            t = template.strip().upper()
            if t in ('SM', 'AM'):
                return 'methylated'
            if t in ('SU', 'AU'):
                return 'unmethylated'
            return None

        from methyl_panel.phase4_bisulfite_ePCR import (
            screen_chromosome, prepare_primer, primer_to_result,
            _resolve_chromosomes, MAX_MM, MAX_PCR_LEN,
        )

        screened = 0
        n_pass = 0
        n_fail = 0
        n_error = 0

        # Batched screening: prepare all primers once, then for each chromosome
        # build the CT/CT_M strand indexes ONCE and screen every still-active
        # primer against them. This is semantically identical to per-pair
        # screen_pair() calls (same early-FAIL deactivation, same product
        # dedup) but avoids rebuilding the 16-mer indexes per primer -- the
        # index build dominates runtime, so this is ~N_primers-fold faster.
        primer_objs = []
        for p in primers_data:
            fp = p.get("left_primer", "")
            rp = p.get("right_primer", "")
            assay_id = p.get("assay_id", "unknown")
            if not fp or not rp:
                primer_objs.append(None)
            else:
                primer_objs.append(prepare_primer(fp, rp, assay_id=assay_id,
                                                  max_mm=MAX_MM))

        chromosomes = _resolve_chromosomes(genome_dir, None)
        t_batch = time.time()
        for chrom in chromosomes:
            chrom_path = os.path.join(genome_dir, f"{chrom}.npy")
            if not os.path.exists(chrom_path):
                continue
            genome = np.load(chrom_path)
            screen_chromosome([po for po in primer_objs if po is not None],
                              genome, chrom, MAX_MM, MAX_PCR_LEN)
            del genome
            if not any(po is not None and po.active for po in primer_objs):
                break

        elapsed = time.time() - t_batch
        for p, po in zip(primers_data, primer_objs):
            template = p.get("template_used", "")
            if po is None:
                p["ePCR_verdict"] = "FAIL"
                p["ePCR_product_count"] = 0
                p["ePCR_products"] = ""
                p["ePCR_error"] = "missing primer sequence"
                p["ePCR_elapsed_sec"] = 0.0
                n_error += 1
                screened += 1
                continue

            result = primer_to_result(po, elapsed)

            if result.error:
                p["ePCR_verdict"] = "FAIL"
                p["ePCR_product_count"] = 0
                p["ePCR_products"] = ""
                p["ePCR_error"] = result.error
                p["ePCR_elapsed_sec"] = round(result.elapsed_sec, 2)
                n_error += 1
            elif result.product_count == 1:
                # Template-aware verdict: check if product is on intended state
                intended = _template_to_state(template)
                prod = result.products[0]
                if intended is not None and prod.state != intended:
                    p["ePCR_verdict"] = "FAIL"
                    p["ePCR_product_count"] = 1
                    p["ePCR_products"] = prod.coordinates()
                    p["ePCR_error"] = f"wrong_genome: product on {prod.state}, template expects {intended}"
                    p["ePCR_elapsed_sec"] = round(result.elapsed_sec, 2)
                    n_fail += 1
                else:
                    p["ePCR_verdict"] = "PASS"
                    p["ePCR_product_count"] = 1
                    p["ePCR_products"] = prod.coordinates()
                    p["ePCR_error"] = ""
                    p["ePCR_elapsed_sec"] = round(result.elapsed_sec, 2)
                    n_pass += 1
            else:
                p["ePCR_verdict"] = "FAIL"
                p["ePCR_product_count"] = result.product_count
                p["ePCR_products"] = "; ".join(
                    prod.coordinates() for prod in result.products)
                p["ePCR_error"] = "" if result.product_count > 1 else "no product"
                p["ePCR_elapsed_sec"] = round(result.elapsed_sec, 2)
                n_fail += 1

            screened += 1

        print(f"  Screened {screened} primer pairs:")
        print(f"    PASS: {n_pass}  FAIL: {n_fail}  ERROR: {n_error}")

    out_path = os.path.join(args.output_dir, "primers.json")
    with open(out_path, 'w') as f:
        json.dump(primers_data, f, indent=2)
    print(f"  ePCR screening complete, updated {out_path}")
    return primers_data


def step5_structure(args, primers_data=None):
    """Step 5: Secondary structure screening.

    Computes hairpin MFE, self-dimer ΔG, heterodimer ΔG, and amplicon
    secondary structure MFE. All metrics are INFORMATIONAL — no primers
    are rejected.
    """
    from methyl_panel.phase5_structure import screen_structure

    print("\n=== Step 5: Secondary structure screening ===")
    if primers_data is None:
        path = os.path.join(args.output_dir, "primers.json")
        with open(path) as f:
            primers_data = json.load(f)

    for p in primers_data:
        result = screen_structure(p["left_primer"], p["right_primer"])
        p["left_structure_mfe"] = result.left_mfe
        p["right_structure_mfe"] = result.right_mfe
        p["left_self_dimer_dg"] = result.left_self_dimer_dg
        p["right_self_dimer_dg"] = result.right_self_dimer_dg
        p["heterodimer_dg"] = result.heterodimer_dg
        p["amplicon_mfe"] = result.amplicon_mfe

    out_path = os.path.join(args.output_dir, "primers.json")
    with open(out_path, 'w') as f:
        json.dump(primers_data, f, indent=2)
    print(f"  Structure screening complete, updated {out_path}")
    return primers_data


def step6_snp(args, primers_data=None):
    """Step 6: Common SNP screening against dbSNP.

    Screens each primer pair for common variants (MAF >= threshold) that
    could cause allele dropout. Also screens the full amplicon for SNPs
    at measured CpG positions (C->T SNP at a CpG would always read as
    unmethylated after bisulfite conversion, confounding methylation calls).

    Primer-level scoring (common_variant_score, 2-digit code):
      00 = clean, 10 = body SNP, 01 = 3' SNP, 11 = both

    Amplicon-level annotations:
      amplicon_snp_count: total common SNPs in the full amplicon
      cpg_snp_count: SNPs at measured CpG positions (could confound methylation)
      cpg_snp_details: per-SNP detail with ref/alt/maf and is_c_to_t flag

    This is an INFORMATIONAL annotation -- it does not reject primers.
    """
    print("\n=== Step 6: Common SNP screening ===")
    if primers_data is None:
        path = os.path.join(args.output_dir, "primers.json")
        with open(path) as f:
            primers_data = json.load(f)

    dbsnp = getattr(args, 'dbsnp', None)
    maf_threshold = getattr(args, 'snp_maf_threshold', 0.01)
    flank = getattr(args, 'flank', 100)

    # Load cpg_sites from dmr_blocks.json for CpG-SNP detection
    cpg_positions_by_seqid = {}
    blocks_path = os.path.join(args.output_dir, "dmr_blocks.json")
    if os.path.exists(blocks_path):
        with open(blocks_path) as f:
            blocks = json.load(f)
        for b in blocks:
            sid = b.get("seq_id", "")
            positions = {c["position"] for c in b.get("cpg_sites", []) if "position" in c}
            if positions:
                cpg_positions_by_seqid[sid] = positions
        print(f"  Loaded CpG positions for {len(cpg_positions_by_seqid)} blocks from dmr_blocks.json")

    if not dbsnp or not os.path.exists(dbsnp):
        print(f"  No dbSNP VCF provided (--dbsnp) -- skipping SNP screening")
        for p in primers_data:
            p["common_variant_score"] = None
            p["snp_count"] = None
            p["amplicon_snp_count"] = None
            p["cpg_snp_count"] = None
            p["cpg_snp_details"] = []
    else:
        print(f"  Screening against dbSNP: {dbsnp}")
        print(f"  MAF threshold: {maf_threshold}")
        from methyl_panel.phase6_snp import batch_screen_snp

        results = batch_screen_snp(
            primers_data, dbsnp, maf_threshold=maf_threshold,
            flank=flank, cpg_positions_by_seqid=cpg_positions_by_seqid)

        screened = 0
        with_snps = 0
        with_cpg_snps = 0
        for p, result in zip(primers_data, results):
            p["common_variant_score"] = result.common_variant_score
            p["snp_count"] = result.snp_count
            p["amplicon_snp_count"] = result.amplicon_snp_count
            p["cpg_snp_count"] = result.cpg_snp_count
            p["cpg_snp_details"] = result.cpg_snp_details
            screened += 1
            if result.snp_count > 0:
                with_snps += 1
            if result.cpg_snp_count > 0:
                with_cpg_snps += 1

        print(f"  Screened {screened} primer pairs: {with_snps} with primer SNPs, "
              f"{with_cpg_snps} with CpG-affecting SNPs")

    out_path = os.path.join(args.output_dir, "primers.json")
    with open(out_path, 'w') as f:
        json.dump(primers_data, f, indent=2)
    print(f"  SNP screening complete, updated {out_path}")
    return primers_data


def step6b_repeats(args, primers_data=None):
    """Step 6b: Segmental duplication and repeat overlap screening.

    Checks primer binding sites AND the full amplicon against UCSC segmental
    duplications and RepeatMasker tracks. INFORMATIONAL only -- does not reject primers.

    Primer-level: segdup_overlap, repeat_overlap, repeat_family
    Amplicon-level: amplicon_segdup_overlap, amplicon_repeat_overlap, amplicon_repeat_family
    """
    print("\n=== Step 6b: Segmental duplication & repeat screening ===")
    if primers_data is None:
        path = os.path.join(args.output_dir, "primers.json")
        with open(path) as f:
            primers_data = json.load(f)

    segdup_bed = getattr(args, 'segdup_bed', None)
    repeat_bed = getattr(args, 'repeat_bed', None)
    flank = getattr(args, 'flank', 100)

    if not segdup_bed and not repeat_bed:
        print("  No segdup/repeat BED files provided -- skipping")
        for p in primers_data:
            p["segdup_overlap"] = None
            p["repeat_overlap"] = None
            p["repeat_family"] = None
            p["amplicon_segdup_overlap"] = None
            p["amplicon_repeat_overlap"] = None
            p["amplicon_repeat_family"] = None
    else:
        from methyl_panel.phase8_repeats import batch_screen_repeats

        results = batch_screen_repeats(
            primers_data, segdup_bed, repeat_bed, flank=flank)

        screened = 0
        segdup_count = 0
        repeat_count = 0
        amp_segdup_count = 0
        amp_repeat_count = 0
        for p, result in zip(primers_data, results):
            p["segdup_overlap"] = result.segdup_overlap
            p["repeat_overlap"] = result.repeat_overlap
            p["repeat_family"] = result.repeat_family
            p["amplicon_segdup_overlap"] = result.amplicon_segdup_overlap
            p["amplicon_repeat_overlap"] = result.amplicon_repeat_overlap
            p["amplicon_repeat_family"] = result.amplicon_repeat_family
            screened += 1
            if result.segdup_overlap:
                segdup_count += 1
            if result.repeat_overlap:
                repeat_count += 1
            if result.amplicon_segdup_overlap:
                amp_segdup_count += 1
            if result.amplicon_repeat_overlap:
                amp_repeat_count += 1

        print(f"  Screened {screened} pairs:")
        print(f"    Primer-level:  {segdup_count} segdup, {repeat_count} repeat overlaps")
        print(f"    Amplicon-level: {amp_segdup_count} segdup, {amp_repeat_count} repeat overlaps")

    out_path = os.path.join(args.output_dir, "primers.json")
    with open(out_path, 'w') as f:
        json.dump(primers_data, f, indent=2)
    print(f"  Repeat screening complete, updated {out_path}")
    return primers_data


def step7_dimer(args, primers_data=None):
    """Step 7: Primer-dimer prediction (DimerDetective).

    Computes end_min_dg (3' end stability), full heterodimer ΔG, and
    homodimer ΔG. All metrics are INFORMATIONAL — no primers are rejected.
    """
    from methyl_panel.phase7_dimer import predict_dimer

    print("\n=== Step 7: DimerDetective primer-dimer prediction ===")
    if primers_data is None:
        path = os.path.join(args.output_dir, "primers.json")
        with open(path) as f:
            primers_data = json.load(f)

    for p in primers_data:
        result = predict_dimer(p["left_primer"], p["right_primer"])
        p["primer_dimer_prediction"] = result.risk_tier
        p["primer_dimer_end_min_dg"] = result.end_min_dg
        p["heterodimer_full_dg"] = result.heterodimer_full_dg
        p["left_homodimer_dg"] = result.left_homodimer_dg
        p["right_homodimer_dg"] = result.right_homodimer_dg

    out_path = os.path.join(args.output_dir, "primers.json")
    with open(out_path, 'w') as f:
        json.dump(primers_data, f, indent=2)
    print(f"  Dimer prediction complete, updated {out_path}")
    return primers_data


def step8_xlsx(args, primers_data=None):
    """Step 8: Generate XLSX output."""
    from methyl_panel.output_xlsx import write_xlsx
    from methyl_panel.phase3_primer3_design import PrimerPair

    print("\n=== Step 8: Generating XLSX output ===")
    if primers_data is None:
        path = os.path.join(args.output_dir, "primers.json")
        with open(path) as f:
            primers_data = json.load(f)

    # Convert dicts back to PrimerPair objects
    pairs = [PrimerPair.from_dict(p) for p in primers_data]

    # Load DMR blocks for WGBS data coupling
    blocks_path = os.path.join(args.output_dir, "dmr_blocks.json")
    dmr_blocks = None
    if os.path.exists(blocks_path):
        with open(blocks_path) as f:
            dmr_blocks = json.load(f)
        _validate_seq_id_uniqueness(dmr_blocks, blocks_path)
        print(f"  Loaded {len(dmr_blocks)} DMR blocks for WGBS data coupling")

    out_path = os.path.join(args.output_dir, "primer_assays.xlsx")
    write_xlsx(pairs, out_path, dmr_blocks=dmr_blocks)
    print(f"  XLSX saved to {out_path}")


def step9_pdf(args, primers_data=None):
    """Step 9: Generate PDF output."""
    from methyl_panel.output_pdf import write_pdf
    from methyl_panel.phase3_primer3_design import PrimerPair

    print("\n=== Step 9: Generating PDF output ===")
    if primers_data is None:
        path = os.path.join(args.output_dir, "primers.json")
        with open(path) as f:
            primers_data = json.load(f)

    pairs = [PrimerPair.from_dict(p) for p in primers_data]

    out_path = os.path.join(args.output_dir, "primer_assays.pdf")
    write_pdf(pairs, output_path=out_path)
    print(f"  PDF saved to {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Cell-Specific Epigenetic dPCR Assay Design Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--steps", default="all",
                        help="Comma-separated step numbers (e.g. 0,2,3 or 1,2,3) or 'all'")
    parser.add_argument("--output-dir", default="results/",
                        help="Output directory for results")
    parser.add_argument("--dmr-xlsx", default="data/WBC_Panel_Top200_v7.9.xlsx",
                        help="Path to DMR per-CpG Excel file (for step 1)")
    parser.add_argument("--genome", default="data/hg19/hg19.fa.gz",
                        help="Path to genome FASTA (can be .gz)")
    parser.add_argument("--settings", help="Primer3Plus settings file")
    parser.add_argument("--cell-type",
                        help="Cell type ID. Supports the 7 merged blood types "
                             "(MONO, BCELL, NK, GRAN, CD3T, CD8T, CD4T) and all "
                             "82 Loyfer atlas groups (e.g. Liver-Hepatocytes, "
                             "Cortex-Neuron, Blood-T-CD4). Required with --discover-dmrs.")
    parser.add_argument("--min-cpg", type=int, default=5, help="Min CpGs per block")
    parser.add_argument("--preferred-min-cpg", type=int, default=7)
    parser.add_argument("--min-tm", type=float, help="Min Tm (°C) — REQUIRED")
    parser.add_argument("--opt-tm", type=float, help="Opt Tm (°C) — REQUIRED")
    parser.add_argument("--max-tm", type=float, help="Max Tm (°C) — REQUIRED")
    parser.add_argument("--salt-preset", default="primer3plus",
                        choices=["primer3plus", "roche_dlc"],
                        help="Salt conditions for Tm calculation (default: primer3plus). "
                             "Use 'roche_dlc' for Roche DLC conditions (100mM Na, 4.5mM Mg).")
    parser.add_argument("--product-size-min", type=int, default=60,
                        help="Minimum product size in bp (default 60)")
    parser.add_argument("--product-size-max", type=int, default=150,
                        help="Maximum product size in bp (default 150)")
    parser.add_argument("--min-cpg-pair-total", type=int, default=4,
                        help="Min total CpGs across both primers (default 4). "
                             "Lower to 3 for sparse-CpG DMRs.")
    parser.add_argument("--min-cpg-per-primer", type=int, default=2,
                        help="Min CpGs in at least one primer (default 2). "
                             "Lower to 1 for sparse-CpG DMRs.")
    parser.add_argument("--cpg-placement-mode", default="strict_terminal",
                        choices=["strict_terminal", "three_prime_half"],
                        help="CpG discrimination region: strict_terminal (last 5 nt, "
                             "default) or three_prime_half (last 50%% of primer length, "
                             "MSP-style)")
    parser.add_argument("--flank", type=int, default=100, help="Flanking bp around DMR")
    parser.add_argument("--max-blocks", type=int, help="Max blocks to process (for testing)")
    parser.add_argument("--epcr-genome-dir", default="data/hg19_bs/",
                        help="Directory with per-chromosome .npy genome files for ePCR screening "
                             "(built by prepare_genomes.py)")
    parser.add_argument("--dbsnp", help="Path to dbSNP VCF (for SNP screening)")
    parser.add_argument("--snp-maf-threshold", type=float, default=0.01,
                        help="Minimum MAF to consider a SNP 'common' (default 0.01 = 1%%)")
    parser.add_argument("--segdup-bed", help="Path to segmental duplications BED (tabix-indexed)")
    parser.add_argument("--repeat-bed", help="Path to RepeatMasker BED (tabix-indexed)")

    # --- Step 0: DMR Discovery arguments ---
    parser.add_argument("--discover-dmrs", action="store_true",
                        help="Run Step 0 (DMR discovery from raw beta files) "
                             "instead of Step 1 (load from Excel)")
    parser.add_argument("--beta-dir", default="data/beta_files/",
                        help="Directory containing .beta files (for --discover-dmrs)")
    parser.add_argument("--blocks-file", default="data/GSE186458_blocks.s205.bed.gz",
                        help="Path to blocks BED file (for --discover-dmrs)")
    parser.add_argument("--groups-csv", default="data/sample_groups.csv",
                        help="Path to full atlas groups CSV (for --discover-dmrs)")
    parser.add_argument("--wgbstools-path",
                        help="Path to wgbstools executable (auto-detected if omitted)")
    parser.add_argument("--threads", type=int, default=2,
                        help="Number of threads for find_markers")
    parser.add_argument("--top-markers", type=int, default=1000,
                        help="Number of top markers per cell type")
    parser.add_argument("--delta-means", type=float, default=0.3,
                        help="Min methylation difference for find_markers (default 0.3)")
    parser.add_argument("--unmeth-mean-thresh", type=float, default=0.15,
                        help="Target mean must be below this (default 0.15)")
    parser.add_argument("--meth-mean-thresh", type=float, default=0.65,
                        help="Background mean must be above this (default 0.65)")
    parser.add_argument("--min-bg-subgroup-meth", type=float, default=0.70,
                        help="Reject blocks where any bg subgroup is below this (default 0.70)")
    parser.add_argument("--max-bg-samples", type=int, default=30,
                        help="Max background samples for per-CpG extraction")
    parser.add_argument("--skip-find-markers", action="store_true",
                        help="Skip find_markers (reuse existing BED output)")
    parser.add_argument("--use-full-atlas", action="store_true",
                        help="Use full atlas (all 207 samples) as find_markers background "
                             "to get more candidates, then apply blood-only per-subgroup "
                             "filter for quality control (hybrid approach)")
    parser.add_argument("--only-hyper", action="store_true",
                        help="[Deprecated — use --direction hyper instead] "
                             "Find hypermethylated markers (target methylated, background "
                             "unmethylated) instead of hypomethylated markers")
    parser.add_argument("--direction", default="both",
                        choices=["hypo", "hyper", "both"],
                        help="DMR direction: hypo (target unmethylated, bg methylated), "
                             "hyper (target methylated, bg unmethylated), or both "
                             "(default: both). Supersedes --only-hyper.")

    args = parser.parse_args()

    # --direction supersedes --only-hyper (kept for backward compatibility)
    if args.only_hyper and args.direction == "both":
        args.direction = "hyper"

    os.makedirs(args.output_dir, exist_ok=True)

    # If --discover-dmrs is set, replace step 1 with step 0
    if args.discover_dmrs:
        if args.steps == "all":
            steps = [0, 2, 3, 4, 5, 6, 61, 7, 8, 9]
        else:
            steps = [int(s.strip()) for s in args.steps.split(",")]
            # If user specified step 1 with --discover-dmrs, replace with 0
            steps = [0 if s == 1 else s for s in steps]
    else:
        # Parse steps normally
        if args.steps == "all":
            steps = [1, 2, 3, 4, 5, 6, 61, 7, 8, 9]
        else:
            steps = [int(s.strip()) for s in args.steps.split(",")]

    # Validate Tm is set for primer design steps
    if any(s in steps for s in [3]) and not all([args.min_tm, args.opt_tm, args.max_tm]):
        if not args.settings:
            parser.error("Tm values (--min-tm, --opt-tm, --max-tm) are REQUIRED for primer design. "
                         "Or provide --settings file with Tm values.")

    # Execute steps
    blocks_data = None
    converted_data = None
    primers_data = None

    for step in steps:
        if step == 0:
            blocks_data = step0_discover_dmrs(args)
        elif step == 1:
            blocks_data = step1_load_dmrs(args)
        elif step == 2:
            converted_data = step2_bisulfite_convert(args, blocks_data)
        elif step == 3:
            primers_data = step3_design_primers(args, converted_data)
        elif step == 4:
            primers_data = step4_epcr(args, primers_data)
        elif step == 5:
            primers_data = step5_structure(args, primers_data)
        elif step == 6:
            primers_data = step6_snp(args, primers_data)
        elif step == 61:
            primers_data = step6b_repeats(args, primers_data)
        elif step == 7:
            primers_data = step7_dimer(args, primers_data)
        elif step == 8:
            step8_xlsx(args, primers_data)
        elif step == 9:
            step9_pdf(args, primers_data)

    print("\n=== Pipeline complete ===")


if __name__ == "__main__":
    main()
