#!/usr/bin/env python3
"""
plant_design.py - Context-aware (CG/CHG/CHH) MSP assay design for plants.

Adapts the mammalian CpG-only pipeline to plant methylomes where
discriminating sites can be in any sequence context (Kawakatsu 2016:
columella specificity is carried almost entirely by CHH).

Template model (per DMR block, discriminating sites = PerCpG rows):
  K_M ("methylated" template)   = CpG C's + discriminating CH sites kept as C
  K_U ("unmethylated" template) = CpG C's only; discriminating CH sites -> T
  SM/AM are built from K_M, SU/AU from K_U (same naming as the mammalian
  pipeline so downstream target-detecting logic is unchanged).

Strand handling (mc_calls '-' rows record the motif's last base; the loader
decodes true cytosine plus-coordinates):
  '+' strand site at plus-coord x  -> top-strand C at x  (kept in CT mask)
  '-' strand site at plus-coord x  -> top-strand G at x   (kept in GA mask);
      for the bottom-strand template the C sits at len-1-x, which
      bisulfite_convert_bottom reaches by passing the "C-analog" x-1.

ePCR genome models (per cell-type comparison, direction-independent):
  state_T (target)     = CpG kept + ALL DMR-block disc sites kept as C
                         (--dmr-xlsx/--cell-type-id). This exactly matches the
                         K_M design template, the plant analog of the mammalian
                         fully-methylated CT_M ePCR state. Without --dmr-xlsx,
                         falls back to WGBS-informed T (CH sites with
                         target-group mean beta >= threshold), which does NOT
                         match the design template and is not recommended.
  state_B (background) = CpG kept + CH sites with any background-group mean beta >= threshold
  state_U (fully converted) = nothing kept (extra conservatism)
  PASS = exactly 1 locus whose best T-state hit has strictly fewer mismatches
  than any B/U hit at that locus.
  Design post-filter: primer pairs must cover >= min-disc-in-primers
  discriminating sites inside their primer spans (a pair whose primers contain
  no disc site cannot discriminate methylation and is discarded).

Usage:
  python -m methyl_panel.plant_design design --dmr-xlsx dmr_blocks.xlsx \
      --cell-type Columella_hyper --fasta-dir /genomes/TAIR10/fa \
      --output-dir run/ [--flank 100] [primer3 settings]
  python -m methyl_panel.plant_design epcr --primers run/primers.json \
      --genome-dir /genomes/TAIR10_bs --samples samples.tsv \
      --cell-type Columella
Then finish with: python -m methyl_panel.pipeline --steps 5,7,8 ...
"""

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

from .phase2_bisulfite_convert import (
    BisulfiteStrands, bisulfite_convert_top, bisulfite_convert_bottom,
    find_cpg_positions, fetch_genomic_sequence, reverse_complement,
)
from .phase3_primer3_design import (
    Primer3PlusConfig, PipelineConfig, design_primers_for_block,
)
from .phase4_bisulfite_ePCR import (
    MAX_MM, MAX_PCR_LEN,
    build_index, search_and_verify, find_products_vec, prepare_primer,
    convert_genome_ct,
)
from .phase0_generic import load_arabidopsis


# ---------------------------------------------------------------------------
# Stage A: context-aware conversion + primer3 design (replaces steps 2+3)
# ---------------------------------------------------------------------------

def build_plant_strands(seq, disc_sites):
    """disc_sites: list of (rel0, strand) for discriminating sites,
    rel0 = 0-based offset in seq. Returns BisulfiteStrands."""
    cpg = find_cpg_positions(seq)
    k_top = set(cpg)
    k_bot = set(cpg)  # convert_bottom expects C-analog positions
    for rel0, strand in disc_sites:
        if strand == "+":
            k_top.add(rel0)
        elif strand == "-":
            # bottom C at len-1-rel0; convert_bottom maps input i -> len-i-2,
            # so pass the C-analog rel0-1 (base before the top-strand G)
            k_bot.add(rel0 - 1)
    k_top = sorted(k_top)
    k_bot = sorted(k_bot)
    return BisulfiteStrands(
        sense=seq,
        sense_methylated=bisulfite_convert_top(seq, k_top, methylated=True),
        sense_unmethylated=bisulfite_convert_top(seq, cpg, methylated=True),
        antisense=reverse_complement(seq),
        antisense_methylated=bisulfite_convert_bottom(seq, k_bot, methylated=True),
        antisense_unmethylated=bisulfite_convert_bottom(seq, cpg, methylated=True),
        cpg_positions=sorted(r for r, _ in disc_sites),
    )


def run_design(args):
    xl = pd.ExcelFile(args.dmr_xlsx)
    bs = xl.parse("Block_Summary")
    bs = bs[bs["cell_type_id"] == args.cell_type]
    if bs.empty:
        raise SystemExit(f"No blocks for cell type {args.cell_type}")
    want = f"PerCpG_{args.cell_type}"[:31]
    if want in xl.sheet_names:
        pc_sheet = want
    else:
        pc_sheet = next(s for s in xl.sheet_names if s.startswith("PerCpG_"))
    pc_all = xl.parse(pc_sheet)

    config = Primer3PlusConfig()
    config.primer_min_tm = args.min_tm
    config.primer_opt_tm = args.opt_tm
    config.primer_max_tm = args.max_tm
    config.primer_product_size_range = [(args.product_size_min,
                                         args.product_size_max)]
    config.validate_tm()
    pipeline_config = PipelineConfig()
    pipeline_config.min_cpg_pair_total = args.min_cpg_pair_total
    pipeline_config.min_cpg_per_primer = args.min_cpg_per_primer
    pipeline_config.cpg_placement_mode = args.cpg_placement_mode

    os.makedirs(args.output_dir, exist_ok=True)
    converted, all_primers, skipped = [], [], 0
    for row in bs.itertuples(index=False):
        seq_id = row.seq_id
        chrom = str(row.Chr)
        start, end = int(row.Start), int(row.End)
        fa = os.path.join(args.fasta_dir, chrom + ".fa")
        if not os.path.exists(fa):
            print(f"  Skipping {seq_id}: no FASTA for chrom {chrom}")
            skipped += 1
            continue
        try:
            seq = fetch_genomic_sequence(fa, chrom, start, end, args.flank)
        except Exception as e:
            print(f"  Skipping {seq_id}: {e}")
            skipped += 1
            continue
        if not seq or "N" * 10 in seq:
            print(f"  Skipping {seq_id}: sequence not available")
            skipped += 1
            continue
        pc = pc_all[pc_all["seq_id"] == seq_id]
        disc = []
        for p in pc.itertuples(index=False):
            rel0 = int(p.cpg_position) - start + args.flank
            if 0 <= rel0 < len(seq):
                disc.append((rel0, str(getattr(p, "strand", "."))))
        if not disc:
            print(f"  Skipping {seq_id}: no discriminating sites in sequence")
            skipped += 1
            continue
        strands = build_plant_strands(seq, disc)
        primers = design_primers_for_block(
            seq, config, pipeline_config,
            f"{seq_id}_001", seq_id, args.cell_type,
            direction=str(row.Direction), chrom=chrom,
            dmr_start=start, dmr_end=end,
            strands_override=strands,
        )
        all_primers.extend(primers)
        converted.append({
            "seq_id": seq_id, "cell_type_id": args.cell_type,
            "direction": str(row.Direction), "chrom": chrom,
            "start": start, "end": end, "genomic_seq": seq,
            "context": str(getattr(row, "context", "CH")),
            "discriminating_sites": [
                {"rel0": r, "strand": s} for r, s in disc],
        })

    with open(os.path.join(args.output_dir, "converted_sequences.json"),
              "w") as f:
        json.dump(converted, f)

    primers_data = [{
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
        "left_gc_percent": p.left_gc_percent,
        "right_gc_percent": p.right_gc_percent,
        "penalty": p.penalty,
        "left_start": p.left_start, "left_len": p.left_len,
        "right_start": p.right_start, "right_len": p.right_len,
        "chrom": p.chrom, "dmr_start": p.dmr_start, "dmr_end": p.dmr_end,
    } for p in all_primers]

    # Post-filter: a pair whose primers cover no discriminating site cannot
    # distinguish methylated from unmethylated template (SM/SU primers would
    # be identical). Require >= min_disc disc sites inside primer spans.
    disc_by_seq = {c["seq_id"]: (c["discriminating_sites"], len(c["genomic_seq"]))
                   for c in converted}
    kept_primers, n_nodisc = [], 0
    for p in primers_data:
        disc, seqlen = disc_by_seq[p["seq_id"]]
        spans = [(p["left_start"], p["left_start"] + p["left_len"]),
                 (p["right_start"], p["right_start"] + p["right_len"])]
        n = 0
        for site in disc:
            rel0, strand = site["rel0"], site.get("strand", "+")
            # sense templates: '+' sites sit at rel0; antisense templates:
            # '-' sites sit at len-1-rel0 ('+' sites are CpGs on the bottom
            # strand and are kept in both SM/SU, so they never discriminate)
            pos = (rel0 if (p["template_used"].startswith("S")
                            and strand == "+")
                   else (seqlen - 1 - rel0
                         if (p["template_used"].startswith("A")
                             and strand == "-") else None))
            if pos is not None and any(a <= pos < b for a, b in spans):
                n += 1
        if n >= args.min_disc_in_primers:
            p["disc_sites_in_primers"] = n
            kept_primers.append(p)
        else:
            n_nodisc += 1
    primers_data = kept_primers

    out_primers = os.path.join(args.output_dir, "primers.json")
    with open(out_primers, "w") as f:
        json.dump(primers_data, f, indent=2)
    print(f"Designed {len(primers_data)} primer pairs for "
          f"{len(converted)} blocks ({skipped} skipped, "
          f"{n_nodisc} dropped: no disc site in primers); {out_primers}")


# ---------------------------------------------------------------------------
# Stage B: masked-genome ePCR (replaces step 4)
# ---------------------------------------------------------------------------

class PlantPrimer:
    """Wraps phase4 PrimerInfo with plant state counters.

    Products accumulate across chromosomes into per-locus records keyed by
    (chrom, strand-tag, fp_pos, rp_pos); a locus found in several methylation
    states counts once, keeping the best T-state and best B/U-state
    mismatches seen for the verdict.
    """
    def __init__(self, info):
        self.info = info
        self.loci = {}   # (chrom, tag, fp_pos, rp_pos) -> record dict

    @property
    def error(self):
        return self.info.error

    def add_product(self, state, prod):
        key = (prod.chrom, prod.strand, prod.fp_pos, prod.rp_pos)
        d = self.loci.setdefault(key, {"t_mm": None, "bu_mm": None,
                                       "prod": None})
        if state == "T":
            if d["t_mm"] is None or prod.total_mm < d["t_mm"]:
                d["t_mm"] = prod.total_mm
            if d["prod"] is None:
                d["prod"] = prod
        else:
            if d["bu_mm"] is None or prod.total_mm < d["bu_mm"]:
                d["bu_mm"] = prod.total_mm

    @property
    def n_T(self):
        return sum(1 for d in self.loci.values() if d["t_mm"] is not None)

    @property
    def n_BU(self):
        return sum(1 for d in self.loci.values() if d["bu_mm"] is not None)


def masked_ct(genome, keep_c):
    conv = genome.copy()
    is_c = conv == 1
    conv[is_c & ~keep_c] = 3
    return conv


def cpg_masks(genome):
    keep_c = np.zeros(len(genome), dtype=bool)
    keep_g = np.zeros(len(genome), dtype=bool)
    cpg = (genome[:-1] == 1) & (genome[1:] == 2)
    keep_c[:-1] |= cpg
    keep_g[1:] |= cpg
    return keep_c, keep_g


def _group_selected_sites(mats, meth_beta, min_cov):
    """mats: list of (pos, meth, tot, ctx, strand) for one group.
    Returns (plus_pos, minus_pos) with group mean beta >= meth_beta."""
    if not mats:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    union = np.unique(np.concatenate([m[0] for m in mats]))
    beta = np.zeros((len(mats), len(union)), dtype=np.float32)
    for i, (pos, meth, tot, _c, _s) in enumerate(mats):
        ok = tot >= min_cov
        b = np.where(ok, meth / np.maximum(tot, 1), 0.0)
        beta[i, np.searchsorted(union, pos)] = b
    mb = beta.mean(axis=0)
    sel = union[mb >= meth_beta]
    smap = dict(zip(mats[0][0].tolist(), np.asarray(mats[0][4])))
    st = np.array([smap.get(int(p), "+") for p in sel])
    return sel[st == "+"], sel[st == "-"]


def build_ch_masks(mc_by_sample, groups, target_group, chroms,
                   meth_beta=0.65, min_cov=3):
    """Per-chrom dict state -> (plus_positions, minus_positions)."""
    sample_ids = list(mc_by_sample)
    tgt_ids = [s for s in sample_ids if groups[s] == target_group]
    bg_groups = sorted({g for s, g in groups.items() if g != target_group})
    masks = {}
    for chrom in chroms:
        out = {}
        # target state
        plus, minus = [], []
        for ctxname in ("CHH", "CHG"):
            key = f"{chrom}:{ctxname}"
            mats = [mc_by_sample[s][key] for s in tgt_ids
                    if key in mc_by_sample[s]]
            p, m = _group_selected_sites(mats, meth_beta, min_cov)
            plus.append(p)
            minus.append(m)
        out["T"] = (np.unique(np.concatenate(plus)) if plus
                    else np.zeros(0, dtype=np.int64),
                    np.unique(np.concatenate(minus)) if minus
                    else np.zeros(0, dtype=np.int64))
        # background state: union over background groups
        plus, minus = [], []
        for ctxname in ("CHH", "CHG"):
            key = f"{chrom}:{ctxname}"
            for g in bg_groups:
                gids = [s for s in sample_ids if groups[s] == g]
                mats = [mc_by_sample[s][key] for s in gids
                        if key in mc_by_sample[s]]
                p, m = _group_selected_sites(mats, meth_beta, min_cov)
                plus.append(p)
                minus.append(m)
        out["B"] = (np.unique(np.concatenate(plus)) if plus
                    else np.zeros(0, dtype=np.int64),
                    np.unique(np.concatenate(minus)) if minus
                    else np.zeros(0, dtype=np.int64))
        masks[chrom] = out
    return masks


def load_block_sites(dmr_xlsx, cell_type_id, chroms):
    """Per-chrom (plus_pos, minus_pos) 0-based plus-coord arrays of ALL
    PerCpG sites of the cell type's DMR blocks (any beta/coverage)."""
    xl = pd.ExcelFile(dmr_xlsx)
    sheet = f"PerCpG_{cell_type_id}"[:31]
    if sheet not in xl.sheet_names:
        sheet = next(s for s in xl.sheet_names if s.startswith("PerCpG_"))
    pc = xl.parse(sheet)
    bs = xl.parse("Block_Summary")
    chrom_of = dict(zip(bs["seq_id"].astype(str), bs["Chr"].astype(str)))
    ch = pc["seq_id"].astype(str).map(chrom_of)
    out = {}
    for chrom in chroms:
        sub = pc[ch == chrom]
        st = sub["strand"].astype(str).str[0]
        out[chrom] = (
            sub.loc[st == "+", "cpg_position"].to_numpy(dtype=np.int64),
            sub.loc[st == "-", "cpg_position"].to_numpy(dtype=np.int64),
        )
    return out


def run_epcr(args):
    with open(args.primers) as f:
        primers_data = json.load(f)
    if not primers_data:
        raise SystemExit("No primers to screen")

    s = pd.read_csv(args.samples, sep="\t")
    groups = dict(zip(s["sample_id"], s["group"]))

    print("Loading mc_calls for CH masks ...", file=sys.stderr)
    mc_by_sample = {}
    for rec in s.itertuples(index=False):
        mc_by_sample[rec.sample_id] = load_arabidopsis(rec.path)

    chroms = sorted(f[:-4] for f in os.listdir(args.genome_dir)
                    if f.endswith(".npy"))
    use_block_mask = (getattr(args, "dmr_xlsx", None) is not None
                      and getattr(args, "cell_type_id", None) is not None)
    if use_block_mask:
        print("T-state = CpG + DMR-block sites (template-consistent mode)",
              file=sys.stderr)
        block_sites = load_block_sites(args.dmr_xlsx, args.cell_type_id,
                                       chroms)
    else:
        print("WARNING: --dmr-xlsx/--cell-type-id not given; falling back to "
              "WGBS-informed T mask which may NOT match the design template",
              file=sys.stderr)
        block_sites = {}
    print("Building CH keep masks ...", file=sys.stderr)
    ch_masks = build_ch_masks(mc_by_sample, groups, args.cell_type, chroms,
                              meth_beta=args.meth_beta, min_cov=args.min_cov)

    # Add block discriminating sites (all of them, regardless of beta) to
    # the target-state mask: the design template keeps C at every block
    # site, so the ePCR target genome must too.
    if args.dmr_xlsx:
        xl = pd.ExcelFile(args.dmr_xlsx)
        ct_id = args.cell_type_id or args.cell_type
        want = f"PerCpG_{ct_id}"[:31]
        sheet = want if want in xl.sheet_names else next(
            s for s in xl.sheet_names if s.startswith("PerCpG_"))
        pc = xl.parse(sheet)
        bs = xl.parse("Block_Summary")
        if bs["cell_type_id"].eq(ct_id).any():
            bs = bs[bs["cell_type_id"] == ct_id]
        for row in bs.itertuples(index=False):
            chrom = str(row.Chr)
            if chrom not in ch_masks:
                continue
            sub = pc[pc["seq_id"] == row.seq_id]
            plus = sub[sub["strand"] == "+"]["cpg_position"].astype(int) - 1
            minus = sub[sub["strand"] == "-"]["cpg_position"].astype(int) - 1
            t_plus, t_minus = ch_masks[chrom]["T"]
            ch_masks[chrom]["T"] = (
                np.union1d(t_plus, plus.to_numpy(dtype=np.int64)),
                np.union1d(t_minus, minus.to_numpy(dtype=np.int64)))
        print(f"  Added block sites from {sheet} to target mask",
              file=sys.stderr)

    prepared = [PlantPrimer(prepare_primer(p["left_primer"],
                                           p["right_primer"],
                                           assay_id=p.get("assay_id", "?")))
                for p in primers_data]

    t0 = time.time()
    error_free = [pr for pr in prepared if pr.error is None]
    for chrom in chroms:
        # early-FAIL: a primer with >1 distinct loci already cannot PASS
        active = [pr for pr in error_free if len(pr.loci) <= 1]
        if not active:
            break
        genome = np.load(os.path.join(args.genome_dir, chrom + ".npy"))
        glen = len(genome)
        comp = np.where(genome < 4, 3 - genome, 4).astype(np.uint8)
        rc = comp[::-1].copy()
        del comp
        keep_c, keep_g = cpg_masks(genome)
        if use_block_mask:
            t_plus, t_minus = block_sites.get(
                chrom, (np.zeros(0, dtype=np.int64),
                        np.zeros(0, dtype=np.int64)))
        else:
            t_plus, t_minus = ch_masks[chrom]["T"]
        b_plus, b_minus = ch_masks[chrom]["B"]
        keep_c_T = keep_c.copy(); keep_c_T[t_plus] = True
        keep_g_T = keep_g.copy(); keep_g_T[t_minus] = True
        keep_c_B = keep_c.copy(); keep_c_B[b_plus] = True
        keep_g_B = keep_g.copy(); keep_g_B[b_minus] = True

        convs = {
            "T": (masked_ct(genome, keep_c_T), masked_ct(rc, keep_g_T[::-1])),
            "B": (masked_ct(genome, keep_c_B), masked_ct(rc, keep_g_B[::-1])),
            "U": (convert_genome_ct(genome), convert_genome_ct(rc)),
        }
        indexes = {k: (build_index(v[0]), build_index(v[1]))
                   for k, v in convs.items()}
        del genome, rc, keep_c, keep_g, keep_c_T, keep_g_T, keep_c_B, keep_g_B

        for pr in active:
            for state in ("T", "B", "U"):
                conv_w, conv_c = convs[state]
                idx_w, idx_c = indexes[state]
                for sgen, idx, tag in ((conv_w, idx_w, "W"),
                                       (conv_c, idx_c, "C")):
                    fp_pos, fp_mm = search_and_verify(
                        pr.info.fp_variants["CT_M"], idx, sgen,
                        pr.info.fp_arr, MAX_MM)
                    rp_pos, rp_mm = search_and_verify(
                        pr.info.rp_variants["CT_M"], idx, sgen,
                        pr.info.rp_rc_arr, MAX_MM)
                    for prod in find_products_vec(
                            fp_pos, fp_mm, rp_pos, rp_mm,
                            pr.info.fp_len, pr.info.rp_len,
                            MAX_PCR_LEN, chrom, state):
                        # report Crick products in forward genome space
                        prod.strand = tag
                        prod.chrom_len = glen
                        pr.add_product(state, prod)
        del convs, indexes

    n_pass = 0
    elapsed = round(time.time() - t0, 2)
    for p, pr in zip(primers_data, prepared):
        if pr.error:
            p["ePCR_verdict"] = "FAIL"
            p["ePCR_product_count"] = 0
            p["ePCR_products"] = ""
            p["ePCR_error"] = pr.error
        else:
            # locus verdict: exactly 1 amplified locus, and its best
            # target-state hit has strictly fewer mismatches than any
            # background/converted hit at that locus
            loci = list(pr.loci.values())
            n_t = pr.n_T
            ok = (len(loci) == 1 and loci[0]["t_mm"] is not None
                  and (loci[0]["bu_mm"] is None
                       or loci[0]["t_mm"] < loci[0]["bu_mm"]))
            if ok:
                p["ePCR_verdict"] = "PASS"
                p["ePCR_product_count"] = 1
                p["ePCR_products"] = loci[0]["prod"].coordinates()
                p["ePCR_error"] = ""
                n_pass += 1
            else:
                p["ePCR_verdict"] = "FAIL"
                p["ePCR_product_count"] = len(loci)
                t_prods = [d["prod"] for d in loci if d["prod"] is not None]
                p["ePCR_products"] = "; ".join(
                    x.coordinates() for x in t_prods[:5])
                p["ePCR_error"] = (f"target_products={n_t}, "
                                   f"background_products={pr.n_BU}")
        p["ePCR_elapsed_sec"] = elapsed

    with open(args.primers, "w") as f:
        json.dump(primers_data, f, indent=2)
    print(f"Plant ePCR done: {n_pass}/{len(prepared)} PASS; updated "
          f"{args.primers}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("design")
    d.add_argument("--dmr-xlsx", required=True)
    d.add_argument("--cell-type", required=True)
    d.add_argument("--fasta-dir", required=True)
    d.add_argument("--output-dir", required=True)
    d.add_argument("--flank", type=int, default=100)
    d.add_argument("--min-tm", type=float, default=58)
    d.add_argument("--opt-tm", type=float, default=60)
    d.add_argument("--max-tm", type=float, default=62)
    d.add_argument("--product-size-min", type=int, default=60)
    d.add_argument("--product-size-max", type=int, default=150)
    d.add_argument("--min-cpg-pair-total", type=int, default=3)
    d.add_argument("--min-cpg-per-primer", type=int, default=1)
    d.add_argument("--min-disc-in-primers", type=int, default=1,
                   help="minimum discriminating sites that must fall inside "
                        "primer spans (pairs with 0 cannot discriminate "
                        "methylation)")
    d.add_argument("--cpg-placement-mode", default="strict_terminal")

    e = sub.add_parser("epcr")
    e.add_argument("--primers", required=True)
    e.add_argument("--genome-dir", required=True)
    e.add_argument("--samples", required=True)
    e.add_argument("--cell-type", required=True,
                   help="Target GROUP name (as in samples.tsv)")
    e.add_argument("--meth-beta", type=float, default=0.65)
    e.add_argument("--min-cov", type=int, default=3)
    e.add_argument("--dmr-xlsx", default=None,
                   help="DMR block workbook; T-state = CpG + ALL block disc "
                        "sites (template-consistent mode)")
    e.add_argument("--cell-type-id", default=None,
                   help="cell_type_id (group_direction) whose blocks define "
                        "the T-state mask")

    args = ap.parse_args()
    if args.cmd == "design":
        run_design(args)
    else:
        run_epcr(args)


if __name__ == "__main__":
    main()
