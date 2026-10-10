#!/usr/bin/env python3
"""
phase0_generic.py - Organism-generic DMR discovery for cross-species
cell-type-specific MSP assay design.

Reconstructed from the session fix log after a sandbox loss; validated
end-to-end by reproducing the arabidopsis v8 discovery results
(Columella_hyper: 200 blocks / 1083 PerCpG sites / 8 primer pairs / 2
ePCR-PASS with identical coordinates).

Loaders (per sample) return dict[key] -> (pos, meth, tot[, ctx][, strand]):
  mouse      load_mouse_cov        Bismark-style .cov.gz
                                    (chrom, start 0-based, end, pct, meth, unmeth)
                                    key = chrom
  arabidopsis load_arabidopsis      Kawakatsu mc_calls TSV
                                    (assembly, position, strand, class, mc, h)
                                    key = "chrom:context"
  zebrafish  load_zebrafish_bigwig score-only per-CpG bigWig (no coverage)
                                    key = chrom

Position convention for arabidopsis mc_calls (validated against TAIR10 base
identity on 70,272 chr5 CHH sites: top[rec]='C' for 100% of '+' rows and
top[rec]='G' for 100% of '-' rows): position is the 0-based plus-strand
coordinate of the RECORDED base -- for '+' rows the cytosine itself, for
'-' rows the paired top-strand base (the motif's last base: G at cytosine+1
for CG, last base at cytosine+2 for CHG/CHH). The decode is therefore the
IDENTITY map: returned pos = recorded position, with a parallel strand
array. Downstream, a '-' site at plus-coord x pairs with the bottom-strand
C at len-1-x.

Methylation state model per sample matrix:
  arabidopsis  files are pre-filtered to methylated calls (mc >= 2), so
               unrecorded sites are effectively unmethylated:
               absent_beta=0.0, absent_cov=nominal coverage (30x).
  mouse        cov files record both methylated and unmethylated calls;
               absent sites have no information: absent_beta=NaN.
  zebrafish    score-only bigWigs carry no coverage track (tot=0 sentinel);
               beta = meth/100 and covm = NaN, so coverage gating falls back
               to "observed in >= half the targets".

samples.tsv columns: sample_id, path, group, role
  One comparison block per target group: rows of that group with
  role=target are the target replicates; rows with role=background (any
  group) are the background replicates. To discover markers for several
  groups from the same files, repeat each file once per comparison with
  comparison-unique sample_ids.

Usage:
  python -m methyl_panel.phase0_generic --organism arabidopsis \
      --samples samples.tsv --out dmr_blocks.xlsx \
      --contexts CG,CHG,CHH --chrom-prefix ""
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Sample loaders: dict[key] -> per-chrom numpy arrays (pos, meth, total[, ...])
# pos is 0-based genomic coordinate of the cytosine (see module docstring for
# the arabidopsis recorded-base convention).
# ---------------------------------------------------------------------------

def load_mouse_cov(path, min_cov=1):
    """Bismark-style cov: chrom, start(0-based), end, pct, meth, unmeth."""
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end", "pct", "meth", "unmeth"],
                     usecols=["chrom", "start", "meth", "unmeth"],
                     dtype={"chrom": "category", "start": "int64",
                            "meth": "int32", "unmeth": "int32"})
    out = {}
    for chrom, g in df.groupby("chrom", observed=True, sort=False):
        pos = g["start"].values
        meth = g["meth"].values
        tot = (g["meth"].values + g["unmeth"].values)
        keep = tot >= min_cov
        pos, meth, tot = pos[keep], meth[keep], tot[keep]
        order = np.argsort(pos, kind="stable")
        out[str(chrom)] = (pos[order], meth[order].astype(np.int32),
                           tot[order].astype(np.int32))
    return out


def load_arabidopsis(path, min_cov=1, chrom_prefix=""):
    """Kawakatsu mc_calls TSV: assembly, position, strand, class, mc, h.

    Position convention (validated against TAIR10 base identity): position
    is the 0-based coordinate of the RECORDED base -- for '+' rows the
    cytosine itself; for '-' rows the paired top-strand base (the LAST base
    of the motif: CG: the G at cytosine+1; CHG/CHH: the last base at
    cytosine+2). The decode is the IDENTITY map (cyto = rec): returned pos
    is the recorded plus-strand coordinate, with a parallel strand array.
    Arrays keyed 'chrom:context'.

    h is TOTAL coverage (mc <= h always); beta = mc/h.
    NOTE: these files are pre-filtered to methylated calls (mc >= 2);
    unrecorded sites are effectively unmethylated.
    """
    df = pd.read_csv(path, sep="\t",
                     usecols=["assembly", "position", "strand", "class", "mc", "h"],
                     dtype={"assembly": "category", "position": "int64",
                            "strand": "category", "class": "category",
                            "mc": "int32", "h": "int32"})
    out = {}
    for (chrom, ctx), g in df.groupby(["assembly", "class"], observed=True,
                                      sort=False):
        rec = g["position"].values  # recorded-base plus-coord (identity decode)
        strand = g["strand"].values.astype(str)
        meth = g["mc"].values
        tot = g["h"].values
        keep = tot >= min_cov
        rec, meth, tot, strand = (rec[keep], meth[keep], tot[keep],
                                  strand[keep])
        order = np.argsort(rec, kind="stable")
        out[chrom_prefix + str(chrom) + ":" + str(ctx)] = (
            rec[order], meth[order].astype(np.int32),
            tot[order].astype(np.int32), str(ctx), strand[order])
    return out


def load_zebrafish_bigwig(path, min_beta_valid=True):
    """Per-CpG methylation scores from bigWig; no coverage track."""
    import pyBigWig
    bw = pyBigWig.open(str(path))
    out = {}
    for c, length in bw.chroms().items():
        if not c.startswith("chr") or c == "chrM":
            continue
        vals = bw.values(c, 0, length, numpy=True)
        idx = np.nonzero(~np.isnan(vals))[0]
        if len(idx) == 0:
            continue
        beta = vals[idx].astype(np.float32)
        # sentinel: total=0 marks "no coverage track"
        out[c] = (idx.astype(np.int64),
                  np.round(beta * 100).astype(np.int32),
                  np.zeros(len(idx), dtype=np.int32))
    bw.close()
    return out


# ---------------------------------------------------------------------------
# Per-chromosome beta matrices
# ---------------------------------------------------------------------------

def chrom_beta(pos_by_sample, key, absent_beta=np.nan, absent_cov=np.nan):
    """Build per-sample beta/coverage matrices over the union of positions.

    pos_by_sample: sample_id -> loader dict; key selects the chromosome
    (mouse/zebrafish) or 'chrom:context' (arabidopsis).

    Zebrafish sentinel (tot=0, score-only): beta = meth/100, covm = NaN.
    Arabidopsis passes absent_beta=0.0 and absent_cov=nominal coverage so
    unrecorded (pre-filtered) sites count as unmethylated with coverage.
    """
    mats = []
    for sid, d in pos_by_sample.items():
        if key in d:
            pos, meth, tot = d[key][:3]
            strand = d[key][4] if len(d[key]) >= 5 else None
            mats.append((sid, pos, meth, tot, strand))
    if not mats:
        return (np.zeros(0, dtype=np.int64), np.zeros((0, 0), dtype=np.float32),
                np.zeros((0, 0), dtype=np.float32), None)
    union = np.unique(np.concatenate([m[1] for m in mats]))
    beta = np.full((len(mats), len(union)), absent_beta, dtype=np.float32)
    covm = np.full((len(mats), len(union)), absent_cov, dtype=np.float32)
    strand_union = np.full(len(union), '.', dtype='U1')
    have_strand = False
    for i, (_sid, pos, meth, tot, strand) in enumerate(mats):
        idx = np.searchsorted(union, pos)
        has_cov = tot > 0
        b = np.where(has_cov, meth / np.maximum(tot, 1), meth / 100.0)
        beta[i, idx] = b
        covm[i, idx] = np.where(has_cov, tot, np.nan)
        if strand is not None:
            have_strand = True
            strand = np.asarray(strand)
            unfilled = strand_union[idx] == '.'
            strand_union[idx[unfilled]] = strand[unfilled]
    return union, beta, covm, (strand_union if have_strand else None)


# ---------------------------------------------------------------------------
# Block-level discrimination score
# ---------------------------------------------------------------------------

def dms_score(tg_mat, bg_mat, direction):
    """Discrimination score in [0, 1] for one DMR block.

    A/B: state separation (target unmethylated + background methylated for
    hypo, and vice versa for hyper); C/D: per-sample consistency within
    target and background; E: fraction of sites observed in >= half the
    target replicates.
    """
    tg_mean = np.nanmean(tg_mat, axis=0)
    bg_mean = np.nanmean(bg_mat, axis=0)
    if tg_mat.size == 0 or bg_mat.size == 0 or np.all(np.isnan(tg_mean)) \
            or np.all(np.isnan(bg_mean)):
        return 0.0
    if direction == "hypo":
        A = max(0.0, float(1.0 - np.nanmax(tg_mean)))
        B = max(0.0, float(np.nanmin(bg_mean)))
    else:
        A = max(0.0, float(np.nanmin(tg_mean)))
        B = max(0.0, float(1.0 - np.nanmax(bg_mean)))
    C = max(0.0, 1.0 - 2.0 * float(np.nanmean(np.nanstd(tg_mat, axis=0))))
    D = max(0.0, 1.0 - 2.0 * float(np.nanmean(np.nanstd(bg_mat, axis=0))))
    E = float(np.mean(np.mean(~np.isnan(tg_mat), axis=0) >= 0.5))
    return float(np.mean([A, B, C, D]) * E)


# ---------------------------------------------------------------------------
# Per-chromosome discovery
# ---------------------------------------------------------------------------

def discover_chrom(union_pos, beta, covm, target_idx, background_idx,
                   cell_type_id, chrom, direction="hypo",
                   min_cpg=4, max_gap=200, delta_site=0.5,
                   delta_block=0.3, unmeth_thresh=0.15,
                   meth_thresh=0.65, min_cov=3, dms_min=0.6,
                   context=None, strand_arr=None):
    """Find DMR blocks on one chromosome/context.

    Candidates: coverage-OK (>= half the targets >= min_cov, or observed in
    >= half the targets when no coverage track exists) AND in the target
    methylation state AND per-site delta >= delta_site. Candidates are
    clustered by max_gap; clusters with >= min_cpg sites and block mean
    delta >= delta_block are emitted with a dms_score >= dms_min.

    strand_arr: parallel array for the arabidopsis recorded-base strands
    ('+'/'-'); must be indexed on the SAME union as beta (the compacted
    candidate columns, not the pre-compaction union).
    """
    tg = beta[target_idx]
    bg = beta[background_idx]
    tg_cov = covm[target_idx]
    has_cov = not np.all(np.isnan(covm))

    n_need = max(1, len(target_idx) // 2)
    if has_cov:
        cov_ok = np.sum(tg_cov >= min_cov, axis=0) >= n_need
    else:
        cov_ok = np.sum(~np.isnan(tg), axis=0) >= n_need

    with np.errstate(invalid="ignore"):
        tg_mean = np.nanmean(tg, axis=0)
        bg_mean = np.nanmean(bg, axis=0)
        delta = bg_mean - tg_mean if direction == "hypo" else tg_mean - bg_mean
        if direction == "hypo":
            state_ok = (tg_mean <= unmeth_thresh) & (bg_mean >= meth_thresh)
        else:
            state_ok = (tg_mean >= meth_thresh) & (bg_mean <= unmeth_thresh)
    cand = cov_ok & state_ok & (delta >= delta_site)
    if not cand.any():
        return [], []
    cand_pos = union_pos[cand]
    cand_delta = delta[cand]
    cand_strand = strand_arr[cand] if strand_arr is not None else None
    if len(cand_pos) < min_cpg:
        return [], []

    # cluster candidates by max_gap
    clusters, cur = [], [0]
    for i in range(1, len(cand_pos)):
        if cand_pos[i] - cand_pos[i - 1] <= max_gap:
            cur.append(i)
        else:
            clusters.append(cur)
            cur = [i]
    clusters.append(cur)

    rows, percpg_rows = [], []
    bid = 0
    for cl in clusters:
        cols = np.array(cl)
        if len(cols) < min_cpg:
            continue
        block_delta = float(np.nanmean(cand_delta[cols]))
        if block_delta < delta_block:
            continue
        tg_mat = tg[:, cols]
        bg_mat = bg[:, cols]
        score = dms_score(tg_mat, bg_mat, direction)
        if score < dms_min:
            continue
        with np.errstate(invalid="ignore"):
            tgm = float(np.nanmean(tg_mean[cols]))
            bgm = float(np.nanmean(bg_mean[cols]))
            delta_means = float(np.nanmean(delta[cols]))
            tg_block = np.nanmean(tg_mat, axis=1)
            bg_block = np.nanmean(bg_mat, axis=1)
        vt = tg_block[~np.isnan(tg_block)]
        vb = bg_block[~np.isnan(bg_block)]
        if len(vt) >= 2 and len(vb) >= 2:
            pval = float(stats.ttest_ind(vt, vb, equal_var=False).pvalue)
        else:
            pval = None
        start, end = int(cand_pos[cl[0]]), int(cand_pos[cl[-1]])
        # include chrom+context so IDs are unique genome-wide
        seq_id = f"{cell_type_id}_{chrom}_{context or 'CG'}_{bid:04d}_{direction[0].upper()}"
        bid += 1
        rows.append({
            "seq_id": seq_id, "cell_type_id": cell_type_id,
            "Chr": chrom, "Start": start + 1, "End": end + 1,
            "Direction": direction, "N_Sites": int(len(cols)),
            "Target_Mean_Meth": float(tgm),
            "Background_Mean_Meth": float(bgm),
            "Delta_Means": float(delta_means),
            "Ttest_Pval": "" if pval is None or np.isnan(pval) else f"{pval:.3g}",
            "dms_score": score, "Target_CellType": cell_type_id,
            "context": context or "CG",
        })
        with np.errstate(invalid="ignore"):
            tg_cpg_mean = np.nanmean(tg_mat, axis=0)
            bg_cpg_mean = np.nanmean(bg_mat, axis=0)
        for j, col in enumerate(cols):
            percpg_rows.append({
                "seq_id": seq_id, "cpg_label": f"CpG_{j+1}",
                "cpg_position": int(cand_pos[col]) + 1,
                "cpg_global_idx": int(cand_pos[col]) + 1,
                "target_mean_beta": float(tg_cpg_mean[j]),
                "background_mean_beta": float(bg_cpg_mean[j]),
                "delta_beta": abs(float(tg_cpg_mean[j] - bg_cpg_mean[j])),
                "context": context or "CG",
                "strand": str(cand_strand[col]) if cand_strand is not None else ".",
            })
    return rows, percpg_rows


def write_block_summary_xlsx(blocks_by_ct, percpg_by_ct, out_path):
    if not blocks_by_ct:
        raise SystemExit("No DMR blocks found for any cell type — "
                         "check thresholds and input data")
    bs_all = pd.concat(blocks_by_ct, ignore_index=True)
    bs_all = bs_all.sort_values(["cell_type_id", "dms_score"],
                                ascending=[True, False])
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as xw:
        bs_all.to_excel(xw, sheet_name="Block_Summary", index=False)
        for ct, pc in percpg_by_ct.items():
            if len(pc):
                sheet = f"PerCpG_{ct}"[:31]
                pc.to_excel(xw, sheet_name=sheet, index=False)
    return out_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", required=True,
                    choices=["mouse", "arabidopsis", "zebrafish"])
    ap.add_argument("--samples", required=True,
                    help="TSV: sample_id, path, group, role (target/background)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--direction", default="both",
                    choices=["both", "hypo", "hyper"])
    ap.add_argument("--contexts", default="CG",
                    help="comma-separated contexts (arabidopsis: CG,CHG,CHH)")
    ap.add_argument("--min-cpg", type=int, default=4)
    ap.add_argument("--max-gap", type=int, default=200)
    ap.add_argument("--delta-site", type=float, default=0.5)
    ap.add_argument("--delta-block", type=float, default=0.3)
    ap.add_argument("--min-cov", type=int, default=3)
    ap.add_argument("--top-n", type=int, default=200)
    ap.add_argument("--dms-min", type=float, default=0.6)
    ap.add_argument("--chrom-prefix", default="",
                    help="Prefix to add to chromosome names (e.g. 'chr')")
    args = ap.parse_args()

    s = pd.read_csv(args.samples, sep="\t")
    contexts = args.contexts.split(",")

    print("Loading samples ...", file=sys.stderr)
    loaders = {
        "mouse": (load_mouse_cov, {}),
        "arabidopsis": (load_arabidopsis, {"chrom_prefix": args.chrom_prefix}),
        "zebrafish": (load_zebrafish_bigwig, {}),
    }
    loader, kwargs = loaders[args.organism]
    pos_by_sample = {}
    for rec in s.itertuples(index=False):
        if rec.sample_id in pos_by_sample:
            continue
        print(f"  {rec.sample_id}: {rec.path}", file=sys.stderr)
        pos_by_sample[rec.sample_id] = loader(rec.path, **kwargs)

    # chromosome universe (arabidopsis keys are 'chrom:context')
    chroms = set()
    for d in pos_by_sample.values():
        for key in d:
            chroms.add(key.split(":")[0])
    chroms = sorted(chroms)

    # arabidopsis absent sites are effectively unmethylated (mc >= 2 pre-filter)
    absent_beta = 0.0 if args.organism == "arabidopsis" else np.nan
    absent_cov = 30.0 if args.organism == "arabidopsis" else np.nan

    directions = (["hypo", "hyper"] if args.direction == "both"
                  else [args.direction])

    blocks_by_ct, percpg_by_ct = [], {}
    for ct, g in s.groupby("group"):
        targets = g[g["role"] == "target"]["sample_id"].tolist()
        backgrounds = s[s["role"] == "background"]["sample_id"].tolist()
        if not backgrounds:
            backgrounds = [sid for sid in pos_by_sample
                           if sid not in set(targets)]
        if not targets or not backgrounds:
            print(f"  Skipping group {ct}: missing targets/backgrounds",
                  file=sys.stderr)
            continue
        target_idx = [list(pos_by_sample).index(sid) for sid in targets]
        background_idx = [list(pos_by_sample).index(sid) for sid in backgrounds]
        for direction in directions:
            cell_type_id = f"{ct}_{direction}"
            ct_rows, ct_percpg = [], []
            for chrom in chroms:
                for ctx in contexts:
                    if args.organism == "arabidopsis":
                        key = f"{chrom}:{ctx}"
                        if not any(key in d for d in pos_by_sample.values()):
                            continue
                        context = ctx
                    else:
                        key = chrom
                        if not any(key in d for d in pos_by_sample.values()):
                            continue
                        context = None
                    union_pos, beta, covm, strand_union = chrom_beta(
                        pos_by_sample, key, absent_beta=absent_beta,
                        absent_cov=absent_cov)
                    if len(union_pos) == 0:
                        continue
                    rows, percpg = discover_chrom(
                        union_pos, beta, covm, target_idx, background_idx,
                        cell_type_id, chrom, direction=direction,
                        min_cpg=args.min_cpg, max_gap=args.max_gap,
                        delta_site=args.delta_site,
                        delta_block=args.delta_block,
                        unmeth_thresh=0.15, meth_thresh=0.65,
                        min_cov=args.min_cov, dms_min=args.dms_min,
                        context=context, strand_arr=strand_union)
                    ct_rows.extend(rows)
                    ct_percpg.extend(percpg)
            # cap per cell type + direction by dms_score
            if ct_rows and args.top_n and len(ct_rows) > args.top_n:
                ct_rows.sort(key=lambda r: -r["dms_score"])
                keep = {r["seq_id"] for r in ct_rows[:args.top_n]}
                ct_rows = [r for r in ct_rows if r["seq_id"] in keep]
                ct_percpg = [p for p in ct_percpg
                             if p["seq_id"] in keep]
            print(f"  {cell_type_id}: {len(ct_rows)} blocks", file=sys.stderr)
            if ct_rows:
                blocks_by_ct.append(pd.DataFrame(ct_rows))
                percpg_by_ct[cell_type_id] = pd.DataFrame(ct_percpg)

    write_block_summary_xlsx(blocks_by_ct, percpg_by_ct, args.out)
    print(f"Wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
