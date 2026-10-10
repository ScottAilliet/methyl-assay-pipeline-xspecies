#!/usr/bin/env python3
"""
validate_epcr.py - Brute-force validation of the hash-based ePCR search engine.

This is an engine-level (machinery) test. For a set of test primers it compares
the hash-based seed search + vectorized verify against an exhaustive brute-force
scan of every position on a converted chromosome. The two must return exactly
the same set of hit positions.

The search machinery is strand-agnostic: the Watson strand and the Crick
(reverse-complement) strand are both screened by converting the sequence with
CT (all C->T, unmethylated) and CT_M (CpG C preserved, non-CpG C->T, methylated)
and running this same seed-lookup + verify path. Validating the engine on CT and
CT_M therefore validates the machinery used on both strands. Product pairing and
strand-coordinate mapping are covered separately by the end-to-end tests.

Usage:
  python validate_epcr.py --genome-dir hg19_bs_test --chrom chr22
"""

import argparse
import sys
import time
import numpy as np

from methyl_panel.phase4_bisulfite_ePCR import (
    SEED_LEN, MAX_MM, CONVERTERS, FP_TEMPLATES, RP_TEMPLATES,
    reverse_complement, build_index, search_seed_positions,
    verify_hits_vectorized, prepare_primer,
)

DECODE = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

# Templates actually used by the screening engine (both strands, both states).
ENGINE_TEMPLATES = sorted(set(FP_TEMPLATES) | set(RP_TEMPLATES))  # ['CT', 'CT_M']


def brute_force_scan(primer_arr: np.ndarray, converted_genome: np.ndarray,
                     max_mm: int) -> np.ndarray:
    """
    Scan every position on the converted genome for full-primer matches.
    Returns array of full-primer start positions where mismatches <= max_mm.
    """
    primer_len = len(primer_arr)
    genome_len = len(converted_genome)
    n = genome_len - primer_len + 1

    if n <= 0:
        return np.array([], dtype=np.int64)

    # Sliding window view: (n, primer_len)
    windows = np.lib.stride_tricks.sliding_window_view(converted_genome, primer_len)

    # Filter out windows containing invalid bases (code >= 4, e.g. N)
    valid_mask = windows.max(axis=1) < 4
    valid_windows = windows[valid_mask]
    valid_positions = np.where(valid_mask)[0]

    # Count mismatches per window
    mismatches = np.count_nonzero(valid_windows != primer_arr, axis=1)

    # Return positions with mismatches <= max_mm
    passed = mismatches <= max_mm
    return valid_positions[passed]


def validate_chromosome(genome_dir: str, chrom: str, max_mm: int = MAX_MM,
                         n_test_primers: int = 20, verbose: bool = True):
    """
    Validate the hash-based ePCR engine against brute force on one chromosome.

    For each test primer, and for the FP and RC(RP) arrays against each engine
    template (CT, CT_M):
      1. Run brute-force scan (every position).
      2. Run hash-based search_seed_positions + verify_hits_vectorized.
      3. The two hit sets must be identical.
    """
    genome_path = f"{genome_dir}/{chrom}.npy"
    genome = np.load(genome_path)
    print(f"Loaded {chrom}: {len(genome):,} bp")

    # Build converted genomes and indexes (only the templates the engine uses).
    indexes = {}
    converted = {}
    for tname in ENGINE_TEMPLATES:
        conv = CONVERTERS[tname](genome)
        converted[tname] = conv
        t0 = time.time()
        indexes[tname] = build_index(conv)
        if verbose:
            print(f"  {tname} index: {len(indexes[tname]):,} entries ({time.time()-t0:.1f}s)")

    # Generate test primers from random genome positions. Both the FP and RP are
    # drawn from the CT-converted genome: the engine searches the FP and RC(RP)
    # against the same CT / CT_M templates on each strand, so RC(RP) must match a
    # CT-converted substring (not a GA one).
    np.random.seed(42)
    valid_mask = genome < 4
    test_primers = []
    for _ in range(n_test_primers * 3):  # try more, keep good ones
        if len(test_primers) >= n_test_primers:
            break
        start = np.random.randint(1000, len(genome) - 300)
        # Require both the FP window and the downstream RP window to be N-free.
        if not valid_mask[start:start + 225].all():
            continue
        # Forward primer: a 25-mer from the CT-converted genome.
        fp_seq = ''.join(DECODE[int(c)] for c in converted['CT'][start:start + 25])
        # Reverse primer: RC of a downstream 25-mer from the CT-converted genome,
        # so that RC(RP) == that CT substring and matches the CT index.
        rp_ct = ''.join(DECODE[int(c)] for c in converted['CT'][start + 200:start + 225])
        rp_seq = reverse_complement(rp_ct)
        test_primers.append((fp_seq, rp_seq, start))

    print(f"\nValidating {len(test_primers)} test primers on {chrom} "
          f"(templates: {', '.join(ENGINE_TEMPLATES)})...")

    all_pass = True
    n_checks = 0
    for i, (fp_seq, rp_seq, _expected_pos) in enumerate(test_primers):
        p = prepare_primer(fp_seq, rp_seq, assay_id=f'val_{i}', max_mm=max_mm)
        if p.error:
            print(f"  Primer {i}: SKIP (error: {p.error})")
            continue

        # Test the FP array and the RC(RP) array, each against its templates.
        roles = [
            ('FP', p.fp_arr, p.fp_variants, FP_TEMPLATES),
            ('RP', p.rp_rc_arr, p.rp_variants, RP_TEMPLATES),
        ]
        for role, primer_arr, variants_map, templates in roles:
            for tname in templates:
                # Brute-force: scan every position.
                t0 = time.time()
                bf_positions = brute_force_scan(primer_arr, converted[tname], max_mm)
                bf_time = time.time() - t0

                # Hash-based: seed search + vectorized verify.
                t0 = time.time()
                hash_positions_raw = search_seed_positions(
                    variants_map[tname], indexes[tname])
                hash_positions, _hash_mm = verify_hits_vectorized(
                    hash_positions_raw, converted[tname], primer_arr, max_mm)
                hash_time = time.time() - t0

                # Convert seed positions (3'-most 16-mer) to full-primer starts.
                offset = len(primer_arr) - SEED_LEN
                hash_full_start = hash_positions - offset

                bf_sorted = np.sort(bf_positions)
                hash_sorted = np.sort(hash_full_start)
                n_checks += 1

                if np.array_equal(bf_sorted, hash_sorted):
                    if verbose:
                        print(f"  Primer {i} {role}/{tname}: PASS "
                              f"({len(bf_sorted)} hits, "
                              f"bf={bf_time:.2f}s hash={hash_time:.3f}s)")
                else:
                    all_pass = False
                    bf_set = set(bf_sorted.tolist())
                    hash_set = set(hash_sorted.tolist())
                    only_bf = sorted(bf_set - hash_set)
                    only_hash = sorted(hash_set - bf_set)
                    print(f"  Primer {i} {role}/{tname}: MISMATCH!")
                    print(f"    Brute-force: {len(bf_sorted)} hits")
                    print(f"    Hash-based:  {len(hash_sorted)} hits")
                    if only_bf:
                        print(f"    Only in brute-force: {only_bf[:10]}")
                    if only_hash:
                        print(f"    Only in hash-based:  {only_hash[:10]}")

    print()
    if all_pass:
        print(f"VALIDATION PASSED: {n_checks} engine checks match brute-force on {chrom}")
    else:
        print("VALIDATION FAILED: Discrepancies found!")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Brute-force validation of the hash-based ePCR search engine"
    )
    parser.add_argument('--genome-dir', required=True,
                        help='Directory containing converted genome .npy files')
    parser.add_argument('--chrom', default='chr22',
                        help='Chromosome to validate (default: chr22)')
    parser.add_argument('--max-mm', type=int, default=MAX_MM,
                        help=f'Max mismatches (default: {MAX_MM})')
    parser.add_argument('--n-primers', type=int, default=20,
                        help='Number of test primers (default: 20)')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Less verbose output')

    args = parser.parse_args()
    validate_chromosome(
        args.genome_dir, args.chrom, args.max_mm,
        args.n_primers, verbose=not args.quiet
    )


if __name__ == '__main__':
    main()
