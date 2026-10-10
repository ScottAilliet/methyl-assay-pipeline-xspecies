#!/usr/bin/env python3
"""
phase4_bisulfite_ePCR.py - BiSearch ePCR primer specificity screening for bisulfite PCR.

Recreates the BiSearch ePCR algorithm (Tusnady et al. 2005 NAR 33:e9;
Aranyi et al. 2006 BMC Bioinformatics 7:431) using 16-mer hash-based
indexing with sorted numpy arrays and vectorized verification.

Checks primer specificity against the bisulfite-converted genome on BOTH
strands and BOTH methylation states:
  - Watson (forward) strand: sense assays (SU/SM) amplify here.
  - Crick (reverse-complement) strand: antisense assays (AU/AM) amplify here.
  Each strand is converted two ways: CT (all C->T, unmethylated) and
  CT_M (CpG C preserved, non-CpG C->T, methylated).

The unconverted genome is NOT checked (bisulfite PCR templates are converted).

Algorithm per primer pair:
  1. Seed = 3'-most 16 bases of the FP and of RC(RP).
  2. Generate all seed variants within mismatch tolerance (<=3 mismatches).
  3. For each strand (Watson, Crick) and each state (CT, CT_M):
     look up variants in the sorted converted-genome index (np.searchsorted),
     verify full primer at ALL seed hits (vectorized numpy), and pair FP with
     RC(RP) hits (correctly oriented, product size <= 1000 bp).
  4. Deduplicate products per (chrom, strand, fp_pos, rp_pos) so a locus found
     in both methylation states counts once; sum Watson and Crick products.
  5. Verdict: PASS = exactly 1 product total genome-wide; FAIL = 0 or >1.

Usage:
  # Single pair
  python phase4_bisulfite_ePCR.py screen --fp ATCG... --rp GGTG... --genome-dir hg19_bs

  # Batch from CSV (columns: assay_id, fp_seq, rp_seq, cell_type, delta_means, ...)
  python phase4_bisulfite_ePCR.py screen --input primers.csv --output results.csv \
      --genome-dir hg19_bs --max-pass-per-ct 20 --screen-per-ct 200

  # Selection from screening results
  python phase4_bisulfite_ePCR.py select --input results.csv --output selected.csv \
      --max-per-celltype 20
"""

import argparse
import csv
import itertools
import json
import os
import sys
import time
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set

import numpy as np

# ============================================================================
# Constants
# ============================================================================

SEED_LEN = 16
MAX_MM = 3
MAX_PCR_LEN = 1000
MM_STRING = "1233333333333333"

BASE_CODES = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
CODE_TO_BASE = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

# Template configurations
# CT:   all C->T, alphabet {A,G,T} (3-letter)
# GA:   all G->A, alphabet {A,C,T} (3-letter)
# CT_M: CpG C preserved, non-CpG C->T, alphabet {A,C,G,T} (4-letter)
# GA_M: CpG G preserved, non-CpG G->A, alphabet {A,C,G,T} (4-letter)
TEMPLATES = {
    'CT':   {'alphabet': [0, 2, 3], 'methylated': False, 'state': 'unmethylated'},
    'GA':   {'alphabet': [0, 1, 3], 'methylated': False, 'state': 'unmethylated'},
    'CT_M': {'alphabet': [0, 1, 2, 3], 'methylated': True, 'state': 'methylated'},
    'GA_M': {'alphabet': [0, 1, 2, 3], 'methylated': True, 'state': 'methylated'},
}

# Both FP and RC(RP) are searched against the same converted template.
# The same seed variant set is reused for the Watson and Crick indexes of a
# given template, so only the CT and CT_M variant sets are precomputed.
FP_TEMPLATES = ['CT', 'CT_M']
RP_TEMPLATES = ['CT', 'CT_M']

# (methylation state, converted-template name) searched on each strand.
STRAND_STATES = [('unmethylated', 'CT'), ('methylated', 'CT_M')]

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


# ============================================================================
# Bisulfite conversion (string-based, matching pipeline/bisulfite.py)
# ============================================================================

def bisulfite_convert_forward(seq: str, methylated: bool = False) -> str:
    """C->T conversion on forward strand."""
    seq = seq.upper()
    if methylated:
        result = []
        for i in range(len(seq)):
            if seq[i] == 'C':
                if i + 1 < len(seq) and seq[i + 1] == 'G':
                    result.append('C')
                else:
                    result.append('T')
            else:
                result.append(seq[i])
        return ''.join(result)
    else:
        return seq.replace('C', 'T')


def bisulfite_convert_reverse(seq: str, methylated: bool = False) -> str:
    """G->A conversion on forward strand (reverse primer perspective)."""
    seq = seq.upper()
    if methylated:
        result = []
        for i in range(len(seq)):
            if seq[i] == 'G':
                if i > 0 and seq[i - 1] == 'C':
                    result.append('G')
                else:
                    result.append('A')
            else:
                result.append(seq[i])
        return ''.join(result)
    else:
        return seq.replace('G', 'A')


def reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))


# ============================================================================
# Sequence encoding
# ============================================================================

def encode_sequence(seq: str) -> Tuple[np.ndarray, bool]:
    """Encode a DNA sequence as numpy uint8 array. Returns (array, all_valid)."""
    arr = np.zeros(len(seq), dtype=np.uint8)
    valid = True
    for i, base in enumerate(seq.upper()):
        code = BASE_CODES.get(base, 4)
        arr[i] = code
        if code >= 4:
            valid = False
    return arr, valid


# ============================================================================
# K-mer encoding
# ============================================================================

def encode_kmer(seq16: str) -> Optional[int]:
    """
    Encode a 16-mer string as uint32.
    Position 0 (5'-most) = bits 30-31, position 15 (3'-most) = bits 0-1.
    """
    if len(seq16) != SEED_LEN:
        return None
    code = 0
    for j in range(SEED_LEN):
        base = BASE_CODES.get(seq16[j].upper(), 4)
        if base >= 4:
            return None
        code |= base << (2 * (15 - j))
    return code


# ============================================================================
# Genome conversion (numpy-based, for index building)
# ============================================================================

def convert_genome_ct(genome: np.ndarray) -> np.ndarray:
    """CT conversion: all C(1) -> T(3)."""
    converted = genome.copy()
    converted[converted == 1] = 3
    return converted


def convert_genome_ga(genome: np.ndarray) -> np.ndarray:
    """GA conversion: all G(2) -> A(0)."""
    converted = genome.copy()
    converted[converted == 2] = 0
    return converted


def convert_genome_ct_m(genome: np.ndarray) -> np.ndarray:
    """CT_M conversion: CpG C preserved, non-CpG C -> T."""
    cpg = np.zeros(len(genome), dtype=bool)
    cpg[:-1] = (genome[:-1] == 1) & (genome[1:] == 2)
    converted = genome.copy()
    is_c = (converted == 1)
    converted[is_c & ~cpg] = 3
    return converted


def convert_genome_ga_m(genome: np.ndarray) -> np.ndarray:
    """GA_M conversion: CpG G preserved, non-CpG G -> A."""
    cpg_g = np.zeros(len(genome), dtype=bool)
    cpg_g[1:] = (genome[:-1] == 1) & (genome[1:] == 2)
    converted = genome.copy()
    is_g = (converted == 2)
    converted[is_g & ~cpg_g] = 0
    return converted


CONVERTERS = {
    'CT': convert_genome_ct,
    'GA': convert_genome_ga,
    'CT_M': convert_genome_ct_m,
    'GA_M': convert_genome_ga_m,
}


# ============================================================================
# Index building
# ============================================================================

def build_index(genome_converted: np.ndarray) -> np.ndarray:
    """
    Build a sorted index of (kmer, position) pairs for a converted chromosome.
    Returns: sorted numpy uint64 array where each entry = (kmer << 32) | position.
    16-mers containing invalid bases (N, etc.) are excluded.
    """
    n = len(genome_converted) - SEED_LEN + 1
    if n <= 0:
        return np.array([], dtype=np.uint64)

    windows = np.lib.stride_tricks.sliding_window_view(genome_converted, SEED_LEN)
    valid_mask = (windows.max(axis=1) < 4)
    if not valid_mask.any():
        return np.array([], dtype=np.uint64)

    valid_windows = windows[valid_mask]
    positions = np.where(valid_mask)[0].astype(np.uint32)

    kmers = np.zeros(len(valid_windows), dtype=np.uint32)
    for j in range(SEED_LEN):
        shift = np.uint32(2 * (15 - j))
        kmers |= (valid_windows[:, j].astype(np.uint32) << shift)

    packed = (kmers.astype(np.uint64) << np.uint64(32)) | positions.astype(np.uint64)
    packed.sort()
    return packed


# ============================================================================
# Variant generation (numpy-optimized)
# ============================================================================

def _mm_string_valid(combo: tuple, mm_string: str) -> bool:
    """
    Check if a combination of mismatch positions satisfies the MM_STRING
    cumulative constraint (BiSearch standard).  Positions are 0-indexed from
    the 3' end of the seed.  For each position p in the sorted combo, the
    cumulative number of mismatches from the 3' end up to and including p
    must not exceed int(mm_string[p]).

    Example: MM_STRING = "1233333333333333"
      - 3'-most base (pos 0): at most 1 cumulative mismatch
      - first 2 bases (pos 0-1): at most 2 cumulative mismatches
      - first 3+ bases: at most 3 cumulative mismatches
    """
    for j, pos in enumerate(combo):
        if (j + 1) > int(mm_string[pos]):
            return False
    return True


def precompute_patterns(alphabet_codes: List[int], max_mm: int = MAX_MM,
                        seed_len: int = SEED_LEN,
                        mm_string: str = MM_STRING) -> Dict:
    """
    Pre-compute numpy arrays for fast variant generation.
    Returns dict keyed by k (1..max_mm) with arrays for vectorized substitution.
    Only combinations satisfying the MM_STRING cumulative constraint are included.
    """
    patterns = {}
    for k in range(1, max_mm + 1):
        clear_masks = []
        sub_bits = []
        shifts = []
        sub_vals = []

        for combo in itertools.combinations(range(seed_len), k):
            if not _mm_string_valid(combo, mm_string):
                continue
            cm = 0xFFFFFFFF
            for p in combo:
                cm &= ~(0b11 << (2 * p))

            for subs in itertools.product(alphabet_codes, repeat=k):
                sb = 0
                for i, p in enumerate(combo):
                    sb |= subs[i] << (2 * p)

                clear_masks.append(cm)
                sub_bits.append(sb)
                shifts.append([2 * p for p in combo])
                sub_vals.append(list(subs))

        patterns[k] = {
            'clear_masks': np.array(clear_masks, dtype=np.uint32),
            'sub_bits': np.array(sub_bits, dtype=np.uint32),
            'shifts': np.array(shifts, dtype=np.uint32),
            'sub_vals': np.array(sub_vals, dtype=np.uint32),
        }
    return patterns


def generate_variants_fast(seed_code: int, patterns: Dict,
                           max_mm: int = MAX_MM) -> np.ndarray:
    """
    Generate all variants of a 16-mer seed within mismatch tolerance.
    Returns: sorted numpy uint32 array of unique variant k-mer codes (including original).
    """
    seed = np.uint32(int(seed_code))
    all_variants = [int(seed_code)]

    for k in range(1, max_mm + 1):
        p = patterns[k]
        base = seed & p['clear_masks']
        candidates = base | p['sub_bits']

        # Filter: each substituted position must differ from original
        valid = np.ones(len(candidates), dtype=bool)
        for j in range(k):
            orig = (seed >> p['shifts'][:, j]) & np.uint32(0b11)
            valid &= (p['sub_vals'][:, j] != orig)

        all_variants.extend(candidates[valid].tolist())

    return np.unique(np.array(all_variants, dtype=np.uint32))


# Pattern cache (built once, reused for all primers)
_PATTERN_CACHE: Dict[Tuple, Dict] = {}


def get_patterns(alphabet_codes: List[int]) -> Dict:
    """Get or create precomputed patterns for an alphabet."""
    key = tuple(sorted(alphabet_codes))
    if key not in _PATTERN_CACHE:
        _PATTERN_CACHE[key] = precompute_patterns(alphabet_codes)
    return _PATTERN_CACHE[key]


# ============================================================================
# Seed lookup (returns numpy array of positions)
# ============================================================================

def search_seed_positions(variant_codes: np.ndarray,
                          index: np.ndarray) -> np.ndarray:
    """
    Look up variant k-mer codes in the sorted index.
    Returns: 1D numpy int64 array of hit positions on the genome.
    Uses vectorized range expansion (cumsum trick) to avoid Python loops.
    """
    if len(variant_codes) == 0 or len(index) == 0:
        return np.array([], dtype=np.int64)

    variants = np.unique(variant_codes)
    packed_low = variants.astype(np.uint64) << np.uint64(32)
    packed_high = packed_low | np.uint64(0xFFFFFFFF)

    left = np.searchsorted(index, packed_low, side='left')
    right = np.searchsorted(index, packed_high, side='right')

    counts = right - left
    has_hits = counts > 0

    if not has_hits.any():
        return np.array([], dtype=np.int64)

    # Vectorized range expansion: build flat indices into the sorted index
    # for all contiguous ranges [left[i], right[i]) in one pass.
    starts = left[has_hits].astype(np.int64)
    sizes = counts[has_hits].astype(np.int64)
    total = int(sizes.sum())

    flat_idx = np.ones(total, dtype=np.int64)
    flat_idx[0] = starts[0]
    # At each range boundary, jump to the next start
    boundaries = np.cumsum(sizes)[:-1]
    flat_idx[boundaries] = starts[1:] - starts[:-1] - sizes[:-1] + 1
    flat_idx = np.cumsum(flat_idx)

    all_hits = index[flat_idx]
    positions = (all_hits & np.uint64(0xFFFFFFFF)).astype(np.int64)

    return positions


# ============================================================================
# Vectorized full-primer verification
# ============================================================================

def verify_hits_vectorized(positions: np.ndarray, converted_genome: np.ndarray,
                           primer_arr: np.ndarray,
                           max_mm: int,
                           chunk_size: int = 500_000) -> Tuple[np.ndarray, np.ndarray]:
    """
    Verify full primer at multiple seed hit positions using vectorized numpy.

    positions: 1D int64 array of seed hit positions (leftmost of 16-mer on genome)
    converted_genome: converted genome as uint8 array
    primer_arr: primer as uint8 array (same orientation as genome)
    max_mm: maximum allowed mismatches
    chunk_size: process at most this many hits per batch to bound peak memory
                (results are identical to processing all at once).

    Returns: (passed_positions, passed_mismatches) as numpy arrays.
    """
    if len(positions) == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int32)

    primer_len = len(primer_arr)
    offset = primer_len - SEED_LEN  # bases upstream of seed (5' extension)
    genome_len = len(converted_genome)

    # Filter positions where full primer fits within genome bounds
    valid = (positions >= offset) & (positions + SEED_LEN <= genome_len)
    valid_positions = positions[valid]

    if len(valid_positions) == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int32)

    col_offsets = np.arange(primer_len, dtype=np.int64)

    def _verify_chunk(chunk_positions):
        # genome_idx[i, j] = chunk_positions[i] - offset + j
        genome_idx = chunk_positions[:, None] - offset + col_offsets[None, :]
        genome_slices = converted_genome[genome_idx]  # (N, primer_len)
        mm = np.count_nonzero(genome_slices != primer_arr, axis=1)
        keep = mm <= max_mm
        return chunk_positions[keep], mm[keep].astype(np.int32)

    if len(valid_positions) <= chunk_size:
        return _verify_chunk(valid_positions)

    # Chunked path: bound the transient (N, primer_len) array for repetitive seeds.
    out_pos, out_mm = [], []
    for start in range(0, len(valid_positions), chunk_size):
        cp, cm = _verify_chunk(valid_positions[start:start + chunk_size])
        if len(cp):
            out_pos.append(cp)
            out_mm.append(cm)
    if out_pos:
        return np.concatenate(out_pos), np.concatenate(out_mm)
    return np.array([], dtype=np.int64), np.array([], dtype=np.int32)


# ============================================================================
# Combined search + verify
# ============================================================================

def search_and_verify(variants: np.ndarray, index: np.ndarray,
                      converted_genome: np.ndarray, primer_arr: np.ndarray,
                      max_mm: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Search seed variants in index and verify full primer at all hits.
    Returns: (passed_positions, passed_mismatches) as numpy arrays.
    """
    positions = search_seed_positions(variants, index)
    return verify_hits_vectorized(positions, converted_genome, primer_arr, max_mm)


# ============================================================================
# Product detection
# ============================================================================

@dataclass
class PCRProduct:
    """A predicted PCR product."""
    chrom: str
    state: str          # 'unmethylated' or 'methylated'
    fp_pos: int         # 0-based leftmost of FP seed on genome
    rp_pos: int         # 0-based leftmost of RP seed on genome
    fp_len: int
    rp_len: int
    size: int           # product size in bp
    fp_mm: int          # FP full-primer mismatches
    rp_mm: int          # RP full-primer mismatches
    total_mm: int       # fp_mm + rp_mm
    strand: str = 'W'   # 'W' (Watson/forward) or 'C' (Crick/reverse-complement)
    chrom_len: int = 0  # chromosome length (needed to map Crick coords to forward)

    def coordinates(self) -> str:
        """Human-readable product coordinates in forward (Watson) genome space."""
        if self.total_mm == 99:
            # Sentinel from the binary-search path: too many pairs to enumerate.
            return f"{self.chrom}:many-products(>1,{self.state},{self.strand})"
        if self.strand == 'C':
            # fp_pos/rp_pos are on the reverse-complement genome; map to forward.
            rc_fp_start = self.fp_pos - self.fp_len + SEED_LEN
            rc_rp_end = self.rp_pos + SEED_LEN - 1
            start = self.chrom_len - 1 - rc_rp_end
            end = self.chrom_len - 1 - rc_fp_start
        else:
            start = self.fp_pos - self.fp_len + SEED_LEN
            end = self.rp_pos + SEED_LEN - 1
        return (f"{self.chrom}:{start+1}-{end+1}"
                f"({self.size}bp,{self.state},{self.strand},mm={self.total_mm})")


def find_products_vec(fp_positions: np.ndarray, fp_mismatches: np.ndarray,
                      rp_positions: np.ndarray, rp_mismatches: np.ndarray,
                      fp_len: int, rp_len: int, max_pcr_len: int,
                      chrom: str, state: str,
                      max_pairs: int = 100_000) -> List[PCRProduct]:
    """
    Find PCR products by pairing FP and RP verified hits.
    Vectorized size computation; product objects created only for valid pairs.
    If potential pairs exceed max_pairs, uses sorted binary search instead of 2D grid.
    """
    if len(fp_positions) == 0 or len(rp_positions) == 0:
        return []

    n_fp = len(fp_positions)
    n_rp = len(rp_positions)
    potential = n_fp * n_rp

    if potential > max_pairs:
        # Too many potential pairs - use sorted binary search to count products
        # without building a 2D grid. Only creates objects if count is small.
        rp_sorted = np.sort(rp_positions)
        total_valid = 0
        for fp_pos in fp_positions:
            # size = rp_pos - fp_pos + fp_len > 0  =>  rp_pos > fp_pos - fp_len
            # size <= max_pcr_len                  =>  rp_pos <= fp_pos + max_pcr_len - fp_len
            lo = fp_pos - fp_len + 1
            hi = fp_pos + max_pcr_len - fp_len
            n = int(np.searchsorted(rp_sorted, hi, side='right') -
                    np.searchsorted(rp_sorted, lo, side='left'))
            total_valid += n
            if total_valid > 1:
                break

        if total_valid > 1:
            # Return 2 dummy products to trigger FAIL
            return [
                PCRProduct(chrom, state, 0, 0, fp_len, rp_len, 999, 99, 99, 99),
                PCRProduct(chrom, state, 1, 1, fp_len, rp_len, 999, 99, 99, 99),
            ]
        elif total_valid == 0:
            return []
        # If exactly 1, fall through to 2D grid to find the actual product

    # 2D grid approach for manageable hit counts
    sizes = rp_positions[:, None] - fp_positions[None, :] + fp_len
    valid = (sizes > 0) & (sizes <= max_pcr_len)
    valid_pairs = np.argwhere(valid)  # (rp_idx, fp_idx) pairs

    products = []
    for rp_idx, fp_idx in valid_pairs:
        products.append(PCRProduct(
            chrom=chrom, state=state,
            fp_pos=int(fp_positions[fp_idx]),
            rp_pos=int(rp_positions[rp_idx]),
            fp_len=fp_len, rp_len=rp_len,
            size=int(sizes[rp_idx, fp_idx]),
            fp_mm=int(fp_mismatches[fp_idx]),
            rp_mm=int(rp_mismatches[rp_idx]),
            total_mm=int(fp_mismatches[fp_idx]) + int(rp_mismatches[rp_idx]),
        ))
    return products


# ============================================================================
# Primer preparation
# ============================================================================

@dataclass
class PrimerInfo:
    """Prepared primer pair for screening."""
    assay_id: str
    fp_seq: str
    rp_seq: str
    fp_arr: np.ndarray
    rp_rc_arr: np.ndarray
    fp_len: int
    rp_len: int
    fp_variants: Dict[str, np.ndarray]  # template_name -> FP seed variant codes
    rp_variants: Dict[str, np.ndarray]  # template_name -> RC(RP) seed variant codes
    products: List[PCRProduct] = field(default_factory=list)
    seen_coords: Dict[Tuple, int] = field(default_factory=dict)  # dedup key -> index in products
    active: bool = True
    error: Optional[str] = None
    row: Optional[dict] = None  # original CSV row (for batch mode)


def prepare_primer(fp_seq: str, rp_seq: str, assay_id: str = "unknown",
                   max_mm: int = MAX_MM,
                   row: Optional[dict] = None) -> PrimerInfo:
    """
    Prepare a primer pair for screening: encode, extract seeds, generate variants.

    FP and RC(RP) seed variants are stored in separate dicts because both are
    searched against the same CT/CT_M templates.
    """
    fp_upper = fp_seq.upper()
    rp_upper = rp_seq.upper()

    for base in fp_upper + rp_upper:
        if base not in BASE_CODES:
            return PrimerInfo(assay_id, fp_seq, rp_seq,
                              np.array([]), np.array([]), 0, 0, {}, {},
                              error=f"Invalid base '{base}' in primer", row=row)

    fp_arr, _ = encode_sequence(fp_upper)
    rp_arr, _ = encode_sequence(rp_upper)
    fp_len = len(fp_arr)
    rp_len = len(rp_arr)

    if fp_len < SEED_LEN or rp_len < SEED_LEN:
        return PrimerInfo(assay_id, fp_seq, rp_seq,
                          fp_arr, np.array([]), fp_len, rp_len, {}, {},
                          error=f"Primer too short (min {SEED_LEN} bp)", row=row)

    # Extract seeds (3'-most 16 bases)
    fp_seed = fp_upper[-SEED_LEN:]
    rp_rc = reverse_complement(rp_upper)
    rp_seed = rp_rc[-SEED_LEN:]

    fp_seed_code = encode_kmer(fp_seed)
    rp_seed_code = encode_kmer(rp_seed)

    if fp_seed_code is None or rp_seed_code is None:
        return PrimerInfo(assay_id, fp_seq, rp_seq,
                          fp_arr, np.array([]), fp_len, rp_len, {}, {},
                          error="Invalid base in seed region", row=row)

    rp_rc_arr, _ = encode_sequence(rp_rc)

    # Generate FP and RC(RP) variants separately for each template
    fp_variants = {}
    for tname in FP_TEMPLATES:
        alphabet = TEMPLATES[tname]['alphabet']
        patterns = get_patterns(alphabet)
        fp_variants[tname] = generate_variants_fast(fp_seed_code, patterns, max_mm)
    rp_variants = {}
    for tname in RP_TEMPLATES:
        alphabet = TEMPLATES[tname]['alphabet']
        patterns = get_patterns(alphabet)
        rp_variants[tname] = generate_variants_fast(rp_seed_code, patterns, max_mm)

    return PrimerInfo(assay_id, fp_seq, rp_seq, fp_arr, rp_rc_arr,
                      fp_len, rp_len, fp_variants, rp_variants, row=row)


# ============================================================================
# Chromosome screening (core engine, shared by single and batch modes)
# ============================================================================

def screen_chromosome(primers: List[PrimerInfo], genome: np.ndarray,
                      chrom: str, max_mm: int, max_pcr_len: int,
                      verbose: bool = False) -> None:
    """
    Screen all active primers against one chromosome, checking BOTH strands.

    Watson (forward) strand carries sense products; Crick (reverse-complement)
    strand carries antisense products. For each strand the CT (unmethylated) and
    CT_M (methylated) converted indexes are built (2 indexes at a time to bound
    peak memory). Both the FP and RC(RP) seeds are searched against the same
    converted template, then paired into products.

    Products are deduplicated per (chrom, strand, fp_pos, rp_pos): a locus found
    in both methylation states counts once, while Watson and Crick products are
    summed. A primer with >1 total products genome-wide is marked inactive
    (early FAIL) and skipped on subsequent strands/chromosomes.
    """
    active = [p for p in primers if p.active and p.error is None]
    if not active:
        return

    if verbose:
        print(f"  {chrom} ({len(genome):,} bp, {len(active)} active primers)...",
              file=sys.stderr)

    glen = len(genome)

    # Build the Crick (reverse-complement) genome once (N/invalid -> 4 preserved).
    comp = np.where(genome < 4, 3 - genome, 4).astype(np.uint8)
    rc_genome = comp[::-1].copy()
    del comp

    for strand, sgen in (('W', genome), ('C', rc_genome)):
        # Stop early if nothing is still active (e.g. all FAILed on Watson).
        if not any(p.active and p.error is None for p in active):
            break

        # Build CT and CT_M indexes for this strand only (peak = 2 indexes).
        indexes: Dict[str, np.ndarray] = {}
        converted: Dict[str, np.ndarray] = {}
        for tname in ('CT', 'CT_M'):
            conv = CONVERTERS[tname](sgen)
            converted[tname] = conv
            indexes[tname] = build_index(conv)

        for p in active:
            if not p.active or p.error is not None:
                continue

            for state, tname in STRAND_STATES:
                fp_pos, fp_mm = search_and_verify(
                    p.fp_variants[tname], indexes[tname],
                    converted[tname], p.fp_arr, max_mm)
                rp_pos, rp_mm = search_and_verify(
                    p.rp_variants[tname], indexes[tname],
                    converted[tname], p.rp_rc_arr, max_mm)

                products = find_products_vec(
                    fp_pos, fp_mm, rp_pos, rp_mm,
                    p.fp_len, p.rp_len, max_pcr_len, chrom, state)

                # Dedup across strands/states; Watson + Crick are summed.
                # A locus seen in both methylation states counts once and is
                # reported in whichever state matched with fewer mismatches.
                for prod in products:
                    prod.strand = strand
                    prod.chrom_len = glen
                    key = (chrom, strand, prod.fp_pos, prod.rp_pos)
                    if key not in p.seen_coords:
                        p.seen_coords[key] = len(p.products)
                        p.products.append(prod)
                    else:
                        idx = p.seen_coords[key]
                        if prod.total_mm < p.products[idx].total_mm:
                            p.products[idx] = prod

            # Early FAIL: >1 products anywhere means non-specific.
            if len(p.products) > 1:
                p.active = False

        del indexes, converted

    del rc_genome


# ============================================================================
# Screening results
# ============================================================================

@dataclass
class ScreenResult:
    """Result of screening a single primer pair."""
    assay_id: str
    fp_seq: str
    rp_seq: str
    verdict: str       # 'PASS' or 'FAIL'
    product_count: int
    products: List[PCRProduct] = field(default_factory=list)
    error: Optional[str] = None
    elapsed_sec: float = 0.0

    def to_dict(self) -> dict:
        return {
            'assay_id': self.assay_id,
            'fp_seq': self.fp_seq,
            'rp_seq': self.rp_seq,
            'verdict': self.verdict,
            'product_count': self.product_count,
            'products': '; '.join(p.coordinates() for p in self.products),
            'error': self.error or '',
            'elapsed_sec': round(self.elapsed_sec, 2),
        }


def primer_to_result(p: PrimerInfo, elapsed: float) -> ScreenResult:
    """Convert a PrimerInfo (after screening) to a ScreenResult."""
    if p.error:
        return ScreenResult(p.assay_id, p.fp_seq, p.rp_seq, "FAIL", 0,
                            error=p.error, elapsed_sec=elapsed)
    verdict = "PASS" if len(p.products) == 1 else "FAIL"
    return ScreenResult(p.assay_id, p.fp_seq, p.rp_seq, verdict,
                        len(p.products), p.products, elapsed_sec=elapsed)


# ============================================================================
# Single-pair screening
# ============================================================================

def screen_pair(fp_seq: str, rp_seq: str, genome_dir: str,
                assay_id: str = "unknown",
                max_mm: int = MAX_MM,
                max_pcr_len: int = MAX_PCR_LEN,
                chromosomes: Optional[List[str]] = None,
                verbose: bool = False) -> ScreenResult:
    """Screen a single primer pair against the bisulfite-converted genome."""
    t0 = time.time()

    primer = prepare_primer(fp_seq, rp_seq, assay_id=assay_id, max_mm=max_mm)
    if primer.error:
        return ScreenResult(assay_id, fp_seq, rp_seq, "FAIL", 0,
                            error=primer.error, elapsed_sec=time.time() - t0)

    chromosomes = _resolve_chromosomes(genome_dir, chromosomes)

    for chrom in chromosomes:
        chrom_path = os.path.join(genome_dir, f"{chrom}.npy")
        if not os.path.exists(chrom_path):
            continue
        genome = np.load(chrom_path)
        screen_chromosome([primer], genome, chrom, max_mm, max_pcr_len, verbose)
        del genome

        if not primer.active:
            break

    return primer_to_result(primer, time.time() - t0)


# ============================================================================
# Batch screening with DMS-ranked early stopping
# ============================================================================

def screen_batch(input_csv: str, output_csv: str, genome_dir: str,
                 max_mm: int = MAX_MM, max_pcr_len: int = MAX_PCR_LEN,
                 chromosomes: Optional[List[str]] = None,
                 max_pass_per_ct: int = 20,
                 screen_per_ct: int = 200,
                 verbose: bool = False) -> None:
    """
    Screen primers from CSV in DMS-ranked order with early stopping.

    1. Read all primers, group by cell_type, sort by delta_means (descending).
    2. Take top screen_per_ct primers per cell type.
    3. Screen across all chromosomes (build indexes once per chromosome).
    4. Select top max_pass_per_ct passing assays per cell type (one per DMR).
    5. Write results to output CSV.

    Input CSV columns: assay_id, fp_seq (or left_primer), rp_seq (or right_primer),
                       cell_type, delta_means, chrom, dmr_start, dmr_end, ...
    """
    t_start = time.time()

    # Read input
    with open(input_csv, newline='') as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        print("No primers found in input CSV.", file=sys.stderr)
        return

    print(f"Read {len(rows)} primers from {input_csv}", file=sys.stderr)

    # Filter to target-detecting primers before triage.
    # Background-detecting primers (M primers on hypo DMRs, U primers on hyper
    # DMRs) are never selected for the panel, so screening them wastes ePCR
    # capacity and displaces target-detecting primers from the top-N triage.
    rows = [r for r in rows if _is_target_detecting(r)]
    print(f"  After target-detecting filter: {len(rows)} primers", file=sys.stderr)

    # Group by cell type and sort by DMS (descending)
    by_ct: Dict[str, List[dict]] = {}
    for r in rows:
        ct = r.get('cell_type', 'unknown')
        by_ct.setdefault(ct, []).append(r)

    def get_dms(row):
        try:
            val = row.get('dms_score')
            if val is not None and val != '':
                return float(val)
            return float(row.get('delta_means', 0) or 0)
        except (ValueError, TypeError):
            return 0.0

    for ct in by_ct:
        by_ct[ct].sort(key=get_dms, reverse=True)

    # Take top screen_per_ct per cell type
    to_screen = []
    ct_counts = {}
    for ct, assays in sorted(by_ct.items()):
        batch = assays[:screen_per_ct]
        to_screen.extend(batch)
        ct_counts[ct] = len(batch)

    total_to_screen = len(to_screen)
    print(f"Screening {total_to_screen} primers "
          f"(top {screen_per_ct} per cell type, {len(ct_counts)} cell types)",
          file=sys.stderr)

    # Prepare all primers (encode, generate variants)
    t_prep = time.time()
    primers: List[PrimerInfo] = []
    n_error = 0
    for row in to_screen:
        assay_id = row.get('assay_id', f'row_{len(primers)}')
        fp_seq = row.get('fp_seq', row.get('left_primer', ''))
        rp_seq = row.get('rp_seq', row.get('right_primer', ''))
        p = prepare_primer(fp_seq, rp_seq, assay_id=assay_id, max_mm=max_mm, row=row)
        if p.error:
            n_error += 1
        primers.append(p)

    print(f"Prepared {len(primers)} primers in {time.time()-t_prep:.1f}s "
          f"({n_error} errors)", file=sys.stderr)

    # Resolve chromosomes
    chromosomes = _resolve_chromosomes(genome_dir, chromosomes)

    # Screen chromosome by chromosome
    for chrom in chromosomes:
        chrom_path = os.path.join(genome_dir, f"{chrom}.npy")
        if not os.path.exists(chrom_path):
            continue

        t_chrom = time.time()
        genome = np.load(chrom_path)
        screen_chromosome(primers, genome, chrom, max_mm, max_pcr_len, verbose)
        del genome

        active = sum(1 for p in primers if p.active and p.error is None)
        n_pass = sum(1 for p in primers if not p.error and len(p.products) == 1)
        if verbose:
            print(f"  {chrom} done in {time.time()-t_chrom:.1f}s | "
                  f"active={active} pass={n_pass}", file=sys.stderr)

    # Determine verdicts
    n_pass = 0
    n_fail = 0
    for p in primers:
        if p.error:
            n_fail += 1
        elif len(p.products) == 1:
            n_pass += 1
        else:
            n_fail += 1

    total_time = time.time() - t_start
    print(f"\nScreening complete: {len(primers)} primers in {total_time:.1f}s",
          file=sys.stderr)
    print(f"  PASS: {n_pass}  FAIL: {n_fail}  ERROR: {n_error}", file=sys.stderr)

    # Write all results
    _write_results(primers, output_csv)

    # Selection: top max_pass_per_ct per cell type
    selected_csv = output_csv.replace('.csv', '_selected.csv')
    select_from_primers(primers, selected_csv, max_per_celltype=max_pass_per_ct)

    print(f"  Results: {output_csv}", file=sys.stderr)
    print(f"  Selected: {selected_csv}", file=sys.stderr)


def _write_results(primers: List[PrimerInfo], output_csv: str) -> None:
    """Write screening results for all primers to CSV."""
    if not primers:
        return

    # Collect all column names from the first primer's row
    input_cols = []
    if primers[0].row:
        input_cols = list(primers[0].row.keys())
    output_cols = input_cols + ['verdict', 'product_count', 'products', 'error', 'elapsed_sec']

    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=output_cols)
        writer.writeheader()
        for p in primers:
            row = dict(p.row) if p.row else {}
            if p.error:
                row['verdict'] = 'FAIL'
                row['product_count'] = 0
                row['products'] = ''
                row['error'] = p.error
            else:
                row['verdict'] = 'PASS' if len(p.products) == 1 else 'FAIL'
                row['product_count'] = len(p.products)
                row['products'] = '; '.join(prod.coordinates() for prod in p.products)
                row['error'] = ''
            row['elapsed_sec'] = ''
            writer.writerow(row)


# ============================================================================
# Selection logic
# ============================================================================

def _is_target_detecting(row: dict) -> bool:
    """
    Determine whether a primer detects the target cell type's methylation state.

    For hypo DMRs (target unmethylated, background methylated):
      - U primers (SU, AU) amplify unmethylated DNA → detect the target.
      - M primers (SM, AM) amplify methylated DNA → detect the background.

    For hyper DMRs (target methylated, background unmethylated):
      - M primers (SM, AM) amplify methylated DNA → detect the target.
      - U primers (SU, AU) amplify unmethylated DNA → detect the background.

    Only target-detecting primers are retained for panel selection.
    """
    try:
        target = float(row.get('target_mean_meth', 0) or 0)
        background = float(row.get('background_mean_meth', 0) or 0)
    except (ValueError, TypeError):
        return False
    template = row.get('template_used', '')
    is_methylated_primer = template in ('SM', 'AM')
    if target < background:  # hypo: target unmethylated
        return not is_methylated_primer  # U primer detects target
    else:                      # hyper: target methylated
        return is_methylated_primer   # M primer detects target


def select_from_primers(primers: List[PrimerInfo], output_csv: str,
                        max_per_celltype: int = 20,
                        one_per_dmr: bool = True) -> None:
    """
    Select passing assays ranked by DMS (delta_means).
    Keeps one best passing assay per DMR, up to max_per_celltype DMRs per cell type.

    Only target-detecting primers are considered: for hypo DMRs (target
    unmethylated), only U primers (SU/AU) are eligible; for hyper DMRs
    (target methylated), only M primers (SM/AM) are eligible.  Background-
    detecting primers are excluded so that every panel assay uniquely
    amplifies the target cell type's DNA after bisulfite conversion.
    """
    passing = [p for p in primers
               if not p.error and len(p.products) == 1
               and _is_target_detecting(p.row or {})]

    def get_dms(p):
        if p.row:
            try:
                val = p.row.get('dms_score')
                if val is not None and val != '':
                    return float(val)
                return float(p.row.get('delta_means', 0) or 0)
            except (ValueError, TypeError):
                return 0.0
        return 0.0

    # Group by cell type
    by_ct: Dict[str, List[PrimerInfo]] = {}
    for p in passing:
        ct = (p.row or {}).get('cell_type', 'unknown')
        by_ct.setdefault(ct, []).append(p)

    selected = []
    for ct, assays in sorted(by_ct.items()):
        assays.sort(key=get_dms, reverse=True)
        count = 0
        if one_per_dmr:
            seen_dmr = set()
            for a in assays:
                dmr_key = (
                    (a.row or {}).get('chrom', ''),
                    (a.row or {}).get('dmr_start', ''),
                    (a.row or {}).get('dmr_end', ''),
                )
                if dmr_key not in seen_dmr:
                    seen_dmr.add(dmr_key)
                    selected.append(a)
                    count += 1
                    if count >= max_per_celltype:
                        break
        else:
            for a in assays[:max_per_celltype]:
                selected.append(a)

    # Write output
    if selected:
        input_cols = list((selected[0].row or {}).keys())
        output_cols = input_cols + ['verdict', 'product_count', 'products', 'error']
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=output_cols)
            writer.writeheader()
            for p in selected:
                row = dict(p.row) if p.row else {}
                row['verdict'] = 'PASS'
                row['product_count'] = len(p.products)
                row['products'] = '; '.join(prod.coordinates() for prod in p.products)
                row['error'] = ''
                writer.writerow(row)

    # Summary
    by_ct_summary = {}
    for p in selected:
        ct = (p.row or {}).get('cell_type', 'unknown')
        by_ct_summary[ct] = by_ct_summary.get(ct, 0) + 1

    print(f"Selection complete: {len(selected)} assays selected", file=sys.stderr)
    for ct, n in sorted(by_ct_summary.items()):
        print(f"  {ct}: {n}", file=sys.stderr)


def select_assays(input_csv: str, output_csv: str,
                  max_per_celltype: int = 20,
                  one_per_dmr: bool = True) -> None:
    """
    Select passing assays from a screening results CSV.
    Ranks by delta_means (descending) within each cell type.

    Only target-detecting primers are considered (see _is_target_detecting).
    """
    with open(input_csv, newline='') as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    passing = [r for r in rows
               if r.get('verdict') == 'PASS' and _is_target_detecting(r)]

    def get_dms(row):
        try:
            val = row.get('dms_score')
            if val is not None and val != '':
                return float(val)
            return float(row.get('delta_means', 0) or 0)
        except (ValueError, TypeError):
            return 0.0

    by_ct: Dict[str, List[dict]] = {}
    for r in passing:
        ct = r.get('cell_type', 'unknown')
        by_ct.setdefault(ct, []).append(r)

    selected = []
    for ct, assays in sorted(by_ct.items()):
        assays.sort(key=get_dms, reverse=True)
        count = 0
        if one_per_dmr:
            seen_dmr = set()
            for a in assays:
                dmr_key = (a.get('chrom', ''), a.get('dmr_start', ''), a.get('dmr_end', ''))
                if dmr_key not in seen_dmr:
                    seen_dmr.add(dmr_key)
                    selected.append(a)
                    count += 1
                    if count >= max_per_celltype:
                        break
        else:
            selected.extend(assays[:max_per_celltype])

    if selected:
        fieldnames = list(selected[0].keys())
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(selected)

    by_ct_summary = {}
    for r in selected:
        ct = r.get('cell_type', 'unknown')
        by_ct_summary[ct] = by_ct_summary.get(ct, 0) + 1

    print(f"Selection complete: {len(selected)} assays selected", file=sys.stderr)
    for ct, n in sorted(by_ct_summary.items()):
        print(f"  {ct}: {n}", file=sys.stderr)
    print(f"  Output: {output_csv}", file=sys.stderr)


# ============================================================================
# Helpers
# ============================================================================

def _resolve_chromosomes(genome_dir: str,
                         chromosomes: Optional[List[str]] = None) -> List[str]:
    """Determine which chromosomes to screen.

    Always scans the directory for .npy files as the source of truth.
    Uses chrom_info.json only for ordering if it exists.
    """
    if chromosomes is not None:
        return chromosomes

    found = []
    for fname in sorted(os.listdir(genome_dir)):
        if fname.startswith('chr') and fname.endswith('.npy'):
            found.append(fname[:-4])

    chrom_info_path = os.path.join(genome_dir, "chrom_info.json")
    if os.path.exists(chrom_info_path):
        try:
            with open(chrom_info_path) as f:
                chrom_info = json.load(f)
            ordered = [c for c in chrom_info['chromosomes'] if c in found]
            extras = [c for c in found if c not in ordered]
            return ordered + extras
        except (KeyError, json.JSONDecodeError):
            pass

    return found


# ============================================================================
# Pipeline integration interface
# ============================================================================

def run_phase4_epcr(
    pairs,
    config,
    top_n: int = 5,
    genome_dir: Optional[str] = None,
    threads: int = 4,
    fallback_to_web: bool = True,
    max_mm: Optional[int] = None,
    max_pcr_len: Optional[int] = None,
):
    """
    Pipeline Step 4 entry point: hash-based ePCR specificity check.

    Uses the 16-mer hash algorithm (BiSearch ePCR) against bisulfite-converted
    genome arrays. Maintains the same interface for pipeline.py.

    Parameters
    ----------
    pairs        : List[PrimerPair] from Phase 3 (passing pairs first)
    config       : PipelineConfig
    top_n        : Only check the top N CpG-passing pairs
    genome_dir   : Directory with converted genome .npy files
    threads      : CPU threads (unused in hash-based version, kept for compat)
    fallback_to_web : Kept for interface compat (hash-based is always local)
    max_mm       : Override config bisearch_max_mismatches
    max_pcr_len  : Override config bisearch_pcr_length

    Returns the same list with bisearch_* fields populated.
    """
    # Resolve parameters from config
    if max_mm is None:
        max_mm = getattr(config, 'bisearch_max_mismatches', MAX_MM)
    if max_pcr_len is None:
        max_pcr_len = getattr(config, 'bisearch_pcr_length', MAX_PCR_LEN)
    if genome_dir is None:
        genome_dir = getattr(config, 'bisearch_genome_dir', None)

    if genome_dir is None or not os.path.isdir(genome_dir):
        msg = (f"Phase 4: Genome directory not found: {genome_dir}\n"
               "  Run prepare_genomes.py to build converted genome files.")
        if fallback_to_web:
            print(f"{msg}\n  Falling back to web BiSearch.", file=sys.stderr)
            return _fallback_web(pairs, config, top_n)
        else:
            print(msg, file=sys.stderr)
            return pairs

    passing = [p for p in pairs if p.cpg_pass]
    to_check = passing[:top_n]

    if not to_check:
        print("Phase 4: No CpG-passing pairs to check.", file=sys.stderr)
        return pairs

    print(f"Phase 4: Hash-based ePCR - checking top {len(to_check)} pairs",
          file=sys.stderr)

    # Prepare all primers
    primers = []
    for pair in to_check:
        p = prepare_primer(pair.fwd_seq, pair.rev_seq,
                           assay_id=pair.assay_id, max_mm=max_mm)
        p._pair = pair  # link back to PrimerPair
        primers.append(p)

    # Resolve chromosomes
    chromosomes = _resolve_chromosomes(genome_dir)

    # Screen chromosome by chromosome
    for chrom in chromosomes:
        chrom_path = os.path.join(genome_dir, f"{chrom}.npy")
        if not os.path.exists(chrom_path):
            continue
        genome = np.load(chrom_path)
        screen_chromosome(primers, genome, chrom, max_mm, max_pcr_len)
        del genome

    # Update PrimerPair fields
    for p in primers:
        pair = getattr(p, '_pair', None)
        if pair is None:
            continue
        if p.error:
            pair.bisearch_amplicons = 0
            pair.bisearch_pass = False
            pair.bisearch_verdict = f"ERROR: {p.error}"
        else:
            pair.bisearch_amplicons = len(p.products)
            pair.bisearch_pass = (len(p.products) == 1)
            if len(p.products) == 1:
                pair.bisearch_verdict = "SPECIFIC"
            elif len(p.products) == 0:
                pair.bisearch_verdict = "NO_PRODUCT"
            else:
                pair.bisearch_verdict = "NON-SPECIFIC"

        # Methylated state fields
        m_products = [prod for prod in p.products if prod.state == 'methylated']
        um_products = [prod for prod in p.products if prod.state == 'unmethylated']
        pair.m_bisearch_amplicons = len(m_products)
        pair.m_bisearch_pass = (len(m_products) <= 1)
        pair.m_bisearch_verdict = (
            "SPECIFIC" if len(m_products) == 1 else
            "NO_PRODUCT" if len(m_products) == 0 else "NON-SPECIFIC"
        )
        pair.m_fp_hits = sum(1 for prod in p.products if prod.state == 'methylated')
        pair.m_rp_hits = pair.m_fp_hits

        print(f"  {pair.assay_id}: {pair.bisearch_verdict} "
              f"({pair.bisearch_amplicons} products)", file=sys.stderr)

    # Re-rank: BiSearch-passing pairs first
    def sort_key(p):
        bs_rank = 0 if p.bisearch_pass else (1 if p.bisearch_pass is None else 2)
        cpg_rank = 0 if p.cpg_pass else 1
        return (bs_rank, cpg_rank, -p.pair_cpg_total)

    pairs.sort(key=sort_key)
    return pairs


def _fallback_web(pairs, config, top_n):
    """Fallback to web BiSearch (imported lazily to avoid circular deps)."""
    try:
        from .phase4_bisearch import run_phase4 as run_web
        return run_web(pairs, config, manual_mode=False)
    except Exception:
        print("Phase 4: Web fallback also unavailable. Skipping.", file=sys.stderr)
        return pairs


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="BiSearch ePCR: bisulfite primer specificity screening"
    )
    subparsers = parser.add_subparsers(dest='command', help='Subcommand')

    # Screen subcommand
    screen_parser = subparsers.add_parser('screen', help='Screen primers for specificity')
    screen_group = screen_parser.add_mutually_exclusive_group(required=True)
    screen_group.add_argument('--fp', help='Forward primer sequence (5\'->3\')')
    screen_group.add_argument('--input', help='Input CSV file (batch mode)')
    screen_parser.add_argument('--rp', help='Reverse primer sequence (5\'->3\')')
    screen_parser.add_argument('--output', help='Output CSV file (batch mode)')
    screen_parser.add_argument('--genome-dir', required=True,
                               help='Directory containing converted genome .npy files')
    screen_parser.add_argument('--assay-id', default='unknown', help='Assay identifier')
    screen_parser.add_argument('--max-mm', type=int, default=MAX_MM,
                               help=f'Max mismatches per primer (default: {MAX_MM})')
    screen_parser.add_argument('--max-pcr-len', type=int, default=MAX_PCR_LEN,
                               help=f'Max PCR product size in bp (default: {MAX_PCR_LEN})')
    screen_parser.add_argument('--chromosomes', nargs='*',
                               help='Chromosomes to screen (default: all in genome dir)')
    screen_parser.add_argument('--max-pass-per-ct', type=int, default=20,
                               help='Max passing assays to select per cell type (default: 20)')
    screen_parser.add_argument('--screen-per-ct', type=int, default=200,
                               help='Max primers to screen per cell type (default: 200)')
    screen_parser.add_argument('--verbose', '-v', action='store_true',
                               help='Verbose progress output')

    # Select subcommand
    select_parser = subparsers.add_parser('select', help='Select passing assays by DMS rank')
    select_parser.add_argument('--input', required=True, help='Screening results CSV')
    select_parser.add_argument('--output', required=True, help='Output selected panel CSV')
    select_parser.add_argument('--max-per-celltype', type=int, default=20,
                               help='Max assays per cell type (default: 20)')
    select_parser.add_argument('--all-per-dmr', action='store_true',
                               help='Allow multiple assays per DMR (default: one per DMR)')

    args = parser.parse_args()

    if args.command == 'screen':
        if args.fp:
            # Single pair mode
            if not args.rp:
                parser.error("--rp is required when using --fp")
            result = screen_pair(
                args.fp, args.rp, args.genome_dir,
                assay_id=args.assay_id,
                max_mm=args.max_mm,
                max_pcr_len=args.max_pcr_len,
                chromosomes=args.chromosomes,
                verbose=args.verbose,
            )
            print(f"Assay:      {result.assay_id}")
            print(f"FP:         {result.fp_seq}")
            print(f"RP:         {result.rp_seq}")
            print(f"Verdict:    {result.verdict}")
            print(f"Products:   {result.product_count}")
            if result.products:
                for p in result.products:
                    print(f"  {p.coordinates()}")
            if result.error:
                print(f"Error:      {result.error}")
            print(f"Time:       {result.elapsed_sec:.2f}s")

        elif args.input:
            # Batch mode
            if not args.output:
                parser.error("--output is required when using --input")
            screen_batch(
                args.input, args.output, args.genome_dir,
                max_mm=args.max_mm,
                max_pcr_len=args.max_pcr_len,
                chromosomes=args.chromosomes,
                max_pass_per_ct=args.max_pass_per_ct,
                screen_per_ct=args.screen_per_ct,
                verbose=args.verbose,
            )

    elif args.command == 'select':
        select_assays(
            args.input, args.output,
            max_per_celltype=args.max_per_celltype,
            one_per_dmr=not args.all_per_dmr,
        )

    else:
        parser.print_help()


if __name__ == '__main__':
    main()
