#!/usr/bin/env python3
"""
prepare_genomes.py - Download and prepare hg19 genome for ePCR screening.

Downloads per-chromosome FASTA files from UCSC, encodes each chromosome
as a numpy uint8 array, and saves to the genome directory.

The screening tool (phase4_bisulfite_ePCR.py) loads these .npy files and performs
bisulfite conversion + index building on-the-fly during screening.

Usage:
  # Download all standard chromosomes
  python prepare_genomes.py --output-dir hg19_bs

  # Download specific chromosomes only (for testing)
  python prepare_genomes.py --output-dir hg19_bs --chromosomes chr22

  # Use existing FASTA files instead of downloading
  python prepare_genomes.py --output-dir hg19_bs --fasta-dir /path/to/fasta
"""

import argparse
import gzip
import json
import os
import sys
import time
import urllib.request

import numpy as np

# Base encoding: A=0, C=1, G=2, T=3, everything else=4
BASE_CODES = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

UCSC_BASE_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes"


def download_chromosome(chrom: str, output_path: str) -> bool:
    """Download a single chromosome FASTA from UCSC and decompress."""
    url = f"{UCSC_BASE_URL}/{chrom}.fa.gz"
    print(f"  Downloading {chrom} from {url}...", file=sys.stderr)
    try:
        urllib.request.urlretrieve(url, output_path + ".gz")
        # Decompress
        with gzip.open(output_path + ".gz", 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                f_out.write(f_in.read())
        os.unlink(output_path + ".gz")
        size_mb = os.path.getsize(output_path) / 1e6
        print(f"  {chrom}: {size_mb:.1f} MB", file=sys.stderr)
        return True
    except Exception as e:
        print(f"  ERROR downloading {chrom}: {e}", file=sys.stderr)
        return False


def read_fasta(path: str) -> tuple:
    """Read a FASTA file. Returns (header, sequence_string)."""
    header = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    break
                header = line[1:]
            else:
                seq_parts.append(line)
    return header, ''.join(seq_parts)


def encode_genome(seq: str) -> np.ndarray:
    """Encode a DNA sequence string as numpy uint8 array."""
    arr = np.zeros(len(seq), dtype=np.uint8)
    for i, base in enumerate(seq.upper()):
        arr[i] = BASE_CODES.get(base, 4)
    return arr


def prepare_genome(output_dir: str,
                   chromosomes: list = None,
                   fasta_dir: str = None) -> None:
    """
    Prepare genome .npy files for ePCR screening.

    Downloads (or reads existing) per-chromosome FASTA files, encodes
    each as a numpy uint8 array, and saves to output_dir.
    Also writes chrom_info.json with chromosome names and lengths.
    """
    if chromosomes is None:
        chromosomes = CHROMOSOMES

    os.makedirs(output_dir, exist_ok=True)

    chrom_info = {'chromosomes': [], 'lengths': {}}
    tmp_dir = os.path.join(output_dir, 'tmp_fasta')
    os.makedirs(tmp_dir, exist_ok=True)

    for chrom in chromosomes:
        npy_path = os.path.join(output_dir, f"{chrom}.npy")

        # Skip if already prepared
        if os.path.exists(npy_path):
            arr = np.load(npy_path)
            print(f"  {chrom}: already prepared ({len(arr):,} bp)", file=sys.stderr)
            chrom_info['chromosomes'].append(chrom)
            chrom_info['lengths'][chrom] = len(arr)
            continue

        # Get FASTA file
        if fasta_dir:
            fa_path = os.path.join(fasta_dir, f"{chrom}.fa")
            if not os.path.exists(fa_path):
                print(f"  WARNING: {fa_path} not found, skipping {chrom}", file=sys.stderr)
                continue
        else:
            fa_path = os.path.join(tmp_dir, f"{chrom}.fa")
            if not os.path.exists(fa_path):
                if not download_chromosome(chrom, fa_path):
                    continue

        # Read and encode
        print(f"  Encoding {chrom}...", file=sys.stderr)
        t0 = time.time()
        header, seq = read_fasta(fa_path)
        arr = encode_genome(seq)
        np.save(npy_path, arr)
        elapsed = time.time() - t0
        print(f"  {chrom}: {len(arr):,} bp encoded in {elapsed:.1f}s", file=sys.stderr)

        chrom_info['chromosomes'].append(chrom)
        chrom_info['lengths'][chrom] = len(arr)

        # Clean up FASTA if we downloaded it
        if not fasta_dir and os.path.exists(fa_path):
            os.unlink(fa_path)

    # Write chromosome info
    info_path = os.path.join(output_dir, "chrom_info.json")
    with open(info_path, 'w') as f:
        json.dump(chrom_info, f, indent=2)

    total_bp = sum(chrom_info['lengths'].values())
    print(f"\nGenome preparation complete:", file=sys.stderr)
    print(f"  Chromosomes: {len(chrom_info['chromosomes'])}", file=sys.stderr)
    print(f"  Total bases: {total_bp:,}", file=sys.stderr)
    print(f"  Output dir:  {output_dir}", file=sys.stderr)

    # Clean up tmp dir
    if os.path.isdir(tmp_dir) and not os.listdir(tmp_dir):
        os.rmdir(tmp_dir)


def main():
    parser = argparse.ArgumentParser(
        description="Prepare hg19 genome for ePCR screening"
    )
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for .npy genome files')
    parser.add_argument('--chromosomes', nargs='*', default=None,
                        help='Chromosomes to prepare (default: chr1-22, chrX, chrY)')
    parser.add_argument('--fasta-dir', default=None,
                        help='Directory with existing per-chromosome FASTA files')
    args = parser.parse_args()

    prepare_genome(
        output_dir=args.output_dir,
        chromosomes=args.chromosomes,
        fasta_dir=args.fasta_dir,
    )


if __name__ == '__main__':
    main()
