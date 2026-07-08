#!/usr/bin/env python3
"""
build_introns_bed.py

Generate an introns BED matching the exact format of
hg38.gencode.basic.v43.introns.bed.gz, from a GENCODE GTF.

Format (tab-separated, 6 columns):
    chrom  start  end  NAME  0  strand
where NAME = {transcript_id}_intron_{n}_0_{chrom}_{start+1}_{strandchar}
  - {n} counts introns in TRANSCRIPT order from the 5' end (0-based)
    (so on minus-strand genes, intron 0 has the HIGHEST genomic coord)
  - start+1 is the 1-based start coordinate
  - strandchar is 'f' for '+' and 'r' for '-'

Introns are derived by taking each transcript's exons, sorting by genomic
coordinate, and taking the gaps between consecutive exons.

Usage:
    python build_introns_bed.py \
        --gtf  /users/dhan30/reference/gencode.vM36.basic.annotation.gtf \
        --out  /users/dhan30/reference/mm39.gencode.basic.vM36.introns.bed.gz
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict


def open_maybe_gz(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def parse_attr(attr_field, key):
    m = re.search(key + r' "([^"]+)"', attr_field)
    return m.group(1) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    # Collect exons per transcript: tx_id -> list of (start0, end, chrom, strand)
    tx_exons = defaultdict(list)
    tx_meta = {}  # tx_id -> (chrom, strand)

    n_lines = 0
    with open_maybe_gz(args.gtf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'exon':
                continue
            chrom = f[0]
            start0 = int(f[3]) - 1   # GTF is 1-based inclusive -> BED 0-based start
            end = int(f[4])          # end stays (1-based inclusive == 0-based exclusive)
            strand = f[6]
            tx_id = parse_attr(f[8], 'transcript_id')
            if tx_id is None:
                continue
            tx_exons[tx_id].append((start0, end))
            tx_meta[tx_id] = (chrom, strand)
            n_lines += 1

    print(f"Parsed {n_lines:,} exon lines across {len(tx_exons):,} transcripts",
          file=sys.stderr)

    n_introns = 0
    with open_maybe_gz(args.out, 'wt') as out:
        for tx_id, exons in tx_exons.items():
            chrom, strand = tx_meta[tx_id]
            # sort exons by genomic coordinate
            exons = sorted(set(exons))
            if len(exons) < 2:
                continue

            # introns = gaps between consecutive exons in genomic order
            genomic_introns = []
            for (s1, e1), (s2, e2) in zip(exons[:-1], exons[1:]):
                istart = e1        # first base after exon end (0-based)
                iend = s2          # exon2 start (0-based) == exclusive end
                if iend > istart:  # skip zero-length / overlapping
                    genomic_introns.append((istart, iend))

            if not genomic_introns:
                continue

            # Number in TRANSCRIPT order (5'->3'):
            #   plus strand: ascending genomic == transcript order
            #   minus strand: descending genomic == transcript order
            if strand == '+':
                ordered = list(enumerate(genomic_introns))
            else:
                ordered = list(enumerate(reversed(genomic_introns)))

            strandchar = 'f' if strand == '+' else 'r'

            # But we still emit rows; the human file emits in genomic order of
            # appearance (descending intron number for minus strand). To match,
            # sort emitted rows by the same order as the source file: the human
            # file lists them in genomic ascending order for the coords BUT with
            # transcript-order numbers. We'll emit in genomic ascending order.
            # Build (istart, iend, intron_num) then sort by istart.
            numbered = []
            for num, (istart, iend) in ordered:
                numbered.append((istart, iend, num))
            numbered.sort(key=lambda x: x[0])

            for istart, iend, num in numbered:
                name = f"{tx_id}_intron_{num}_0_{chrom}_{istart + 1}_{strandchar}"
                out.write(f"{chrom}\t{istart}\t{iend}\t{name}\t0\t{strand}\n")
                n_introns += 1

    print(f"Wrote {n_introns:,} introns -> {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()