#!/usr/bin/env python3
"""
pairwise_splice_site_spanning.py

Measures 5'ss vs 3'ss continuous coverage on the LATER-SPLICED (retained) intron
of each informative read pair, keyed to the same pairs used to call splicing
order. Human-first (hg38), designed to run on the informative_pairs BAM.

WHY PAIRWISE (not per-intron):
The frozen-intermediate question is about ORDER. When the downstream intron
splices first, the still-retained upstream intron is the pre-mRNA side, and its
splicing STEP is what we probe: an intron frozen after step 1 (branching) but
before step 2 has an intact 3'ss (intron still joined to 3' exon) but a broken
5'ss (5' exon already released). So a 3'ss > 5'ss spanning skew on the retained
intron, within the pairs where it is the later-spliced one, is the signature.
A per-intron aggregate would average over pairing contexts and wash this out.

DISJOINT READS (avoids circularity):
In each informative pair, the JUNCTION read (CIGAR has 'N') identifies which
intron already spliced -> makes the order call. The MATE read overlaps the
retained intron. We measure spanning on the MATE read only. Order call and
spanning phenotype therefore come from physically different reads in the pair,
mirroring the disjoint-read-set principle used for the skipping analysis.

This script mirrors analyze_splicing_order.py's pair logic (same fuzzy_match,
same junction extraction, same evidence direction convention) so the pairs it
scores are exactly the informative pairs the order call is built on.

OUTPUT (one row per intron pair, keyed like splicing_order_pooled.tsv):
  chr gene_id
  retained_intron_start retained_intron_end retained_strand
  spliced_intron_start  spliced_intron_end
  direction                      # 'upstream' or 'downstream' (which spliced 1st)
  fivep_span fivep_over threep_span threep_over    # on retained intron, mate reads
  fivep_frac threep_frac span_ratio_3p_over_5p

Pooling across samples is done by pool_pairwise_spanning (sums counts per
(pair, direction) key, then recomputes fractions from pooled sums).

Usage:
  python3 pairwise_splice_site_spanning.py \
      --bam SAMPLE_informative_pairs.bam \
      --intron-bed /users/dhan30/reference/hg38.gencode.basic.v43.introns.bed.gz \
      --output SAMPLE_pairwise_spanning.tsv \
      --min-mapq 10 --tolerance 10 --k 10
"""
import argparse
import gzip
import sys
from collections import defaultdict

import pysam

_MATCH = {0, 7, 8}   # M, =, X consume query+ref as aligned
_REF_SKIP = 3        # N


# --- annotation loading: same shape as analyze_splicing_order, plus strand ---
def load_introns(intron_bed_file):
    """Return intron_dict{(chr,s,e):gene}, gene_dict{gene:[(chr,s,e)]},
    strand_dict{(chr,s,e):strand}."""
    intron_dict = {}
    gene_dict = defaultdict(list)
    strand_dict = {}
    opener = gzip.open if intron_bed_file.endswith(".gz") else open
    n = 0
    with opener(intron_bed_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            gene_info = fields[3]
            if "_" in gene_info:
                gene_id = gene_info.split("_")[0]
            elif "|" in gene_info:
                gene_id = gene_info.split("|")[0]
            else:
                gene_id = gene_info
            strand = fields[5] if len(fields) >= 6 else "+"
            key = (chrom, start, end)
            intron_dict[key] = gene_id
            gene_dict[gene_id].append(key)
            strand_dict[key] = strand
            n += 1
    print(f"Loaded {n:,} introns from {len(gene_dict):,} genes", file=sys.stderr)
    return intron_dict, gene_dict, strand_dict


def fuzzy_match(query_coord, intron_dict, tolerance):
    if query_coord in intron_dict:
        return query_coord, intron_dict[query_coord]
    q_chr, q_start, q_end = query_coord
    for so in range(-tolerance, tolerance + 1):
        for eo in range(-tolerance, tolerance + 1):
            test = (q_chr, q_start + so, q_end + eo)
            if test in intron_dict:
                return test, intron_dict[test]
    return None, None


def extract_junctions(read):
    if read.is_unmapped or not read.cigartuples:
        return []
    junctions = []
    pos = read.reference_start
    chrom = read.reference_name
    for op, length in read.cigartuples:
        if op == _REF_SKIP:
            junctions.append((chrom, pos, pos + length))
        if op in (0, 2, 3, 7, 8):
            pos += length
    return junctions


def read_spans_boundary(read, boundary, k):
    """True if read has matched bases on both sides of `boundary` within k bp,
    with no N block crossing the boundary. Same logic as splice_site_spanning."""
    left_lo, left_hi = boundary - k, boundary
    right_lo, right_hi = boundary, boundary + k
    has_left = has_right = False
    pos = read.reference_start
    for op, length in read.cigartuples:
        if op in _MATCH:
            seg_start, seg_end = pos, pos + length
            if seg_start < left_hi and seg_end > left_lo:
                has_left = True
            if seg_start < right_hi and seg_end > right_lo:
                has_right = True
            pos += length
        elif op == _REF_SKIP:
            if pos <= boundary < pos + length:
                return False       # spliced here, not spanning
            pos += length
        elif op == 2:              # D consumes ref
            pos += length
    return has_left and has_right


def read_overlaps_boundary(read, boundary):
    if read.is_unmapped or read.reference_end is None:
        return False
    return read.reference_start < boundary < read.reference_end


def which_is_junction_read(read1, read2):
    """Return (junction_read, mate_read) or (None, None). The junction read is
    the one whose CIGAR carries an N (it proves an intron spliced). If both or
    neither carry N, we can't cleanly separate -> skip for disjoint measurement.
    """
    n1 = read1.cigartuples is not None and any(op == _REF_SKIP for op, _ in read1.cigartuples)
    n2 = read2.cigartuples is not None and any(op == _REF_SKIP for op, _ in read2.cigartuples)
    if n1 and not n2:
        return read1, read2
    if n2 and not n1:
        return read2, read1
    return None, None


def score_pair(read1, read2, intron_dict, gene_dict, strand_dict, tolerance, k):
    """
    Reproduce analyze_pair_minimal's classification, then measure spanning on
    the retained intron using the MATE (non-junction) read only.

    Yields dicts, one per evidence event.
    """
    out = []
    junction_read, mate_read = which_is_junction_read(read1, read2)
    if junction_read is None:
        return out   # can't do disjoint measurement on this pair

    junctions = extract_junctions(junction_read)
    if not junctions:
        return out

    matched = []
    for j in junctions:
        mi, gid = fuzzy_match(j, intron_dict, tolerance)
        if mi:
            matched.append((mi, gid))
    if not matched:
        return out

    mate_span = None
    if not mate_read.is_unmapped and mate_read.reference_end is not None:
        mate_span = (mate_read.reference_name,
                     mate_read.reference_start, mate_read.reference_end)
    if mate_span is None:
        return out

    for (splice_chr, splice_start, splice_end), gene_id in matched:
        for intron_chr, intron_start, intron_end in gene_dict.get(gene_id, []):
            if intron_start == splice_start and intron_end == splice_end:
                continue  # same intron as the junction
            # retained intron must be overlapped by the MATE read
            if mate_span[0] != splice_chr:
                continue
            if not (mate_span[1] < intron_end and mate_span[2] > intron_start):
                continue

            # direction convention matches analyze_splicing_order.py:
            # retained intron downstream of spliced -> 'downstream' spliced first
            if splice_start < intron_start:
                direction = "upstream"   # upstream (spliced) first; downstream retained
            else:
                direction = "downstream" # downstream (spliced) first; upstream retained

            strand = strand_dict.get((intron_chr, intron_start, intron_end), "+")

            # 5'ss / 3'ss boundary on the RETAINED intron, strand-aware
            if strand == "-":
                fivep_boundary = intron_end
                threep_boundary = intron_start
            else:
                fivep_boundary = intron_start
                threep_boundary = intron_end

            # measure spanning on the MATE read only (disjoint from order call)
            fo = read_overlaps_boundary(mate_read, fivep_boundary)
            to = read_overlaps_boundary(mate_read, threep_boundary)
            fs = fo and read_spans_boundary(mate_read, fivep_boundary, k)
            ts = to and read_spans_boundary(mate_read, threep_boundary, k)

            # only emit if the mate actually reaches at least one boundary
            if not (fo or to):
                continue

            out.append({
                "chr": intron_chr, "gene_id": gene_id,
                "retained_intron_start": intron_start,
                "retained_intron_end": intron_end,
                "retained_strand": strand,
                "spliced_intron_start": splice_start,
                "spliced_intron_end": splice_end,
                "direction": direction,
                "fivep_span": int(bool(fs)), "fivep_over": int(bool(fo)),
                "threep_span": int(bool(ts)), "threep_over": int(bool(to)),
            })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True, help="informative_pairs BAM")
    ap.add_argument("--intron-bed", required=True,
                    help="hg38 intron BED with strand in col 6 (.bed/.bed.gz)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--tolerance", type=int, default=10)
    ap.add_argument("--k", type=int, default=10, help="bp each side of a junction")
    args = ap.parse_args()

    intron_dict, gene_dict, strand_dict = load_introns(args.intron_bed)

    # group reads by name (same as analyze_splicing_order)
    pairs = defaultdict(list)
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.mapping_quality < args.min_mapq:
                continue
            pairs[read.query_name].append(read)

    # accumulate per (pair, direction) key
    acc = defaultdict(lambda: {"fivep_span": 0, "fivep_over": 0,
                               "threep_span": 0, "threep_over": 0,
                               "gene_id": None, "retained_strand": "+"})
    n_pairs = 0
    for read_name, rl in pairs.items():
        if len(rl) != 2:
            continue
        n_pairs += 1
        for ev in score_pair(rl[0], rl[1], intron_dict, gene_dict,
                              strand_dict, args.tolerance, args.k):
            key = (ev["chr"], ev["retained_intron_start"],
                   ev["retained_intron_end"], ev["spliced_intron_start"],
                   ev["spliced_intron_end"], ev["direction"])
            a = acc[key]
            a["fivep_span"] += ev["fivep_span"]
            a["fivep_over"] += ev["fivep_over"]
            a["threep_span"] += ev["threep_span"]
            a["threep_over"] += ev["threep_over"]
            a["gene_id"] = ev["gene_id"]
            a["retained_strand"] = ev["retained_strand"]

    with open(args.output, "w") as out:
        out.write("\t".join([
            "chr", "gene_id",
            "retained_intron_start", "retained_intron_end", "retained_strand",
            "spliced_intron_start", "spliced_intron_end", "direction",
            "fivep_span", "fivep_over", "threep_span", "threep_over",
            "fivep_frac", "threep_frac", "span_ratio_3p_over_5p",
        ]) + "\n")
        for key, a in sorted(acc.items()):
            chrom, ris, rie, sis, sie, direction = key
            ffrac = a["fivep_span"] / a["fivep_over"] if a["fivep_over"] else float("nan")
            tfrac = a["threep_span"] / a["threep_over"] if a["threep_over"] else float("nan")
            ratio = (tfrac / ffrac) if (a["fivep_over"] and a["threep_over"]
                                        and ffrac > 0) else float("nan")
            out.write("\t".join(str(x) for x in [
                chrom, a["gene_id"], ris, rie, a["retained_strand"],
                sis, sie, direction,
                a["fivep_span"], a["fivep_over"], a["threep_span"], a["threep_over"],
                "%.4f" % ffrac if a["fivep_over"] else "NA",
                "%.4f" % tfrac if a["threep_over"] else "NA",
                "%.4f" % ratio if ratio == ratio else "NA",
            ]) + "\n")

    print(f"Scored {n_pairs:,} read pairs -> {len(acc):,} (pair,direction) rows",
          file=sys.stderr)
    print(f"Written: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()