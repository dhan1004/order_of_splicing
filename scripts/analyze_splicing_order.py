#!/usr/bin/env python3

import pysam
import argparse
from collections import defaultdict
import sys
import gzip

def load_introns_simple(intron_bed_file):
    """
    Load introns into a dict
    Returns: (intron_dict, gene_dict)
    intron_dict: {(chr, start, end): gene_id}
    gene_dict: {gene_id: [(chr, start, end), ...]}
    """
    intron_dict = {}
    gene_dict = defaultdict(list)
    
    opener = gzip.open if intron_bed_file.endswith(".gz") else open
    
    n_introns = 0
    with opener(intron_bed_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            
            # parse gene ID
            gene_info = fields[3]
            if '_' in gene_info:
                gene_id = gene_info.split('_')[0]
            elif '|' in gene_info:
                gene_id = gene_info.split('|')[0]
            else:
                gene_id = gene_info
            
            key = (chrom, start, end)
            intron_dict[key] = gene_id
            gene_dict[gene_id].append(key)
            n_introns += 1
    
    print(f"Loaded {n_introns:,} introns from {len(gene_dict):,} genes", file=sys.stderr)
    return intron_dict, gene_dict

def fuzzy_match(query_coord, intron_dict, tolerance):
    """
    fuzzy matching within tolerance
    """
    q_chr, q_start, q_end = query_coord
    
    # check exact match first
    if query_coord in intron_dict:
        return query_coord, intron_dict[query_coord]
    
    # check nearby coordinates within tolerance
    for start_offset in range(-tolerance, tolerance + 1):
        for end_offset in range(-tolerance, tolerance + 1):
            test_key = (q_chr, q_start + start_offset, q_end + end_offset)
            if test_key in intron_dict:
                return test_key, intron_dict[test_key]
    
    return None, None

def extract_junctions_from_cigar(read):
    """
    Extract splice junctions from CIGAR string
    Returns: [(chr, start, end), ...]
    """
    if read.is_unmapped or not read.cigartuples:
        return []
    
    junctions = []
    pos = read.reference_start
    chrom = read.reference_name
    
    for op, length in read.cigartuples:
        if op == 3:  # N = splice junction
            junctions.append((chrom, pos, pos + length))
        if op in [0, 2, 3, 7, 8]: 
            pos += length
    
    return junctions

def get_read_span(read):
    """
    Get genomic span of read
    """
    if read.is_unmapped:
        return None
    return (read.reference_name, read.reference_start, read.reference_end)

def analyze_pair_minimal(read1, read2, intron_dict, gene_dict, tolerance):
    """
    find introns and count evidence for splicing order in this read pair
    """
    evidence = []
    
    # extract junctions from both reads
    junctions = extract_junctions_from_cigar(read1) + extract_junctions_from_cigar(read2)
    if not junctions:
        return evidence
    
    # match junctions to annotated introns
    matched_junctions = []
    for junction in junctions:
        matched_intron, gene_id = fuzzy_match(junction, intron_dict, tolerance)
        if matched_intron:
            matched_junctions.append((matched_intron, gene_id))
    
    if not matched_junctions:
        return evidence
    
    # Get read spans
    spans = []
    for read in [read1, read2]:
        span = get_read_span(read)
        if span:
            spans.append(span)
    
    # For each matched junction, find retained introns in the same gene
    for (splice_chr, splice_start, splice_end), gene_id in matched_junctions:
        
        # Get all introns for this gene
        gene_introns = gene_dict.get(gene_id, [])
        
        # Check which ones overlap our read spans
        for span_chr, span_start, span_end in spans:
            if span_chr != splice_chr:
                continue
            
            for intron_chr, intron_start, intron_end in gene_introns:
                # Skip if it's the same intron
                if intron_start == splice_start and intron_end == splice_end:
                    continue
                
                # Check if span overlaps this intron
                if span_start < intron_end and span_end > intron_start:
                    # Found evidence! Determine order
                    if splice_start < intron_start:
                        # Upstream spliced, downstream retained
                        evidence.append((
                            gene_id,
                            (splice_start, splice_end),
                            (intron_start, intron_end),
                            'upstream'
                        ))
                    else:
                        # Downstream spliced, upstream retained
                        evidence.append((
                            gene_id,
                            (intron_start, intron_end),
                            (splice_start, splice_end),
                            'downstream'
                        ))
    
    return evidence

def main():
    parser = argparse.ArgumentParser(
        description='Zero-redundancy splicing order counter - trusts pre-filtered BAM'
    )
    parser.add_argument('--bam', required=True, help='Informative pairs BAM (pre-filtered)')
    parser.add_argument('--intron-bed', required=True, help='Intron annotations BED')
    parser.add_argument('--output', required=True, help='Output TSV')
    parser.add_argument('--min-reads', type=int, default=5, help='Min read support')
    parser.add_argument('--min-mapq', type=int, default=10, help='Min MAPQ filter')
    parser.add_argument('--tolerance', type=int, default=10, help='Junction match tolerance')
    parser.add_argument('--max-intron-length', type=int, default=None, 
                    help='Maximum intron length (bp) - filter intron pairs')
    args = parser.parse_args()
    
    # Load annotations
    print("Loading intron annotations...", file=sys.stderr)
    intron_dict, gene_dict = load_introns_simple(args.intron_bed)
    
    # Process BAM - single pass, group by read name
    print(f"Processing BAM: {args.bam}", file=sys.stderr)
    pairs = defaultdict(list)
    n_reads = 0
    n_filtered = 0
    
    with pysam.AlignmentFile(args.bam, 'rb') as bam:
        for read in bam:
            n_reads += 1
            
            # Basic quality filter
            if read.is_unmapped or read.mapping_quality < args.min_mapq:
                n_filtered += 1
                continue
            
            pairs[read.query_name].append(read)
            
            if n_reads % 100000 == 0:
                print(f"  {n_reads:,} reads, {len(pairs):,} pairs", file=sys.stderr)
    
    print(f"  Total reads: {n_reads:,}", file=sys.stderr)
    print(f"  Filtered (unmapped/low MAPQ): {n_filtered:,}", file=sys.stderr)
    print(f"  Valid pairs: {len(pairs):,}", file=sys.stderr)
    
    # Count evidence
    print("Counting evidence...", file=sys.stderr)
    evidence_counts = defaultdict(lambda: {'upstream': 0, 'downstream': 0})
    
    n_pairs = 0
    n_evidence = 0
    
    for read_name, read_list in pairs.items():
        # Skip if not proper pair
        if len(read_list) != 2:
            continue
        
        n_pairs += 1
        if n_pairs % 10000 == 0:
            print(f"  {n_pairs:,} pairs analyzed, {n_evidence:,} evidence events", 
                  file=sys.stderr)
        
        # Analyze this pair
        evidence = analyze_pair_minimal(
            read_list[0], read_list[1],
            intron_dict, gene_dict,
            args.tolerance
        )
        
        # Count evidence
        for gene_id, intron1, intron2, direction in evidence:
            # Use first read's chr for output
            if args.max_intron_length:
                intron1_len = abs(intron1[1] - intron1[0])  # end - start
                intron2_len = abs(intron2[1] - intron2[0])
                if intron1_len > args.max_intron_length or intron2_len > args.max_intron_length:
                    continue
            chrom = read_list[0].reference_name
            key = f"{chrom}:{gene_id}:{intron1[0]}:{intron1[1]}:{intron2[0]}:{intron2[1]}"
            evidence_counts[key][direction] += 1
            n_evidence += 1
    
    print(f"  Total pairs analyzed: {n_pairs:,}", file=sys.stderr)
    print(f"  Total evidence events: {n_evidence:,}", file=sys.stderr)
    print(f"  Unique intron pairs: {len(evidence_counts):,}", file=sys.stderr)
    
    # Write output
    print(f"Writing results to {args.output}...", file=sys.stderr)
    n_written = 0
    
    with open(args.output, 'w') as out:
        # Header
        out.write("chr\tgene_id\tintron1_start\tintron1_end\t"
                 "intron2_start\tintron2_end\t"
                 "upstream\tdownstream\ttotal\tfraction_downstream\n")
        
        # Data
        for key in sorted(evidence_counts.keys()):
            parts = key.split(':')
            chrom = parts[0]
            gene_id = parts[1]
            coords = list(map(int, parts[2:]))
            
            upstream = evidence_counts[key]['upstream']
            downstream = evidence_counts[key]['downstream']
            total = upstream + downstream
            
            # Filter by minimum reads
            if total >= args.min_reads:
                fraction_down = downstream / total if total > 0 else 0.5
                out.write(f"{chrom}\t{gene_id}\t"
                         f"{coords[0]}\t{coords[1]}\t{coords[2]}\t{coords[3]}\t"
                         f"{upstream}\t{downstream}\t{total}\t{fraction_down:.4f}\n")
                n_written += 1
    
    print(f"Done! Wrote {n_written:,} intron pairs with ≥{args.min_reads} reads", 
          file=sys.stderr)
    print(f"Output: {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()