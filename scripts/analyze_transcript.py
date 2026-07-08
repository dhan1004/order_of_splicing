#!/usr/bin/env python3
"""
Unified RNA Structure Analysis for Splicing Order

Metric conventions:
  - intron1 = UPSTREAM intron (lower genomic coordinates)
  - intron2 = DOWNSTREAM intron (higher genomic coordinates)
  - All difference metrics = downstream - upstream (intron2 - intron1)
      Positive value → downstream intron has MORE of that feature
  - H-bond counts normalized by intron length (_per_nt columns)
  - Accessibility metrics are already fractions [0,1], no normalization needed

Exon flanks are computed proportionally per intron:
  flank = clamp(0.10 * intron_length, min=200bp, max=500bp)

Author: Squash
"""

import argparse
import sys
from pathlib import Path
import subprocess
import tempfile
from collections import defaultdict
import re
import pandas as pd
import pysam

def fold_sequence_rnafold(sequence, temperature=37.0):
    """
    Fold sequence using RNAfold via string stdin
    """
    try:
        input_str = f">seq\n{sequence}\n"
        result = subprocess.run(
            ['RNAfold', '--noPS', f'--temp={temperature}'],
            input=input_str,
            capture_output=True,
            text=True,
            check=True
        )
        lines = result.stdout.strip().split('\n')
        # RNAfold output: line 0 = ">seq", line 1 = sequence, line 2 = structure + MFE
        if len(lines) >= 3:
            match = re.match(r'([.()\[\]]+)\s+\(\s*(-?\d+\.\d+)\s*\)', lines[2])
            if match:
                return match.group(1), float(match.group(2))
        print(f"WARNING: RNAfold parsing failed for sequence length {len(sequence)}: {result.stdout[:200]}",
              file=sys.stderr)
        return None, None
    except subprocess.CalledProcessError as e:
        print(f"ERROR: RNAfold failed: {e.stderr[:200]}", file=sys.stderr)
        return None, None
    except Exception as e:
        print(f"ERROR: fold_sequence_rnafold: {e}", file=sys.stderr)
        return None, None


def parse_dot_bracket(structure):
    stack, pairs = [], []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')' and stack:
            pairs.append((stack.pop(), i))
    return pairs


def count_hydrogen_bonds_from_pairs(sequence, pairs):
    hbond_count = gc_pairs = au_pairs = gu_pairs = 0
    for i, j in pairs:
        pair = ''.join(sorted([sequence[i].upper(), sequence[j].upper()]))
        if pair == 'CG':
            hbond_count += 3; gc_pairs += 1
        elif pair == 'AU':
            hbond_count += 2; au_pairs += 1
        elif pair == 'GT':
            hbond_count += 2; gu_pairs += 1
    return hbond_count, gc_pairs, au_pairs, gu_pairs


# SLIDING WINDOW

def sliding_window_fold(sequence, window_size=300, step_size=150):
    all_pairs = set()
    window_mfes = []
    attempted = succeeded = 0

    for start in range(0, len(sequence), step_size):
        end = min(start + window_size, len(sequence))
        if end - start < 30:
            break
        attempted += 1
        structure, mfe = fold_sequence_rnafold(sequence[start:end])
        if structure:
            succeeded += 1
            window_mfes.append(mfe)
            for i, j in parse_dot_bracket(structure):
                all_pairs.add((start + i, start + j))
        elif attempted == 1:
            print(f"  WARNING: RNAfold failed for first window {start}-{end}", file=sys.stderr)

    if attempted > 0 and succeeded / attempted < 0.5:
        print(f"  WARNING: Low RNAfold success: {succeeded}/{attempted}", file=sys.stderr)

    return all_pairs, window_mfes


# Exon flank calculation

def compute_exon_flank(intron_length, min_flank=100, max_flank=1000, fraction=0.10):
    """
    Proportional exon flank: clamp(fraction * intron_length, min_flank, max_flank).
    Each intron gets its own flank based on its own length.
    """
    return int(max(min_flank, min(max_flank, fraction * intron_length)))


# boundary bond analysis

def identify_boundary_hbonds(sequence, pairs, region_boundaries):
    """
    Count hydrogen bonds crossing each exon-intron boundary.

    Junction bonds per intron
      upstream intron (intron1):   upstream_exon/intron1  +  intron1/middle_exon
      downstream intron (intron2): middle_exon/intron2    +  intron2/downstream_exon
    """
    ue_end  = region_boundaries['upstream_exon_end']
    i1_end  = region_boundaries['intron1_end']
    me_end  = region_boundaries['middle_exon_end']
    i2_end  = region_boundaries['intron2_end']

    def get_region(pos):
        if pos < ue_end:   return 'upstream_exon'
        elif pos < i1_end: return 'intron1'
        elif pos < me_end: return 'middle_exon'
        elif pos < i2_end: return 'intron2'
        else:              return 'downstream_exon'

    counts = defaultdict(int)
    for i, j in pairs:
        pair = ''.join(sorted([sequence[i].upper(), sequence[j].upper()]))
        if pair == 'CG':           hbonds = 3
        elif pair in ('AU', 'GT'): hbonds = 2
        else: continue
        key = tuple(sorted([get_region(i), get_region(j)]))
        if key == ('intron1', 'upstream_exon'):      counts['intron1_5ss'] += hbonds
        elif key == ('intron1', 'middle_exon'):      counts['intron1_3ss'] += hbonds
        elif key == ('intron2', 'middle_exon'):      counts['intron2_5ss'] += hbonds
        elif key == ('downstream_exon', 'intron2'):  counts['intron2_3ss'] += hbonds
        elif key == ('intron1', 'intron2'):          counts['intron1_intron2'] += hbonds

    intron1_junction = counts['intron1_5ss'] + counts['intron1_3ss']
    intron2_junction = counts['intron2_5ss'] + counts['intron2_3ss']

    return {
        # Individual boundary counts (raw)
        'intron1_5ss_hbonds':      counts['intron1_5ss'],
        'intron1_3ss_hbonds':      counts['intron1_3ss'],
        'intron1_junction_hbonds': intron1_junction,
        'intron2_5ss_hbonds':      counts['intron2_5ss'],
        'intron2_3ss_hbonds':      counts['intron2_3ss'],
        'intron2_junction_hbonds': intron2_junction,
        'intron1_intron2_hbonds':  counts['intron1_intron2'],
        # Legacy names for backwards compatibility
        'intron1_exon_hbonds':              intron1_junction,
        'intron2_exon_hbonds':              intron2_junction,
        'upstream_exon_intron1_hbonds':     counts['intron1_5ss'],
        'intron1_middle_exon_hbonds':       counts['intron1_3ss'],
        'middle_exon_intron2_hbonds':       counts['intron2_5ss'],
        'intron2_downstream_exon_hbonds':   counts['intron2_3ss'],
    }


# internal structure analysis

def analyze_internal_structure(sequence, pairs, intron_start, intron_end,
                                region_start, exon_flank):
    seq_start  = intron_start - region_start
    seq_end    = intron_end   - region_start
    intron_len = seq_end - seq_start
    intron_seq = sequence[seq_start:seq_end]

    intron_pairs = [(i, j) for i, j in pairs
                    if seq_start <= i < seq_end and seq_start <= j < seq_end]
    local_pairs  = [(i - seq_start, j - seq_start) for i, j in intron_pairs]

    total_hbonds, gc_p, au_p, gu_p = count_hydrogen_bonds_from_pairs(intron_seq, local_pairs)

    paired_positions = {p for pair in pairs for p in pair}
    ss5_pos = seq_start
    ss3_pos = seq_end - 1

    def accessibility(pos, window):
        s = max(0, pos - window)
        e = min(len(sequence), pos + window + 1)
        unpaired = sum(1 for p in range(s, e) if p not in paired_positions)
        return unpaired / (e - s) if e > s else float('nan')

    gc_content       = (intron_seq.count('G') + intron_seq.count('C')) / intron_len if intron_len > 0 else float('nan')
    pairing_fraction = len(intron_pairs) * 2 / intron_len if intron_len > 0 else float('nan')

    return {
        'length':                 intron_len,
        'gc_content':             gc_content,
        'pairing_fraction':       pairing_fraction,
        'total_hbonds':           total_hbonds,
        'gc_pairs':               gc_p,
        'au_pairs':               au_p,
        'gu_pairs':               gu_p,
        'total_pairs':            len(intron_pairs),
        'ss5_accessibility_10nt': accessibility(ss5_pos, 10),
        'ss5_accessibility_25nt': accessibility(ss5_pos, 25),
        'ss5_accessibility_50nt': accessibility(ss5_pos, 50),
        'ss5_paired':             int(ss5_pos in paired_positions),
        'ss3_accessibility_10nt': accessibility(ss3_pos, 10),
        'ss3_accessibility_25nt': accessibility(ss3_pos, 25),
        'ss3_accessibility_50nt': accessibility(ss3_pos, 50),
        'ss3_paired':             int(ss3_pos in paired_positions),
    }

def extract_premrna_sequence(genome_fasta, chrom, start, end, strand='+'):
    bed_string = f"{chrom}\t{start}\t{end}\t.\t.\t{strand}\n"
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write(bed_string)
        bed_file = f.name
    try:
        result = subprocess.run(
            ['bedtools', 'getfasta', '-fi', genome_fasta, '-bed', bed_file, '-s', '-tab'],
            capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split('\n')
        if lines:
            parts = lines[0].split('\t')
            if len(parts) == 2:
                return parts[1].upper().replace('T', 'U')
        return None
    finally:
        Path(bed_file).unlink(missing_ok=True)

def safe_div(a, b):
    try:
        return a / b if b and b > 0 else float('nan')
    except Exception:
        return float('nan')


# main processing

def process_intron_pair(row, genome_fasta, window_size=200, step_size=50, exon_flank=None):
    """
    Process one intron pair.

    intron1 = upstream intron  (lower genomic coords)
    intron2 = downstream intron (higher genomic coords)

    Exon flanks: clamp(0.10 * intron_length, 200, 500) computed per intron.
    The exon_flank argument is accepted for backwards compatibility but ignored.

    All diff_* columns = downstream - upstream (intron2 - intron1).
    All _per_nt columns = raw count / intron length.
    """
    chrom         = row['chr']
    intron1_start = int(row['intron1_start'])
    intron1_end   = int(row['intron1_end'])
    intron2_start = int(row['intron2_start'])
    intron2_end   = int(row['intron2_end'])
    strand        = row.get('strand', '+')

    intron1_len = intron1_end - intron1_start
    intron2_len = intron2_end - intron2_start

    # Proportional flanks - each intron gets its own
    upstream_flank   = compute_exon_flank(intron1_len)
    downstream_flank = compute_exon_flank(intron2_len)

    region_start = intron1_start - upstream_flank
    region_end   = intron2_end   + downstream_flank

    sequence = extract_premrna_sequence(genome_fasta, chrom, region_start, region_end, strand)
    if not sequence or len(sequence) < 100:
        return None

    middle_exon_len = intron2_start - intron1_end
    region_boundaries = {
        'upstream_exon_end': upstream_flank,
        'intron1_end':       upstream_flank + intron1_len,
        'middle_exon_end':   upstream_flank + intron1_len + middle_exon_len,
        'intron2_end':       upstream_flank + intron1_len + middle_exon_len + intron2_len,
    }

    all_pairs, window_mfes = sliding_window_fold(sequence, window_size, step_size)

    if len(window_mfes) == 0:
        print(f"  WARNING: No windows folded for {chrom}:{region_start}-{region_end}", file=sys.stderr)

    boundary = identify_boundary_hbonds(sequence, all_pairs, region_boundaries)
    intron1_internal = analyze_internal_structure(
        sequence, all_pairs, intron1_start, intron1_end, region_start, upstream_flank)
    intron2_internal = analyze_internal_structure(
        sequence, all_pairs, intron2_start, intron2_end, region_start, downstream_flank)

    # Length-normalized H-bond metrics
    i1_hbonds_per_nt = safe_div(intron1_internal['total_hbonds'], intron1_len)
    i2_hbonds_per_nt = safe_div(intron2_internal['total_hbonds'], intron2_len)
    i1_junc_per_nt   = safe_div(boundary['intron1_junction_hbonds'], intron1_len)
    i2_junc_per_nt   = safe_div(boundary['intron2_junction_hbonds'], intron2_len)
    i1_5ss_per_nt    = safe_div(boundary['intron1_5ss_hbonds'],      intron1_len)
    i1_3ss_per_nt    = safe_div(boundary['intron1_3ss_hbonds'],      intron1_len)
    i2_5ss_per_nt    = safe_div(boundary['intron2_5ss_hbonds'],      intron2_len)
    i2_3ss_per_nt    = safe_div(boundary['intron2_3ss_hbonds'],      intron2_len)

    results = dict(row)
    results = {
        # Identifiers
        'sample_id':       row.get('sample_id', 'unknown'),
        'chr':             chrom,
        'gene_id':         row.get('gene_id', 'unknown'),
        'region_start':    region_start,
        'region_end':      region_end,
        'region_length':   len(sequence),
        'intron1_start':   intron1_start,
        'intron1_end':     intron1_end,
        'intron2_start':   intron2_start,
        'intron2_end':     intron2_end,

        # Splicing order outcome
        'upstream_count':      row.get('upstream', 0),
        'downstream_count':    row.get('downstream', 0),
        'total_reads':         row.get('total', 0),
        'fraction_downstream': row.get('fraction_downstream', float('nan')),

        # Lengths and flanks
        'intron1_length':           intron1_len,
        'intron2_length':           intron2_len,
        'upstream_exon_flank_bp':   upstream_flank,
        'downstream_exon_flank_bp': downstream_flank,

        # Folding metadata
        'num_windows':            len(window_mfes),
        'avg_window_mfe':         sum(window_mfes) / len(window_mfes) if window_mfes else float('nan'),
        'total_pairs_identified': len(all_pairs),
    }

    # Per-intron internal metrics (raw)
    for key, value in intron1_internal.items():
        results[f'intron1_{key}'] = value
    for key, value in intron2_internal.items():
        results[f'intron2_{key}'] = value

    # Per-intron internal metrics (length-normalized)
    results['intron1_hbonds_per_nt'] = i1_hbonds_per_nt
    results['intron2_hbonds_per_nt'] = i2_hbonds_per_nt

    # Boundary metrics (raw + legacy names)
    results.update(boundary)

    # Boundary metrics (length-normalized)
    results['intron1_junction_hbonds_per_nt'] = i1_junc_per_nt
    results['intron1_5ss_hbonds_per_nt']      = i1_5ss_per_nt
    results['intron1_3ss_hbonds_per_nt']       = i1_3ss_per_nt
    results['intron2_junction_hbonds_per_nt'] = i2_junc_per_nt
    results['intron2_5ss_hbonds_per_nt']      = i2_5ss_per_nt
    results['intron2_3ss_hbonds_per_nt']       = i2_3ss_per_nt

    for w in ['10nt', '25nt', '50nt']:
        results[f'diff_ss5_accessibility_{w}'] = (
            intron2_internal[f'ss5_accessibility_{w}'] - intron1_internal[f'ss5_accessibility_{w}'])
        results[f'diff_ss3_accessibility_{w}'] = (
            intron2_internal[f'ss3_accessibility_{w}'] - intron1_internal[f'ss3_accessibility_{w}'])

    results['diff_hbonds_per_nt']          = i2_hbonds_per_nt  - i1_hbonds_per_nt
    results['diff_junction_hbonds_per_nt'] = i2_junc_per_nt    - i1_junc_per_nt
    results['diff_5ss_hbonds_per_nt']      = i2_5ss_per_nt     - i1_5ss_per_nt
    results['diff_3ss_hbonds_per_nt']      = i2_3ss_per_nt     - i1_3ss_per_nt
    results['diff_junction_hbonds_raw']    = boundary['intron2_junction_hbonds'] - boundary['intron1_junction_hbonds']
    results['diff_hbonds_raw']             = intron2_internal['total_hbonds']    - intron1_internal['total_hbonds']

    # Legacy names
    results['ss5_accessibility_difference'] = results['diff_ss5_accessibility_25nt']
    results['ss3_accessibility_difference'] = results['diff_ss3_accessibility_25nt']
    results['internal_hbond_difference']    = results['diff_hbonds_raw']
    results['exon_hbond_difference']        = results['diff_junction_hbonds_raw']
    results['internal_hbond_ratio']         = safe_div(intron1_internal['total_hbonds'],
                                                       intron2_internal['total_hbonds'])

    return results

###############################################################################

def main():
    parser = argparse.ArgumentParser(
        description='RNA structure analysis for splicing order',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Metric conventions:
  intron1_* / intron2_* = per-intron raw values (intron1=upstream, intron2=downstream)
  *_per_nt              = raw count / intron length
  diff_*                = downstream - upstream (intron2 - intron1)
  diff_*_accessibility  = positive → downstream more accessible → positive r expected
  diff_*_hbonds_per_nt  = positive → downstream more occluded  → negative r expected

Exon flanks: clamp(0.10 * intron_length, min=200, max=500) per intron
        """
    )
    parser.add_argument('input_tsv',    help='Input TSV with intron pair coordinates')
    parser.add_argument('genome_fasta', help='Genome FASTA file (indexed)')
    parser.add_argument('-o', '--output',      required=True)
    parser.add_argument('-w', '--window-size', type=int, default=300)
    parser.add_argument('-s', '--step-size',   type=int, default=150)
    parser.add_argument('-e', '--exon-flank',  type=int, default=None,
                        help='Ignored — flanks computed as clamp(0.10 * intron_length, 100, 1000)')
    parser.add_argument('--max-rows',          type=int, default=None)
    args = parser.parse_args()

    for tool in ['RNAfold', 'bedtools']:
        try:
            subprocess.run([tool, '--version'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print(f"ERROR: {tool} not found.", file=sys.stderr); sys.exit(1)

    print(f"Reading {args.input_tsv}...")
    df = pd.read_csv(args.input_tsv, sep='\t')
    if args.max_rows:
        df = df.head(args.max_rows)
    print(f"Processing {len(df):,} intron pairs...")

    results, failed = [], 0
    for idx, row in df.iterrows():
        if idx % 100 == 0:
            print(f"  {idx:,}/{len(df):,}", end='\r')
        try:
            r = process_intron_pair(row, args.genome_fasta, args.window_size, args.step_size)
            if r: results.append(r)
            else: failed += 1
        except Exception as e:
            print(f"\nWarning: row {idx} failed: {e}", file=sys.stderr)
            failed += 1

    print(f"\nProcessed {len(results):,}/{len(df):,} pairs ({failed} failed)")

    if not results:
        print("ERROR: no results", file=sys.stderr); sys.exit(1)

    out_df = pd.DataFrame(results)
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f"Written to {args.output}")

    zero_windows = (out_df['num_windows'] == 0).sum()
    if zero_windows > 0:
        print(f"WARNING: {zero_windows} rows have num_windows=0 — RNAfold may be failing", file=sys.stderr)
    else:
        print(f"SUCCESS: avg {out_df['num_windows'].mean():.1f} windows/pair, "
              f"avg MFE {out_df['avg_window_mfe'].mean():.1f} kcal/mol")

    diff_cols = [c for c in out_df.columns if c.startswith('diff_')]
    if diff_cols:
        try:
            from scipy.stats import pearsonr
            print("\nCorrelation preview (diff_* vs fraction_downstream):")
            for col in diff_cols:
                valid = out_df[['fraction_downstream', col]].dropna()
                if len(valid) > 50:
                    r, p = pearsonr(valid['fraction_downstream'], valid[col])
                    print(f"  {col}: r={r:.3f} p={p:.2e} n={len(valid):,}")
        except ImportError:
            pass


if __name__ == '__main__':
    main()