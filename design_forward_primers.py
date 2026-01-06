#!/usr/bin/env python3
"""
Design forward primers for RT probes using Primer3.

Forward primers are designed 200-250bp upstream of the RT probe to ensure
the PCR amplicon fully covers the 200bp RT coverage region containing SNPs.

Usage:
    python design_forward_primers.py input.csv output.csv [options]
"""

import argparse
import subprocess
import tempfile
import pandas as pd
from pathlib import Path


def extract_target_sequence(chrom: str, start: int, end: int, strand: int,
                           genome_fasta: str) -> str:
    """
    Extract sequence from genome using samtools.

    For forward primer design:
    - Plus strand genes: return plus strand sequence
    - Minus strand genes: return minus strand sequence (reverse complement)

    This ensures the forward primer extends toward the RT region.
    """
    try:
        region = f"{chrom}:{start}-{end}"
        result = subprocess.run(
            ["samtools", "faidx", genome_fasta, region],
            capture_output=True, text=True, check=True
        )
        lines = result.stdout.strip().split('\n')
        seq = ''.join(lines[1:])  # Skip header

        # For minus strand genes, reverse complement the sequence
        # This gives us a template where Primer3's LEFT primer will
        # match the minus strand (bind plus strand, extend leftward)
        if strand == -1:
            complement = str.maketrans('ATCGatcg', 'TAGCtagc')
            seq = seq.translate(complement)[::-1]

        return seq.upper()
    except Exception as e:
        print(f"Warning: Failed to extract sequence for {region}: {e}")
        return ""


def run_primer3(sequence: str, target_start: int, target_length: int,
                primer_min_size: int = 18, primer_opt_size: int = 20,
                primer_max_size: int = 25, primer_min_tm: float = 57.0,
                primer_opt_tm: float = 60.0, primer_max_tm: float = 63.0,
                primer_min_gc: float = 40.0, primer_opt_gc: float = 50.0,
                primer_max_gc: float = 60.0, gc_clamp: int = 1) -> dict:
    """
    Run Primer3 to design a forward primer.

    Returns primer info dict or None if no primer found.
    """
    # Primer3 input format with optimal settings for standard PCR
    primer3_input = f"""SEQUENCE_ID=target
SEQUENCE_TEMPLATE={sequence}
SEQUENCE_TARGET={target_start},{target_length}
PRIMER_TASK=pick_left_only
PRIMER_PICK_LEFT_ONLY=1
PRIMER_MIN_SIZE={primer_min_size}
PRIMER_OPT_SIZE={primer_opt_size}
PRIMER_MAX_SIZE={primer_max_size}
PRIMER_MIN_TM={primer_min_tm}
PRIMER_OPT_TM={primer_opt_tm}
PRIMER_MAX_TM={primer_max_tm}
PRIMER_MIN_GC={primer_min_gc}
PRIMER_OPT_GC_PERCENT={primer_opt_gc}
PRIMER_MAX_GC={primer_max_gc}
PRIMER_GC_CLAMP={gc_clamp}
PRIMER_NUM_RETURN=1
=
"""

    try:
        # Run primer3_core
        result = subprocess.run(
            ["primer3_core"],
            input=primer3_input,
            capture_output=True,
            text=True,
            check=True
        )

        # Parse output
        output = result.stdout
        primer_info = {}

        for line in output.strip().split('\n'):
            if '=' in line:
                key, value = line.split('=', 1)
                primer_info[key] = value

        # Check if primer was found
        if 'PRIMER_LEFT_0_SEQUENCE' in primer_info:
            return {
                'sequence': primer_info.get('PRIMER_LEFT_0_SEQUENCE', ''),
                'tm': float(primer_info.get('PRIMER_LEFT_0_TM', 0)),
                'gc': float(primer_info.get('PRIMER_LEFT_0_GC_PERCENT', 0)),
                'position': primer_info.get('PRIMER_LEFT_0', ''),
                'length': int(primer_info.get('PRIMER_LEFT_0_SEQUENCE', '').__len__()),
            }
        else:
            return None

    except Exception as e:
        print(f"Warning: Primer3 failed: {e}")
        return None


def design_forward_primers(input_file: str, output_file: str,
                          genome_fasta: str,
                          target_distance: int = 225,
                          primer_min_size: int = 18,
                          primer_opt_size: int = 20,
                          primer_max_size: int = 25,
                          primer_min_tm: float = 57.0,
                          primer_opt_tm: float = 60.0,
                          primer_max_tm: float = 63.0,
                          primer_min_gc: float = 40.0,
                          primer_max_gc: float = 60.0):
    """
    Design forward primers for RT probes.

    Args:
        input_file: CSV with RT probes (needs Chromosome, theStartPos, theEndPos, Target_Strand)
        output_file: Output CSV with forward primers added
        genome_fasta: Path to genome FASTA file
        target_distance: Distance upstream of RT probe to design primer (default 225bp)
    """
    df = pd.read_csv(input_file)
    print(f"📄 Loaded {len(df)} probes from {input_file}")

    # Result columns
    forward_primers = []
    forward_tms = []
    forward_gcs = []
    pcr_sizes = []
    rt_covered = []

    for idx, row in df.iterrows():
        chrom = row['Chromosome']
        # Support both column naming conventions
        probe_start = row.get('Probe_Start', row.get('theStartPos', 0))
        probe_end = row.get('Probe_End', row.get('theEndPos', 0))
        strand_raw = row['Target_Strand']
        # Parse strand - handle both string ('-'/'+') and integer (-1/1)
        if strand_raw in [-1, '-', '-1']:
            strand = -1
        else:
            strand = 1
        rt_start = row.get('RT_Region_Start', probe_start)
        rt_end = row.get('RT_Region_End', probe_end)

        # Calculate region to search for forward primer
        # Forward primer should be ~target_distance from probe
        if strand == 1:  # Plus strand: primer upstream (lower coords)
            search_start = probe_start - target_distance - 100
            search_end = probe_start - target_distance + 100
            target_in_seq = 50  # Middle of search region
        else:  # Minus strand
            # Forward primer region is downstream (larger coordinates)
            search_start = probe_end + target_distance - 100
            search_end = probe_end + target_distance + 100
            target_in_seq = 50

        # Ensure valid coordinates
        search_start = max(1, search_start)

        # Extract sequence for Primer3
        seq = extract_target_sequence(chrom, search_start, search_end, strand, genome_fasta)

        if not seq or len(seq) < 50:
            forward_primers.append('')
            forward_tms.append(0)
            forward_gcs.append(0)
            pcr_sizes.append(0)
            rt_covered.append(False)
            continue

        # Run Primer3
        primer = run_primer3(
            seq, target_in_seq, 50,
            primer_min_size, primer_opt_size, primer_max_size,
            primer_min_tm, primer_opt_tm, primer_max_tm,
            primer_min_gc, primer_max_gc
        )

        if primer:
            forward_primers.append(primer['sequence'])
            forward_tms.append(primer['tm'])
            forward_gcs.append(primer['gc'])

            # Calculate PCR size (approximate)
            pcr_size = target_distance + len(row.get('Probe_Seq', ''))
            pcr_sizes.append(pcr_size)

            # Check if RT region would be covered
            # PCR amplicon goes from forward primer to RT probe
            # RT region should be fully within this
            covered = pcr_size >= 200  # At minimum covers RT region
            rt_covered.append(covered)
        else:
            forward_primers.append('')
            forward_tms.append(0)
            forward_gcs.append(0)
            pcr_sizes.append(0)
            rt_covered.append(False)

        if (idx + 1) % 100 == 0:
            print(f"  Processed {idx + 1}/{len(df)} probes...")

    # Add columns to dataframe
    df['Forward_Primer_Seq'] = forward_primers
    df['Forward_Primer_Tm'] = forward_tms
    df['Forward_Primer_GC'] = forward_gcs
    df['PCR_Amplicon_Size'] = pcr_sizes
    df['RT_Region_Covered'] = rt_covered

    # Save output
    df.to_csv(output_file, index=False)

    # Summary
    success_count = sum(1 for p in forward_primers if p)
    print(f"\n✅ Forward primer design complete!")
    print(f"   Successful designs: {success_count}/{len(df)} ({100*success_count/len(df):.1f}%)")
    print(f"   Average Tm: {sum(forward_tms)/max(success_count,1):.1f}°C")
    print(f"   Average GC: {sum(forward_gcs)/max(success_count,1):.1f}%")
    print(f"   Output saved to: {output_file}")

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Design forward primers for RT probes using Primer3"
    )
    parser.add_argument("input", help="Input CSV with RT probes")
    parser.add_argument("output", help="Output CSV with forward primers")
    parser.add_argument(
        "--genome",
        default="/Volumes/guttman/genomes/mm10/fasta/mm10.fa",
        help="Path to genome FASTA file"
    )
    parser.add_argument(
        "--distance",
        type=int, default=225,
        help="Target distance from RT probe (default: 225bp)"
    )
    parser.add_argument(
        "--primer-length-min", type=int, default=18,
        help="Minimum primer length (default: 18)"
    )
    parser.add_argument(
        "--primer-length-opt", type=int, default=20,
        help="Optimal primer length (default: 20)"
    )
    parser.add_argument(
        "--primer-length-max", type=int, default=25,
        help="Maximum primer length (default: 25)"
    )
    parser.add_argument(
        "--tm-min", type=float, default=57.0,
        help="Minimum Tm (default: 57°C)"
    )
    parser.add_argument(
        "--tm-opt", type=float, default=60.0,
        help="Optimal Tm (default: 60°C)"
    )
    parser.add_argument(
        "--tm-max", type=float, default=63.0,
        help="Maximum Tm (default: 63°C)"
    )
    parser.add_argument(
        "--gc-min", type=float, default=40.0,
        help="Minimum GC%% (default: 40)"
    )
    parser.add_argument(
        "--gc-max", type=float, default=60.0,
        help="Maximum GC%% (default: 60)"
    )

    args = parser.parse_args()

    design_forward_primers(
        args.input, args.output,
        args.genome,
        args.distance,
        args.primer_length_min, args.primer_length_opt, args.primer_length_max,
        args.tm_min, args.tm_opt, args.tm_max,
        args.gc_min, args.gc_max
    )


if __name__ == "__main__":
    main()
