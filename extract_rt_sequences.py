#!/usr/bin/env python3

import pandas as pd
import subprocess
from pathlib import Path
from Bio.Seq import Seq

def extract_sequence_samtools(chromosome, start, end, strand, genome_fasta_path):
    """Extract sequence using samtools, handling strand properly"""
    try:
        region = f"{chromosome}:{start}-{end}"
        cmd = ["samtools", "faidx", genome_fasta_path, region]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        lines = result.stdout.strip().split('\n')
        if len(lines) < 2:
            return None

        genomic_sequence = ''.join(lines[1:]).upper()

        # Handle strand orientation
        if strand == '+':
            return genomic_sequence
        else:  # strand == '-'
            return str(Seq(genomic_sequence).reverse_complement())

    except Exception as e:
        print(f"Sequence extraction error: {e}")
        return None

def extract_rt_sequences(input_file, output_file):
    """Extract RT coverage sequences and validate target sequences"""

    # Genome path
    genome_fasta_path = "/home/gmgao/Desktop/Central_Guttman/genomes/mm10/GRCm3868.fa"

    print(f"Reading input file: {input_file}")

    # Read input file
    if input_file.endswith('.xlsx'):
        df = pd.read_excel(input_file)
    else:
        df = pd.read_csv(input_file)

    print(f"Loaded {len(df)} probes")

    # Initialize RT_seq column
    df['RT_seq'] = ''

    # Process each probe
    for idx, row in df.iterrows():
        print(f"Processing {idx+1}/{len(df)}: {row['GeneName']}")

        # Extract target sequence for validation
        target_seq = extract_sequence_samtools(
            row['Chromosome'], 
            row['theStartPos'], 
            row['theEndPos'], 
            row['Target_Strand'],
            genome_fasta_path
        )

        # Validate target sequence
        if target_seq and target_seq == row['Target_Seq']:
            print(f"  ✓ Target sequence validated")
        else:
            print(f"  ⚠ Target sequence mismatch for {row['GeneName']}")
            print(f"    Expected: {row['Target_Seq']}")
            print(f"    Got:      {target_seq}")

        # Extract RT coverage sequence
        rt_seq = extract_sequence_samtools(
            row['Chromosome'],
            row['RT_Coverage_Start'],
            row['RT_Coverage_End'],
            row['RT_Coverage_Strand'],
            genome_fasta_path
        )

        if rt_seq:
            df.at[idx, 'RT_seq'] = rt_seq
            print(f"  ✓ RT sequence extracted ({len(rt_seq)} bp)")
        else:
            print(f"  ✗ Failed to extract RT sequence")

    # Select the specified columns plus RT_seq
    output_columns = [
        'GeneName',
        'Species', 
        'Chromosome',
        'RegionType',
        'Seq',
        'Target_Seq',
        'ProbeSize',
        'theStartPos',
        'theEndPos',
        'RT_Coverage_Start',
        'RT_Coverage_End',
        'RT_seq'
    ]

    # Create output dataframe
    output_df = df[output_columns].copy()

    # Save to CSV
    output_df.to_csv(output_file, index=False)
    print(f"\nOutput saved to: {output_file}")
    print(f"Total probes processed: {len(output_df)}")
    print(f"RT sequences extracted: {len(output_df[output_df['RT_seq'] != ''])}")

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage: python extract_rt_sequences.py <input_file> <output_file>")
        print("Example: python extract_rt_sequences.py probes.xlsx probes_with_rt.csv")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    extract_rt_sequences(input_file, output_file)
