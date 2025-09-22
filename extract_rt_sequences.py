#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


def load_genome(genome_fasta_path):
    """Load genome sequences into memory"""
    print("Loading genome sequences...")
    genome = {}
    for record in SeqIO.parse(genome_fasta_path, "fasta"):
        # Handle both "chr1" and "1" chromosome naming
        chrom_name = record.id
        if chrom_name.startswith("chr"):
            chrom_name = chrom_name[3:]  # Remove "chr" prefix
        genome[chrom_name] = str(record.seq).upper()
    print(f"Loaded {len(genome)} chromosomes")
    return genome


def extract_sequence(chromosome, start, end, strand, genome):
    """Extract sequence from genome dictionary"""
    try:
        if chromosome not in genome:
            print(f"  ✗ Chromosome {chromosome} not found in genome")
            return None

        chrom_seq = genome[chromosome]

        # Convert to 0-based indexing for Python
        seq_start = start - 1
        seq_end = end

        if seq_start < 0 or seq_end > len(chrom_seq):
            print(f"  ✗ Coordinates out of bounds: {start}-{end}")
            return None

        genomic_sequence = chrom_seq[seq_start:seq_end]

        # Handle strand orientation
        if strand == "+":
            return genomic_sequence
        else:  # strand == '-'
            return str(Seq(genomic_sequence).reverse_complement())

    except Exception as e:
        print(f"  ✗ Sequence extraction error: {e}")
        return None


def extract_rt_sequences(input_file, output_file):
    """Extract RT coverage sequences and validate target sequences"""

    # Genome path
    genome_fasta_path = "/home/gmgao/Desktop/Central_Guttman/genomes/mm10/GRCm38_68.fa"

    print(f"Reading input file: {input_file}")

    # Read input file
    if input_file.endswith(".xlsx"):
        df = pd.read_excel(input_file)
    else:
        df = pd.read_csv(input_file)

    print(f"Loaded {len(df)} probes")

    # Load genome sequences
    genome = load_genome(genome_fasta_path)

    # Initialize RT_seq column
    df["RT_seq"] = ""

    # Process each probe
    for idx, row in df.iterrows():
        print(f"Processing {idx+1}/{len(df)}: {row['GeneName']}")

        # Extract target sequence for validation
        target_seq = extract_sequence(
            row["Chromosome"],
            row["theStartPos"],
            row["theEndPos"],
            row["Target_Strand"],
            genome,
        )

        # Validate target sequence
        if target_seq and target_seq == row["Target_Seq"]:
            print(f"  ✓ Target sequence validated")
        else:
            print(f"  ⚠ Target sequence mismatch for {row['GeneName']}")
            if target_seq:
                print(f"    Expected: {row['Target_Seq']}")
                print(f"    Got:      {target_seq}")

        # Extract RT coverage sequence
        rt_seq = extract_sequence(
            row["Chromosome"],
            row["RT_Coverage_Start"],
            row["RT_Coverage_End"],
            row["RT_Coverage_Strand"],
            genome,
        )

        if rt_seq:
            df.at[idx, "RT_seq"] = rt_seq
            print(f"  ✓ RT sequence extracted ({len(rt_seq)} bp)")
        else:
            print(f"  ✗ Failed to extract RT sequence")

    # Select the specified columns plus RT_seq
    output_columns = [
        "GeneName",
        "Species",
        "Chromosome",
        "RegionType",
        "Seq",
        "Target_Seq",
        "ProbeSize",
        "theStartPos",
        "theEndPos",
        "RT_Coverage_Start",
        "RT_Coverage_End",
        "RT_seq",
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
