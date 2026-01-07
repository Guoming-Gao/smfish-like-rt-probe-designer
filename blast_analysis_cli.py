#!/usr/bin/env python3
# blast_analysis_cli.py - CLI-based BLAST Analysis Tool

import argparse
import pandas as pd
import re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_blast_results(blast_text):
    """
    Parse BLAST result text for top probes and extract required information
    Returns DataFrame with columns: ProbeName, ProbeSequence, PercentAlignment,
    NumberOfHits, UniqueHitName, Start, End, GeneName
    """
    # Split by queries - handle the new format
    queries = re.split(r"Query #\d+:\s*", blast_text)[1:]  # skip first empty split
    records = []

    print(f"  Found {len(queries)} queries in BLAST results")

    for i, query in enumerate(queries):
        # Extract probe name from new format: "Atrx_probe_24 Atrx | probe_only | size:27nt..."
        probe_name_search = re.search(r"^(\S+)\s+(\w+)\s+\|", query)
        if probe_name_search:
            full_probe_name = probe_name_search.group(1)  # e.g., "Atrx_probe_24"
            gene_name = probe_name_search.group(2)  # e.g., "Atrx"
        else:
            # Fallback to extract just the first word
            probe_name_search = re.search(r"^(\S+)", query)
            full_probe_name = (
                probe_name_search.group(1) if probe_name_search else f"Query_{i+1}"
            )
            gene_name = None

        # Count individual High Scoring Pairs (HSPs) instead of just subjects
        num_hits = len(re.findall(r" Score =", query))

        # Extract unique hit name if only one hit
        unique_hit_name = None
        if num_hits == 1:
            unique_hit_name = alignment_headers[0].strip()

        # Extract start and end positions (from first alignment)
        start_end_search = re.search(r"Range 1: (\d+) to (\d+)", query)
        start = int(start_end_search.group(1)) if start_end_search else None
        end = int(start_end_search.group(2)) if start_end_search else None

        # Extract percentage identity
        perc_identity_search = re.search(r"Identities:\s*\d+/\d+\((\d+)%\)", query)
        perc_identity = (
            int(perc_identity_search.group(1)) if perc_identity_search else None
        )

        # Extract probe sequence (longest query sequence from alignment)
        query_seqs = re.findall(r"Query\s+\d+\s+([A-Z]+)\s+\d+", query)
        probe_seq = max(query_seqs, key=len) if query_seqs else None

        records.append(
            {
                "ProbeName": full_probe_name,
                "GeneName": gene_name,
                "ProbeSequence": probe_seq,
                "PercentAlignment": perc_identity,
                "NumberOfHits": num_hits,
                "UniqueHitName": unique_hit_name,
                "Start": start,
                "End": end,
            }
        )

    return pd.DataFrame(records)


def create_probe_identifier(row):
    """Create a probe identifier from CSV row data for matching"""
    # Use gene name + sequence as identifier since we don't have probe numbers
    gene_name = str(row.get("GeneName", "")).strip()
    sequence = str(row.get("Seq", "")).strip()
    return f"{gene_name}_{sequence}"


def main():
    parser = argparse.ArgumentParser(
        description="Analyze BLAST results and filter for unique hits",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Output Files:
  *_analysis.csv      - Parsed BLAST results
  *_with_BLAST.csv    - Probe data merged with BLAST results
  *_UNIQUE_HITS.csv   - Probes with exactly one BLAST hit
  *_UNIQUE_HITS.fasta - FASTA of unique hits

Example:
  python blast_analysis_cli.py -b blast_results.txt -c FISH_RT_probes_TOP10.csv
  python blast_analysis_cli.py -b blast.txt -c probes.csv -o ./output
"""
    )

    parser.add_argument(
        "-b", "--blast",
        required=True,
        help="BLAST result text file"
    )
    parser.add_argument(
        "-c", "--csv",
        required=True,
        help="Probes CSV file (e.g., TOP10.csv)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        default=None,
        help="Output directory (default: same as CSV file)"
    )

    args = parser.parse_args()

    # Validate input files
    blast_path = Path(args.blast)
    csv_path = Path(args.csv)

    if not blast_path.exists():
        print(f"Error: BLAST file not found: {args.blast}")
        return 1
    if not csv_path.exists():
        print(f"Error: CSV file not found: {args.csv}")
        return 1

    output_dir = Path(args.output_dir) if args.output_dir else csv_path.parent

    print("BLAST Analysis for Top Probes (CLI)")
    print("=" * 50)

    # Read and parse BLAST results
    print(f"\nReading BLAST results: {blast_path.name}")
    with open(blast_path, "r") as f:
        blast_text = f.read()

    print("Parsing BLAST results...")
    blast_df = parse_blast_results(blast_text)
    print(f"  Parsed {len(blast_df)} probe results")

    # Read CSV file
    print(f"\nReading CSV file: {csv_path.name}")
    csv_df = pd.read_csv(csv_path)

    # Check expected columns exist
    expected_cols = ["GeneName", "Seq"]
    missing_cols = [col for col in expected_cols if col not in csv_df.columns]
    if missing_cols:
        print(f"Error: Missing expected columns: {missing_cols}")
        print(f"Available columns: {list(csv_df.columns)}")
        return 1

    print(f"  Loaded {len(csv_df)} probes from {csv_df['GeneName'].nunique()} genes")

    # Create matching keys for merging
    print("\nPreparing data for merging...")

    # Create probe identifiers for CSV data
    csv_df["ProbeKey"] = csv_df.apply(create_probe_identifier, axis=1)

    # Create matching keys for BLAST data using sequence
    blast_df["ProbeKey"] = blast_df.apply(
        lambda row: (
            f"{row['GeneName']}_{row['ProbeSequence']}"
            if pd.notna(row["GeneName"]) and pd.notna(row["ProbeSequence"])
            else None
        ),
        axis=1,
    )

    # Debug: Check key matching
    csv_keys = set(csv_df["ProbeKey"].dropna())
    blast_keys = set(blast_df["ProbeKey"].dropna())
    common_keys = len(csv_keys & blast_keys)
    print(f"  Common keys found: {common_keys}")

    if common_keys == 0:
        print("\nWarning: No matching keys found between BLAST and CSV data")
        print("  Check that gene names and sequences match between files")

    # Merge dataframes
    print("Merging BLAST results with probe data...")
    merged_df = csv_df.merge(
        blast_df[
            [
                "ProbeKey",
                "ProbeName",
                "PercentAlignment",
                "NumberOfHits",
                "UniqueHitName",
                "Start",
                "End",
            ]
        ],
        on="ProbeKey",
        how="left",
    )

    # Drop the temporary ProbeKey column
    merged_df.drop(columns=["ProbeKey"], inplace=True)

    # Filter to keep only unique hits
    print("Filtering for unique hits only...")
    unique_df = merged_df[merged_df["NumberOfHits"] == 1].copy()

    # Generate summary statistics
    print(f"\nBLAST Analysis Summary:")
    print(f"  Total probes in CSV: {len(csv_df)}")
    print(f"  Probes with BLAST results: {len(merged_df.dropna(subset=['NumberOfHits']))}")
    print(f"  Probes with unique hits: {len(unique_df)}")
    print(f"  Probes removed (multiple hits): {len(merged_df) - len(unique_df)}")

    # Hit distribution
    if "NumberOfHits" in merged_df.columns:
        print(f"\nHit Count Distribution:")
        hit_counts = merged_df["NumberOfHits"].value_counts().sort_index()
        for hits, count in hit_counts.items():
            if pd.notna(hits):
                print(f"  {int(hits)} hit(s): {count} probes")

        no_results = merged_df["NumberOfHits"].isna().sum()
        if no_results > 0:
            print(f"  No BLAST results: {no_results} probes")

    # Gene-level summary
    print(f"\nGene-level Summary:")
    gene_summary = (
        merged_df.groupby("GeneName")
        .agg({"NumberOfHits": ["count", lambda x: (x == 1).sum()]})
        .round(2)
    )
    gene_summary.columns = ["Total_Probes", "Unique_Hits"]

    for gene, row in gene_summary.iterrows():
        total = int(row["Total_Probes"])
        unique = int(row["Unique_Hits"])
        percentage = unique / total * 100 if total > 0 else 0
        print(f"  {gene}: {unique}/{total} unique hits ({percentage:.1f}%)")

    # Save results
    print(f"\nSaving results to: {output_dir}")

    # Save BLAST analysis results
    blast_csv_path = output_dir / f"{blast_path.stem}_analysis.csv"
    blast_df.to_csv(blast_csv_path, index=False)
    print(f"  BLAST analysis: {blast_csv_path.name}")

    # Save merged results (all probes with BLAST data)
    merged_csv_path = output_dir / f"{csv_path.stem}_with_BLAST.csv"
    merged_df.to_csv(merged_csv_path, index=False)
    print(f"  Merged data: {merged_csv_path.name}")

    # Save unique hits only
    unique_csv_path = output_dir / f"{csv_path.stem}_UNIQUE_HITS.csv"
    unique_df.to_csv(unique_csv_path, index=False)
    print(f"  Unique hits: {unique_csv_path.name}")

    # Save unique hits as FASTA file
    if len(unique_df) > 0:
        unique_fasta_path = output_dir / f"{csv_path.stem}_UNIQUE_HITS.fasta"

        records = []
        for idx, row in unique_df.iterrows():
            gene_name = row.get("GeneName", "Unknown")
            sequence = row.get("Seq", "")
            snp_count = row.get("SNPs_Covered_Count", 0)
            pnas_score = row.get("NbOfPNAS", 0)
            percent_identity = row.get("PercentAlignment", 0)

            # Create descriptive FASTA header
            header = f"{gene_name}_unique_{idx}"
            description = f"Gene:{gene_name} | SNPs:{snp_count} | PNAS:{pnas_score} | Identity:{percent_identity}% | BLAST_unique"

            record = SeqRecord(Seq(sequence), id=header, description=description)
            records.append(record)

        # Write FASTA file
        SeqIO.write(records, unique_fasta_path, "fasta")
        print(f"  Unique hits FASTA: {unique_fasta_path.name} ({len(records)} sequences)")
    else:
        print("  No unique hits to export as FASTA")

    # Quality assessment
    print(f"\nQuality Assessment:")

    if len(unique_df) > 0:
        # Percent identity distribution for unique hits
        print(f"  Percent Identity (Unique Hits):")
        print(f"    Average: {unique_df['PercentAlignment'].mean():.1f}%")
        print(f"    Range: {unique_df['PercentAlignment'].min():.0f}% - {unique_df['PercentAlignment'].max():.0f}%")

        # Perfect matches
        perfect_matches = (unique_df["PercentAlignment"] == 100).sum()
        print(f"    Perfect matches (100%): {perfect_matches} probes")

        # SNP coverage for unique hits
        if "SNPs_Covered_Count" in unique_df.columns:
            avg_snps = unique_df["SNPs_Covered_Count"].mean()
            print(f"  Average SNPs covered (unique hits): {avg_snps:.1f}")
    else:
        print("  No unique hits found for quality assessment")

    print(f"\nAnalysis complete!")

    # Recommendations
    print(f"\nRecommendations:")
    if len(unique_df) > 0:
        print(f"  Use probes from '{unique_csv_path.name}' for experiments")
        print(f"  These {len(unique_df)} probes have unique genomic targets")
    else:
        print("  No probes with unique hits found")
        print("  Check the debug output above to troubleshoot matching issues")

    return 0


if __name__ == "__main__":
    exit(main())
