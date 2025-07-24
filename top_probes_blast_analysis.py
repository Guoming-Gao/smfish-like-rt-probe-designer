#!/usr/bin/env python3
# top_probes_blast_analysis.py - BLAST Analysis for Top Probes (FIXED)

import pandas as pd
import re
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import os
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

    print(f"🔍 Found {len(queries)} queries in BLAST results")

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
            print(f"⚠️  Could not parse gene name from query {i+1}: {query[:100]}...")

        # Count alignment sections starting with ">"
        alignment_headers = re.findall(r"^>([^\n]+)", query, re.MULTILINE)
        num_hits = len(alignment_headers)

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

        # Debug: Print first few examples
        if i < 3:
            print(
                f"  Debug Query {i+1}: gene={gene_name}, seq_len={len(probe_seq) if probe_seq else 0}, hits={num_hits}"
            )

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
    # Create GUI window (hidden)
    root = Tk()
    root.withdraw()  # Hide main window

    print("🧬 BLAST Analysis for Top Probes")
    print("=" * 50)

    # Step 1: Select BLAST result text file
    print("\n📁 Please select your BLAST result text file...")
    blast_file_path = askopenfilename(
        title="Select BLAST result text file",
        filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
    )

    if not blast_file_path:
        print("❌ No BLAST file selected. Exiting.")
        return

    # Read and parse BLAST results
    print(f"📖 Reading BLAST results from: {os.path.basename(blast_file_path)}")
    with open(blast_file_path, "r") as f:
        blast_text = f.read()

    print("🔄 Parsing BLAST results...")
    blast_df = parse_blast_results(blast_text)
    print(f"✅ Parsed {len(blast_df)} probe results")

    # Debug: Check BLAST parsing results
    print(f"🔍 BLAST data sample:")
    print(f"  Genes found: {blast_df['GeneName'].dropna().unique()[:5]}")
    print(
        f"  Sequences (first 3): {[seq[:20] + '...' if seq and len(seq) > 20 else seq for seq in blast_df['ProbeSequence'].dropna().head(3).tolist()]}"
    )

    # Step 2: Select CSV file with top probes
    print("\n📁 Please select your TOP probes CSV file...")
    csv_file_path = askopenfilename(
        title="Select TOP probes CSV file",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
    )

    if not csv_file_path:
        print("❌ No CSV file selected. Exiting.")
        return

    # Read CSV file
    print(f"📖 Reading CSV file from: {os.path.basename(csv_file_path)}")
    csv_df = pd.read_csv(csv_file_path)

    # Check expected columns exist
    expected_cols = ["GeneName", "Seq"]
    missing_cols = [col for col in expected_cols if col not in csv_df.columns]
    if missing_cols:
        print(f"❌ Error: Missing expected columns: {missing_cols}")
        print(f"📋 Available columns: {list(csv_df.columns)}")
        return

    print(
        f"✅ CSV loaded: {len(csv_df)} probes from {csv_df['GeneName'].nunique()} genes"
    )

    # Debug: Check CSV data
    print(f"🔍 CSV data sample:")
    print(f"  Genes found: {csv_df['GeneName'].unique()[:5]}")
    print(
        f"  Sequences (first 3): {[seq[:20] + '...' if len(seq) > 20 else seq for seq in csv_df['Seq'].head(3).tolist()]}"
    )

    # Step 3: Create matching keys for merging
    print("\n🔄 Preparing data for merging...")

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
    print(f"🔍 Key matching debug:")
    print(f"  CSV keys (first 3): {list(csv_keys)[:3]}")
    print(f"  BLAST keys (first 3): {list(blast_keys)[:3]}")
    print(f"  Common keys: {len(csv_keys & blast_keys)}")
    print(f"  CSV only: {len(csv_keys - blast_keys)}")
    print(f"  BLAST only: {len(blast_keys - csv_keys)}")

    # Step 4: Merge dataframes
    print("🔗 Merging BLAST results with probe data...")
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

    # Step 5: Filter to keep only unique hits
    print("🎯 Filtering for unique hits only...")
    unique_df = merged_df[merged_df["NumberOfHits"] == 1].copy()

    # Step 6: Generate summary statistics
    print("\n📊 BLAST Analysis Summary:")
    print(f"📌 Total probes in CSV: {len(csv_df)}")
    print(
        f"📌 Probes with BLAST results: {len(merged_df.dropna(subset=['NumberOfHits']))}"
    )
    print(f"📌 Probes with unique hits: {len(unique_df)}")
    print(f"📌 Probes removed (multiple hits): {len(merged_df) - len(unique_df)}")

    # Hit distribution
    if "NumberOfHits" in merged_df.columns:
        print(f"\n📈 Hit Count Distribution:")
        hit_counts = merged_df["NumberOfHits"].value_counts().sort_index()
        for hits, count in hit_counts.items():
            if pd.notna(hits):
                print(f"   {int(hits)} hit(s): {count} probes")

        no_results = merged_df["NumberOfHits"].isna().sum()
        if no_results > 0:
            print(f"   No BLAST results: {no_results} probes")

    # Gene-level summary (FIXED: Handle division by zero)
    print(f"\n🧬 Gene-level Summary:")
    gene_summary = (
        merged_df.groupby("GeneName")
        .agg({"NumberOfHits": ["count", lambda x: (x == 1).sum()]})
        .round(2)
    )
    gene_summary.columns = ["Total_Probes", "Unique_Hits"]

    for gene, row in gene_summary.iterrows():
        total = int(row["Total_Probes"])
        unique = int(row["Unique_Hits"])
        # FIXED: Handle division by zero
        percentage = unique / total * 100 if total > 0 else 0
        print(f"   {gene}: {unique}/{total} unique hits ({percentage:.1f}%)")

    # Step 7: Save results automatically
    print(f"\n💾 Saving results...")

    # Save BLAST analysis results
    blast_csv_path = re.sub(
        r"\.txt$", "_analysis.csv", blast_file_path, flags=re.IGNORECASE
    )
    blast_df.to_csv(blast_csv_path, index=False)
    print(f"✅ BLAST analysis saved: {os.path.basename(blast_csv_path)}")

    # Save merged results (all probes with BLAST data)
    base_csv_name = os.path.splitext(csv_file_path)[0]
    merged_csv_path = f"{base_csv_name}_with_BLAST.csv"
    merged_df.to_csv(merged_csv_path, index=False)
    print(f"✅ Merged data saved: {os.path.basename(merged_csv_path)}")

    # Save unique hits only
    unique_csv_path = f"{base_csv_name}_UNIQUE_HITS.csv"
    unique_df.to_csv(unique_csv_path, index=False)
    print(f"✅ Unique hits saved: {os.path.basename(unique_csv_path)}")

    # NEW: Save unique hits as FASTA file
    if len(unique_df) > 0:
        unique_fasta_path = f"{base_csv_name}_UNIQUE_HITS.fasta"

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
        print(
            f"✅ Unique hits FASTA saved: {os.path.basename(unique_fasta_path)} ({len(records)} sequences)"
        )
    else:
        print(f"⚠️  No unique hits to export as FASTA")

    # Step 8: Quality assessment
    print(f"\n🎯 Quality Assessment:")

    if len(unique_df) > 0:
        # Percent identity distribution for unique hits
        print(f"📊 Percent Identity (Unique Hits):")
        print(f"   Average: {unique_df['PercentAlignment'].mean():.1f}%")
        print(
            f"   Range: {unique_df['PercentAlignment'].min():.0f}% - {unique_df['PercentAlignment'].max():.0f}%"
        )

        # Perfect matches
        perfect_matches = (unique_df["PercentAlignment"] == 100).sum()
        print(f"   Perfect matches (100%): {perfect_matches} probes")

        # SNP coverage for unique hits
        if "SNPs_Covered_Count" in unique_df.columns:
            avg_snps = unique_df["SNPs_Covered_Count"].mean()
            print(f"📈 Average SNPs covered (unique hits): {avg_snps:.1f}")
    else:
        print(f"⚠️  No unique hits found for quality assessment")

    print(f"\n🎉 Analysis complete!")
    print(f"📁 Output files saved in: {os.path.dirname(csv_file_path)}")

    # Step 9: Recommendations
    print(f"\n💡 Recommendations:")
    if len(unique_df) > 0:
        print(
            f"✅ Use probes from '{os.path.basename(unique_csv_path)}' for experiments"
        )
        print(f"✅ These {len(unique_df)} probes have unique genomic targets")
    else:
        print(f"⚠️  No probes with unique hits found")
        print(f"🔍 Check the debug output above to troubleshoot matching issues")
        print(f"🔧 Possible issues:")
        print(f"   - Gene name extraction from BLAST headers")
        print(f"   - Sequence extraction from BLAST alignments")
        print(f"   - Case sensitivity in gene names")

    root.destroy()


if __name__ == "__main__":
    main()
