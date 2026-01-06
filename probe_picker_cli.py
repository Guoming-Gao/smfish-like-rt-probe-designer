#!/usr/bin/env python3
# probe_picker_cli.py - CLI-based Top Probe Selection Tool

import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def select_top_probes(df, default_probes_per_gene, probes_per_gene_override=None):
    """
    Sort and select top N probes per gene.

    Args:
        df: DataFrame with probe data
        default_probes_per_gene: Default number of probes to select per gene
        probes_per_gene_override: Dict mapping gene names to specific counts
                                  e.g., {"Xist": 20, "Tsix": 10}
    """
    if probes_per_gene_override is None:
        probes_per_gene_override = {}

    # Sort probes by selection criteria
    df_sorted = df.sort_values(
        ["SNP_Count", "NbOfPNAS", "dGScore"],
        ascending=[False, False, False],
    )

    # Select top probes per gene with per-gene override support
    selected_probes = []
    for gene_name, gene_df in df_sorted.groupby("GeneName"):
        n_probes = probes_per_gene_override.get(gene_name, default_probes_per_gene)
        selected_probes.append(gene_df.head(n_probes))

    top_probes = pd.concat(selected_probes, ignore_index=True)

    # Sort final output by GeneName, then SNP_Count (desc), then NbOfPNAS (desc)
    top_probes = top_probes.sort_values(
        ["GeneName", "SNP_Count", "NbOfPNAS"],
        ascending=[True, False, False],
    )

    return top_probes


def generate_fasta(df, output_file):
    """Generate FASTA file from CSV data"""
    records = []

    for idx, row in df.iterrows():
        gene_name = row["GeneName"]

        # Use probe sequence only (no RTBC for BLAST)
        sequence = row["Probe_Seq"]
        seq_type = "probe_only"

        # Create descriptive FASTA header
        probe_size = row.get("Probe_Length", len(row["Probe_Seq"]))
        snp_count = row.get("SNP_Count", 0)
        region_type = row.get("RegionType", "unknown")
        pnas_score = row.get("NbOfPNAS", 0)

        description = f"{gene_name} | {seq_type} | size:{probe_size}nt | SNPs:{snp_count} | region:{region_type} | PNAS:{pnas_score}/5"

        record = SeqRecord(
            Seq(sequence), id=f"{gene_name}_probe_{idx}", description=description
        )
        records.append(record)

    SeqIO.write(records, output_file, "fasta")
    return len(records)


def main():
    parser = argparse.ArgumentParser(
        description="Select top N probes per gene from FISH-RT probe CSV file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Selection Criteria (in order):
  1. SNP_Count (descending) - Higher SNP coverage first
  2. NbOfPNAS (descending) - Better PNAS score first
  3. dGScore (descending) - Better thermodynamic score first

Example:
  python probe_picker_cli.py -i FISH_RT_probes_FINAL_SELECTION.csv -n 10
  python probe_picker_cli.py -i candidates.csv -o MY_TOP_PROBES -n 5
"""
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input CSV file with probe data"
    )
    parser.add_argument(
        "-o", "--output-prefix",
        default="FISH_RT_probes_TOP10",
        help="Output file prefix (default: FISH_RT_probes_TOP10)"
    )
    parser.add_argument(
        "-n", "--probes-per-gene",
        type=int,
        default=10,
        help="Number of probes to select per gene (default: 10)"
    )

    args = parser.parse_args()

    # Validate input file
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {args.input}")
        return 1

    print(f"FISH-RT Probe Picker (CLI)")
    print("=" * 50)

    # Read CSV file
    print(f"\nReading CSV file: {input_path.name}")
    df = pd.read_csv(input_path)

    print(f"  Total probes: {len(df)}")
    print(f"  Genes found: {df['GeneName'].nunique()}")

    # Select top probes
    print(f"\nSorting probes by SNP coverage, PNAS score, and dG score...")
    print(f"Selecting top {args.probes_per_gene} probes per gene...")
    top_probes = select_top_probes(df, args.probes_per_gene)

    # Report selection statistics
    print(f"\nSelected probes per gene:")
    gene_counts = top_probes["GeneName"].value_counts().sort_index()
    for gene, count in gene_counts.items():
        avg_snps = top_probes[top_probes["GeneName"] == gene]["SNP_Count"].mean()
        avg_pnas = top_probes[top_probes["GeneName"] == gene]["NbOfPNAS"].mean()
        print(f"  {gene}: {count} probes (avg SNPs: {avg_snps:.1f}, avg PNAS: {avg_pnas:.1f}/5)")

    print(f"\nTotal selected: {len(top_probes)} probes from {top_probes['GeneName'].nunique()} genes")

    # Save outputs
    output_dir = input_path.parent
    csv_output = output_dir / f"{args.output_prefix}.csv"
    fasta_output = output_dir / f"{args.output_prefix}.fasta"

    top_probes.to_csv(csv_output, index=False)
    print(f"\nSaved CSV: {csv_output}")

    num_seqs = generate_fasta(top_probes, fasta_output)
    print(f"Saved FASTA: {fasta_output} ({num_seqs} sequences)")

    # Summary statistics
    print(f"\nSUMMARY STATISTICS:")
    print(f"  Input probes: {len(df)}")
    print(f"  Selected probes: {len(top_probes)}")
    print(f"  Reduction: {(1 - len(top_probes)/len(df))*100:.1f}%")
    print(f"  Average SNPs per probe: {top_probes['SNP_Count'].mean():.1f}")
    print(f"  Average PNAS score: {top_probes['NbOfPNAS'].mean():.1f}/5")

    # SNP coverage distribution
    snp_dist = top_probes["SNP_Count"].value_counts().sort_index()
    print(f"\n  SNP coverage distribution:")
    for snps, count in snp_dist.items():
        print(f"    {snps} SNPs: {count} probes")

    print(f"\nProbe selection completed!")
    return 0


if __name__ == "__main__":
    exit(main())
