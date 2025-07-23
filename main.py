#!/usr/bin/env python3
# main.py - SIMPLIFIED (No BLAST Integration)

import os
import sys
import argparse
import pandas as pd
from pathlib import Path
from rich.progress import track
from rich.console import Console

# Import our modules (BLAST module removed)
from config import FISH_RT_CONFIG, TEST_GENES_21
from gene_fetcher import GeneSequenceFetcher
from utils.oligostan_core import design_fish_probes
from snp_analyzer import SNPCoverageAnalyzer
from output_generator import OutputGenerator

console = Console()


def setup_output_directory(output_path):
    """Create output directory structure"""
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    # Create subdirectories (removed blast_results)
    (output_path / "gene_sequences").mkdir(exist_ok=True)
    (output_path / "snp_analysis").mkdir(exist_ok=True)
    (output_path / "fasta_for_blast").mkdir(exist_ok=True)

    return output_path


def main():
    """Main FISH-RT probe design pipeline (Manual BLAST workflow)"""

    console.print(
        "[bold blue]🧬 smfish-like-rt-probe-designer (Manual BLAST)[/bold blue]"
    )
    console.print("=" * 60)
    console.print(
        "[cyan]Workflow: Design → Filter → Export FASTA → Manual NCBI BLAST[/cyan]"
    )

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Design FISH probes for RT-PCR with manual BLAST workflow"
    )
    parser.add_argument("--config", help="Path to custom config file")
    parser.add_argument("--output", help="Output directory path")
    parser.add_argument("--genes", nargs="+", help="List of gene names to process")
    parser.add_argument("--test", action="store_true", help="Run with test gene set")
    args = parser.parse_args()

    # Load configuration
    config = FISH_RT_CONFIG.copy()

    # Update output directory
    if args.output:
        config["output_directory"] = args.output

    # Determine gene list
    if args.test:
        gene_list = TEST_GENES_21
        console.print(f"[yellow]Using test gene set: {len(gene_list)} genes[/yellow]")
    elif args.genes:
        gene_list = args.genes
        console.print(
            f"[green]Processing user-specified genes: {len(gene_list)} genes[/green]"
        )
    else:
        gene_list = config["gene_list"]
        console.print(f"[blue]Using config gene list: {len(gene_list)} genes[/blue]")

    # Setup output directory
    output_path = setup_output_directory(config["output_directory"])
    console.print(f"[green]Output directory: {output_path}[/green]")

    # Check configuration status (BLAST removed)
    console.print("\n[bold]Configuration Status:[/bold]")
    console.print(
        f"  Local genome FASTA: {'✅ Available' if config['local_genome_fasta_path'] else '❌ Not configured'}"
    )
    console.print(
        f"  SNP file: {'✅ Available' if config['snp_file_path'] else '❌ Not configured'}"
    )
    console.print(f"  Manual BLAST: ✅ FASTA files will be generated")

    # Initialize components (BLAST analyzer removed)
    console.print("\n[bold]Initializing pipeline components...[/bold]")

    gene_fetcher = GeneSequenceFetcher(config)
    snp_analyzer = SNPCoverageAnalyzer(config) if config["snp_file_path"] else None
    output_generator = OutputGenerator(config)

    # Step 1: Fetch gene sequences
    console.print("\n[bold cyan]Step 1: Fetching gene sequences...[/bold cyan]")

    all_gene_data = []
    failed_genes = []

    for gene_name in track(gene_list, description="Fetching sequences..."):
        try:
            gene_data = gene_fetcher.fetch_gene_sequence(gene_name)
            if gene_data:
                all_gene_data.append(gene_data)
                console.print(
                    f"✅ {gene_name}: {len(gene_data['sequence'])} bp ({gene_data['source']})"
                )
            else:
                failed_genes.append(gene_name)
                console.print(f"❌ {gene_name}: Failed to fetch")
        except Exception as e:
            failed_genes.append(gene_name)
            console.print(f"❌ {gene_name}: Error - {str(e)}")

    if failed_genes:
        console.print(
            f"[red]Failed to fetch {len(failed_genes)} genes: {failed_genes}[/red]"
        )

    if not all_gene_data:
        console.print("[red]No genes successfully fetched. Exiting.[/red]")
        return

    # Save gene sequences as FASTA files
    console.print("\n[bold cyan]Saving gene sequences...[/bold cyan]")
    gene_fetcher.save_gene_sequences(all_gene_data, output_path / "gene_sequences")

    # Step 2: Design FISH probes
    console.print("\n[bold cyan]Step 2: Designing FISH probes...[/bold cyan]")

    all_probes = []
    for gene_data in track(all_gene_data, description="Designing probes..."):
        try:
            probes = design_fish_probes(gene_data, config)
            if probes:
                all_probes.extend(probes)
                console.print(f"✅ {gene_data['gene_name']}: {len(probes)} probes")
            else:
                console.print(f"⚠️  {gene_data['gene_name']}: No probes generated")
        except Exception as e:
            console.print(f"❌ {gene_data['gene_name']}: Error - {str(e)}")

    if not all_probes:
        console.print("[red]No probes generated. Exiting.[/red]")
        return

    console.print(f"[green]Total probes generated: {len(all_probes)}[/green]")

    # Step 3: SNP coverage analysis
    if snp_analyzer:
        console.print("\n[bold cyan]Step 3: SNP coverage analysis...[/bold cyan]")
        try:
            all_probes = snp_analyzer.analyze_probes(all_probes)

            # Report SNP coverage statistics
            snp_counts = [p.get("SNPs_Covered_Count", 0) for p in all_probes]
            avg_snps = sum(snp_counts) / len(snp_counts) if snp_counts else 0
            max_snps = max(snp_counts) if snp_counts else 0
            probes_with_snps = sum(1 for c in snp_counts if c > 0)

            console.print(f"✅ SNP coverage analysis completed")
            console.print(f"   Average SNPs per probe: {avg_snps:.1f}")
            console.print(f"   Maximum SNPs per probe: {max_snps}")
            console.print(
                f"   Probes with SNPs: {probes_with_snps}/{len(all_probes)} ({probes_with_snps/len(all_probes)*100:.1f}%)"
            )

        except Exception as e:
            console.print(f"❌ SNP analysis failed: {str(e)}")
    else:
        console.print(
            "\n[yellow]Step 3: SNP analysis skipped (no SNP file configured)[/yellow]"
        )

    # Step 4: Generate outputs with FASTA for manual BLAST
    console.print(
        "\n[bold cyan]Step 4: Generating outputs + BLAST FASTA files...[/bold cyan]"
    )

    try:
        output_files = output_generator.generate_outputs(all_probes, output_path)

        console.print("[bold green]✅ Pipeline completed successfully![/bold green]")
        console.print(f"[green]Generated files:[/green]")
        for file_path in output_files:
            console.print(f"  📄 {file_path}")

        # Print summary statistics
        csv_files = [f for f in output_files if str(f).endswith(".csv")]
        if csv_files:
            df_all = pd.read_csv(csv_files[0])  # First CSV file
            df_filt = pd.read_csv(csv_files[1]) if len(csv_files) > 1 else df_all

            console.print(f"\n[bold]📊 Final Summary:[/bold]")
            console.print(f"  Total probes designed: {len(df_all)}")
            console.print(f"  Probes after quality filtering: {len(df_filt)}")
            console.print(f"  Genes processed: {len(set(df_all['GeneName']))}")
            console.print(f"  Output directory: {output_path}")

            if "SNPs_Covered_Count" in df_all.columns:
                avg_snps = df_all["SNPs_Covered_Count"].mean()
                console.print(f"  Average SNPs covered per probe: {avg_snps:.1f}")

    except Exception as e:
        console.print(f"[red]❌ Output generation failed: {str(e)}[/red]")
        return


if __name__ == "__main__":
    main()
