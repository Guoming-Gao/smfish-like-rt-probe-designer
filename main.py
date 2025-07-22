#!/usr/bin/env python3
# main.py - smfish-like-rt-probe-designer Main Pipeline

import os
import sys
import argparse
import pandas as pd
from pathlib import Path
from rich.progress import track
from rich.console import Console

# Import our modules
from config import FISH_RT_CONFIG, TEST_GENES_21
from gene_fetcher import GeneSequenceFetcher
from utils.oligostan_core import design_fish_probes
from snp_analyzer import SNPCoverageAnalyzer
from blast_analyzer import BLASTSpecificityAnalyzer
from output_generator import OutputGenerator

console = Console()


def setup_output_directory(output_path):
    """Create output directory structure"""
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    # Create subdirectories
    (output_path / "gene_sequences").mkdir(exist_ok=True)
    (output_path / "blast_results").mkdir(exist_ok=True)
    (output_path / "snp_analysis").mkdir(exist_ok=True)

    return output_path


def main():
    """Main FISH-RT probe design pipeline"""

    console.print("[bold blue]🧬 smfish-like-rt-probe-designer[/bold blue]")
    console.print("=" * 60)

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Design FISH probes for RT-PCR with SNP analysis"
    )
    parser.add_argument("--config", help="Path to custom config file")
    parser.add_argument("--output", help="Output directory path", required=True)
    parser.add_argument("--genes", nargs="+", help="List of gene names to process")
    parser.add_argument("--test", action="store_true", help="Run with test gene set")
    args = parser.parse_args()

    # Load configuration
    config = FISH_RT_CONFIG.copy()
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

    # Check if required paths are configured
    missing_paths = []
    if not config["local_genome_directory"]:
        missing_paths.append("local_genome_directory")
    if not config["local_snp_file1"]:
        missing_paths.append("local_snp_file1")
    if not config["local_snp_file2"]:
        missing_paths.append("local_snp_file2")

    if missing_paths:
        console.print(f"[red]Warning: Missing paths in config: {missing_paths}[/red]")
        console.print("[yellow]Some features may be disabled[/yellow]")

    # Initialize components
    console.print("\n[bold]Initializing pipeline components...[/bold]")

    gene_fetcher = GeneSequenceFetcher(config)
    snp_analyzer = SNPCoverageAnalyzer(config) if config["local_snp_file1"] else None
    blast_analyzer = (
        BLASTSpecificityAnalyzer(config) if config["blast_database_path"] else None
    )
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
                console.print(f"✅ {gene_name}: {len(gene_data['sequence'])} bp")
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
        console.print("\n[bold cyan]Step 3: Analyzing SNP coverage...[/bold cyan]")
        try:
            all_probes = snp_analyzer.analyze_probes(all_probes)
            console.print("✅ SNP coverage analysis completed")
        except Exception as e:
            console.print(f"❌ SNP analysis failed: {str(e)}")
    else:
        console.print(
            "\n[yellow]Step 3: SNP analysis skipped (no SNP files configured)[/yellow]"
        )

    # Step 4: BLAST specificity analysis
    if blast_analyzer and config["run_blast_analysis"]:
        console.print("\n[bold cyan]Step 4: BLAST specificity analysis...[/bold cyan]")
        try:
            all_probes = blast_analyzer.analyze_probes(
                all_probes, output_path / "blast_results"
            )
            unique_probes = [p for p in all_probes if p.get("NumberOfHits") == 1]
            console.print(
                f"✅ BLAST completed: {len(unique_probes)}/{len(all_probes)} probes are unique"
            )
        except Exception as e:
            console.print(f"❌ BLAST analysis failed: {str(e)}")
    else:
        console.print(
            "\n[yellow]Step 4: BLAST analysis skipped (not configured)[/yellow]"
        )

    # Step 5: Generate final output
    console.print("\n[bold cyan]Step 5: Generating output files...[/bold cyan]")

    try:
        output_files = output_generator.generate_outputs(all_probes, output_path)

        console.print("[bold green]✅ Pipeline completed successfully![/bold green]")
        console.print(f"[green]Generated files:[/green]")
        for file_path in output_files:
            console.print(f"  📄 {file_path}")

        # Print summary statistics
        df_all = pd.read_csv(output_files[0])  # ALL file
        df_filt = pd.read_csv(output_files[1]) if len(output_files) > 1 else df_all

        console.print(f"\n[bold]📊 Summary Statistics:[/bold]")
        console.print(f"  Total probes designed: {len(df_all)}")
        console.print(f"  Probes after filtering: {len(df_filt)}")
        console.print(f"  Genes processed: {len(set(df_all['GeneName']))}")

        if "SNPs_Covered_Count" in df_all.columns:
            avg_snps = df_all["SNPs_Covered_Count"].mean()
            console.print(f"  Average SNPs covered per probe: {avg_snps:.1f}")

        if "NumberOfHits" in df_all.columns:
            unique_pct = (df_all["NumberOfHits"] == 1).mean() * 100
            console.print(f"  Probes with unique BLAST hits: {unique_pct:.1f}%")

    except Exception as e:
        console.print(f"[red]❌ Output generation failed: {str(e)}[/red]")
        return


if __name__ == "__main__":
    main()
