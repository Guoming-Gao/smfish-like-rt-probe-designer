#!/usr/bin/env python3

# main.py - smfish-like-rt-probe-designer

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
from utils.sequence_utils import validate_probe_strand_complementarity
from snp_analyzer import SNPCoverageAnalyzer
from output_generator import OutputGenerator

console = Console()


def setup_output_directory(output_path):
    """Create output directory structure"""
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    (output_path / "gene_sequences").mkdir(exist_ok=True)
    return output_path


def apply_early_filtering(probes, config):
    """Apply stringent filtering early to focus on high-quality probes only"""
    console.print(
        f"[cyan]Applying stringent filtering to {len(probes)} probes...[/cyan]"
    )

    filtered_probes = []
    filter_stats = {
        "total": len(probes),
        "gc_pass": 0,
        "pnas_pass": 0,
        "dustmasker_pass": 0,
        "all_filters_pass": 0,
    }

    for probe in probes:
        # Check individual filters
        gc_pass = probe["GCFilter"] == 1
        pnas_pass = probe["PNASFilter"] == 1
        dustmasker_pass = probe.get("MaskedFilter", 1) == 1  # Default pass if not set

        # Update statistics
        if gc_pass:
            filter_stats["gc_pass"] += 1
        if pnas_pass:
            filter_stats["pnas_pass"] += 1
        if dustmasker_pass:
            filter_stats["dustmasker_pass"] += 1

        # Apply ALL filters (stringent)
        if gc_pass and pnas_pass and dustmasker_pass:
            filtered_probes.append(probe)
            filter_stats["all_filters_pass"] += 1

    # Report filtering statistics
    console.print(f"[green]Filtering Results:[/green]")
    console.print(f" Total probes: {filter_stats['total']}")
    console.print(f" GC filter pass: {filter_stats['gc_pass']}")
    console.print(f" PNAS filter pass: {filter_stats['pnas_pass']}")
    console.print(f" Dustmasker pass: {filter_stats['dustmasker_pass']}")
    console.print(f" All filters pass: {filter_stats['all_filters_pass']}")
    console.print(
        f" Retention rate: {filter_stats['all_filters_pass']/filter_stats['total']*100:.1f}%"
    )

    return filtered_probes


def validate_probe_orientations(probes, gene_data_list, config):
    """
    NEW: Validate probe strand orientations against gene strands
    """
    console.print(f"[cyan]🔍 Validating probe strand orientations...[/cyan]")

    # Create gene lookup
    gene_lookup = {gene["gene_name"]: gene for gene in gene_data_list}

    validation_results = []
    issues_found = 0

    for probe in probes:
        gene_name = probe.get("GeneName")
        if gene_name in gene_lookup:
            gene_data = gene_lookup[gene_name]
            validation = validate_probe_strand_complementarity(probe, gene_data)
            validation_results.append(validation)

            if not validation["overall_validation"]:
                issues_found += 1
                console.print(
                    f"[red]❌ {gene_name}: Gene strand {validation['gene_strand']}, "
                    f"Target {validation['target_strand']}, Probe {validation['probe_strand']} "
                    f"(Expected: Target {validation['expected_target']}, Probe {validation['expected_probe']})[/red]"
                )

    # Summary
    total_validated = len(validation_results)
    if issues_found == 0:
        console.print(
            f"[green]✅ All {total_validated} probes have correct strand orientations![/green]"
        )
    else:
        console.print(
            f"[red]❌ Found {issues_found}/{total_validated} probes with incorrect orientations![/red]"
        )
        console.print(
            f"[red]   This indicates a serious bug in sequence preparation.[/red]"
        )

    return validation_results


def main():
    """Main FISH-RT probe design pipeline (FIXED: Added validation)"""
    console.print(
        "[bold blue]🧬 smfish-like-rt-probe-designer (FIXED: Strand validation)[/bold blue]"
    )
    console.print("=" * 60)
    console.print("[cyan]Mode: Local GTF + FASTA files only + Strand validation[/cyan]")

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Design high-quality FISH probes using local files"
    )
    parser.add_argument(
        "--output", help="Output directory path (optional, uses config default)"
    )
    parser.add_argument("--genes", nargs="+", help="List of gene names to process")
    parser.add_argument("--test", action="store_true", help="Run with test gene set")
    parser.add_argument(
        "--validate-only", action="store_true", help="Only run strand validation"
    )
    args = parser.parse_args()

    # Load configuration
    config = FISH_RT_CONFIG.copy()

    # Use config default if --output not provided
    if args.output:
        config["output_directory"] = args.output
        console.print(f"[yellow]Using custom output directory: {args.output}[/yellow]")
    else:
        console.print(f"[green]Using config default output directory[/green]")

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

    # Show configuration
    console.print(f"\n[bold]Configuration:[/bold]")
    console.print(f" PNAS rules: {config['pnas_filter_rules']}")
    console.print(
        f" Dustmasker: {'✅ Enabled' if config['use_dustmasker'] else '❌ Disabled'}"
    )
    console.print(f" Min SNP coverage: {config['min_snp_coverage_for_final']}")
    console.print(f" RT coverage: {config['rt_coverage_downstream']} nt downstream")

    # Initialize components
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
                    f"✅ {gene_name}: {len(gene_data['sequence'])} bp ({gene_data['source']}) strand {gene_data['strand']}"
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

    # Save gene sequences (for validation)
    gene_fetcher.save_gene_sequences(all_gene_data, output_path / "gene_sequences")

    # NEW: Validate gene sequence orientations
    validation_results = gene_fetcher.validate_strand_orientation(all_gene_data)

    if args.validate_only:
        console.print("[yellow]Validation-only mode complete. Exiting.[/yellow]")
        return

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
                console.print(f"⚠️ {gene_data['gene_name']}: No probes generated")
        except Exception as e:
            console.print(f"❌ {gene_data['gene_name']}: Error - {str(e)}")

    if not all_probes:
        console.print("[red]No probes generated. Exiting.[/red]")
        return

    console.print(f"[green]Total probes designed: {len(all_probes)}[/green]")

    # NEW: Step 2.5: Validate probe strand orientations
    console.print("\n[bold cyan]Step 2.5: Validating probe orientations...[/bold cyan]")
    probe_validations = validate_probe_orientations(all_probes, all_gene_data, config)

    # Step 3: Early stringent filtering
    console.print("\n[bold cyan]Step 3: Early stringent filtering...[/bold cyan]")
    high_quality_probes = apply_early_filtering(all_probes, config)

    if not high_quality_probes:
        console.print(
            "[red]No probes passed stringent filtering. Consider relaxing filter parameters.[/red]"
        )
        return

    console.print(f"[green]High-quality probes: {len(high_quality_probes)}[/green]")

    # Step 4: SNP analysis (only on filtered probes) - UPDATED WITH PROGRESS BAR
    if snp_analyzer:
        console.print(
            f"\n[bold cyan]Step 4: SNP analysis (filtered probes only)...[/bold cyan]"
        )
        try:
            # UPDATED: Add progress bar for SNP analysis
            updated_probes = []
            for probe in track(
                high_quality_probes, description="Analyzing SNP coverage..."
            ):
                # Try individual probe analysis first
                try:
                    if hasattr(snp_analyzer, "analyze_probe"):
                        updated_probe = snp_analyzer.analyze_probe(probe)
                    else:
                        # Fallback: analyze single probe using analyze_probes
                        updated_probe = snp_analyzer.analyze_probes([probe])[0]
                    updated_probes.append(updated_probe)
                except Exception as probe_error:
                    console.print(
                        f"⚠️ SNP analysis failed for probe: {str(probe_error)}"
                    )
                    updated_probes.append(
                        probe
                    )  # Keep original probe if analysis fails

            high_quality_probes = updated_probes

            # Report SNP coverage statistics
            snp_counts = [p.get("SNPs_Covered_Count", 0) for p in high_quality_probes]
            avg_snps = sum(snp_counts) / len(snp_counts) if snp_counts else 0
            max_snps = max(snp_counts) if snp_counts else 0
            high_snp_count = sum(
                1
                for count in snp_counts
                if count >= config["min_snp_coverage_for_final"]
            )

            console.print(f"✅ SNP coverage analysis completed")
            console.print(f" Average SNPs per probe: {avg_snps:.1f}")
            console.print(f" Maximum SNPs per probe: {max_snps}")
            console.print(
                f" Probes with ≥{config['min_snp_coverage_for_final']} SNPs: {high_snp_count}"
            )

        except Exception as e:
            console.print(f"❌ SNP analysis failed: {str(e)}")
    else:
        console.print(
            "\n[yellow]Step 4: SNP analysis skipped (no SNP file configured)[/yellow]"
        )

    # Step 5: Generate focused output files
    console.print(
        f"\n[bold cyan]Step 5: Generating focused output files...[/bold cyan]"
    )

    try:
        output_files = output_generator.generate_focused_outputs(
            high_quality_probes, output_path
        )

        console.print("[bold green]✅ Pipeline completed successfully![/bold green]")
        console.print(f"[green]Generated files:[/green]")
        for file_path in output_files:
            console.print(f" 📄 {file_path}")

        # Final summary
        if output_files:
            csv_files = [f for f in output_files if str(f).endswith(".csv")]
            if csv_files:
                df_filt = pd.read_csv(csv_files[0])
                console.print(f"\n[bold]📊 Final Summary:[/bold]")
                console.print(f" High-quality probes: {len(df_filt)}")
                console.print(f" Genes processed: {len(set(df_filt['GeneName']))}")
                console.print(f" Output directory: {output_path}")

                if "SNPs_Covered_Count" in df_filt.columns:
                    avg_snps = df_filt["SNPs_Covered_Count"].mean()
                    high_snp_final = (
                        df_filt["SNPs_Covered_Count"]
                        >= config["min_snp_coverage_for_final"]
                    ).sum()
                    console.print(f" Average SNPs covered: {avg_snps:.1f}")
                    console.print(
                        f" Probes with ≥{config['min_snp_coverage_for_final']} SNPs: {high_snp_final}"
                    )

    except Exception as e:
        console.print(f"[red]❌ Output generation failed: {str(e)}[/red]")


if __name__ == "__main__":
    main()
