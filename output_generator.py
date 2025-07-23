# output_generator.py - Enhanced for Manual BLAST Workflow (No Local BLAST)

import pandas as pd
import os
from pathlib import Path
from rich.console import Console
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

console = Console()


class OutputGenerator:
    """Generate outputs optimized for manual NCBI BLAST workflow"""

    def __init__(self, config):
        self.config = config
        self.add_rtbc = config.get("add_rtbc_barcode", True)
        self.rtbc_sequence = config.get("rtbc_sequence", "/5Phos/TGACTTGAGGAT")
        self.generate_blast_fasta = config.get("generate_blast_fasta", True)
        self.fasta_description_format = config.get(
            "fasta_description_format", "detailed"
        )

    def generate_outputs(self, probe_list, output_dir):
        """
        Generate all output files optimized for manual BLAST workflow
        Returns list of generated file paths
        """
        output_dir = Path(output_dir)
        generated_files = []

        try:
            console.print(
                "[cyan]Generating outputs for manual BLAST workflow...[/cyan]"
            )

            # Add RTBC barcodes to all probes
            probe_list_with_rtbc = self._add_rtbc_barcodes(probe_list)

            # Create comprehensive DataFrame
            df_all = pd.DataFrame(probe_list_with_rtbc)

            # Sort by gene name, then by SNP coverage (descending), then by quality
            df_all = df_all.sort_values(
                ["GeneName", "SNPs_Covered_Count", "NbOfPNAS", "GCpc"],
                ascending=[True, False, False, False],
            )

            # Generate ALL probes file
            all_file = output_dir / "FISH_RT_probes_ALL.csv"
            df_all.to_csv(all_file, index=False)
            generated_files.append(all_file)
            console.print(f"✅ Generated: {all_file}")

            # Generate FILTERED probes file (GC + PNAS only, no BLAST)
            filtered_df = self._apply_quality_filters(df_all)
            filt_file = output_dir / "FISH_RT_probes_FILTERED.csv"
            filtered_df.to_csv(filt_file, index=False)
            generated_files.append(filt_file)
            console.print(f"✅ Generated: {filt_file}")

            # Generate summary statistics
            summary_file = self._generate_summary(df_all, filtered_df, output_dir)
            generated_files.append(summary_file)

            # Generate FASTA files for manual BLAST
            if self.generate_blast_fasta:
                fasta_files = self._generate_blast_fasta_files(
                    df_all, filtered_df, output_dir
                )
                generated_files.extend(fasta_files)

            console.print(
                f"[green]Generated {len(generated_files)} output files[/green]"
            )
            return generated_files

        except Exception as e:
            console.print(f"[red]Error generating outputs: {str(e)}[/red]")
            return []

    def _add_rtbc_barcodes(self, probe_list):
        """Add RTBC barcodes to probe sequences"""

        if not self.add_rtbc:
            # Just add empty column if RTBC is disabled
            for probe in probe_list:
                probe["RTBC_5Prime_Sequence"] = probe["Seq"]
            return probe_list

        console.print(f"[cyan]Adding RTBC barcode: {self.rtbc_sequence}[/cyan]")

        for probe in probe_list:
            # Add RTBC sequence to 5' end of probe
            rtbc_sequence = self.rtbc_sequence + probe["Seq"]
            probe["RTBC_5Prime_Sequence"] = rtbc_sequence

        return probe_list

    def _apply_quality_filters(self, df_all):
        """Apply quality filters (GC + PNAS only, no BLAST)"""

        # Start with all probes
        filtered = df_all.copy()

        console.print(f"Starting with {len(filtered)} probes")

        # Filter 1: Basic quality filters (GC + PNAS)
        basic_filter = (filtered["GCFilter"] == 1) & (filtered["PNASFilter"] == 1)
        filtered = filtered[basic_filter]
        console.print(f"After GC+PNAS filters: {len(filtered)} probes")

        # Filter 2: dustmasker (if enabled)
        if "MaskedFilter" in filtered.columns and self.config.get(
            "use_dustmasker", False
        ):
            masked_filter = filtered["MaskedFilter"] == 1
            filtered = filtered[masked_filter]
            console.print(f"After dustmasker filter: {len(filtered)} probes")

        return filtered

    def _generate_summary(self, df_all, df_filtered, output_dir):
        """Generate summary statistics file"""

        summary_file = output_dir / "FISH_RT_analysis_summary.txt"

        try:
            with open(summary_file, "w") as f:
                f.write("FISH-RT Probe Design Analysis Summary\n")
                f.write("=" * 50 + "\n\n")

                # Overall statistics
                f.write("OVERALL STATISTICS:\n")
                f.write(f"  Total probes designed: {len(df_all)}\n")
                f.write(f"  Probes after quality filtering: {len(df_filtered)}\n")
                f.write(
                    f"  Quality filtering efficiency: {len(df_filtered)/len(df_all)*100:.1f}%\n\n"
                )

                f.write("MANUAL BLAST WORKFLOW:\n")
                f.write(f"  1. Use the generated FASTA files for NCBI BLAST\n")
                f.write(f"  2. Submit to: https://blast.ncbi.nlm.nih.gov/Blast.cgi\n")
                f.write(f"  3. Use 'Mouse genomic + transcript' database\n")
                f.write(
                    f"  4. Filter results for unique hits (1 significant hit only)\n\n"
                )

                # Gene-level statistics
                f.write("GENE-LEVEL STATISTICS:\n")
                gene_stats = (
                    df_all.groupby("GeneName")
                    .agg(
                        {
                            "Seq": "count",
                            "SNPs_Covered_Count": "mean",
                            "GCFilter": "sum",
                            "PNASFilter": "sum",
                        }
                    )
                    .round(2)
                )
                gene_stats.columns = [
                    "Total_Probes",
                    "Avg_SNPs_Covered",
                    "GC_Pass",
                    "PNAS_Pass",
                ]

                f.write(gene_stats.to_string())
                f.write("\n\n")

                # Quality metrics
                if len(df_all) > 0:
                    f.write("QUALITY METRICS:\n")
                    f.write(f"  Average GC content: {df_all['GCpc'].mean():.3f}\n")
                    f.write(
                        f"  Average probe size: {df_all['ProbeSize'].mean():.1f} nt\n"
                    )
                    f.write(f"  Average dG37 score: {df_all['dGScore'].mean():.3f}\n")

                    if "SNPs_Covered_Count" in df_all.columns:
                        avg_snps = df_all["SNPs_Covered_Count"].mean()
                        f.write(f"  Average SNPs covered per probe: {avg_snps:.1f}\n")

                    f.write("\n")

                # Filter breakdown
                f.write("FILTER BREAKDOWN:\n")
                f.write(
                    f"  GC content filter pass: {(df_all['GCFilter'] == 1).sum()}/{len(df_all)}\n"
                )
                f.write(
                    f"  PNAS filter pass: {(df_all['PNASFilter'] == 1).sum()}/{len(df_all)}\n"
                )

                if "MaskedFilter" in df_all.columns:
                    masked_pass = (df_all["MaskedFilter"] == 1).sum()
                    f.write(f"  dustmasker filter pass: {masked_pass}/{len(df_all)}\n")

            console.print(f"✅ Generated summary: {summary_file}")
            return summary_file

        except Exception as e:
            console.print(f"[red]Error generating summary: {str(e)}[/red]")
            return None

    def _generate_blast_fasta_files(self, df_all, df_filtered, output_dir):
        """Generate FASTA files optimized for NCBI BLAST submission"""
        generated_files = []

        try:
            console.print("[cyan]Generating FASTA files for manual BLAST...[/cyan]")

            # FASTA file 1: ALL probes (probe sequences only, no RTBC)
            all_blast_fasta = output_dir / "ALL_probes_for_BLAST.fasta"
            self._write_blast_fasta(
                df_all, all_blast_fasta, "All probes for BLAST specificity check"
            )
            generated_files.append(all_blast_fasta)

            # FASTA file 2: FILTERED probes (probe sequences only, no RTBC)
            filt_blast_fasta = output_dir / "FILTERED_probes_for_BLAST.fasta"
            self._write_blast_fasta(
                df_filtered, filt_blast_fasta, "Quality-filtered probes for BLAST"
            )
            generated_files.append(filt_blast_fasta)

            # FASTA file 3: FILTERED probes with RTBC (for ordering)
            filt_rtbc_fasta = output_dir / "FILTERED_probes_with_RTBC.fasta"
            self._write_rtbc_fasta(df_filtered, filt_rtbc_fasta)
            generated_files.append(filt_rtbc_fasta)

            console.print(f"✅ Generated {len(generated_files)} FASTA files for BLAST")

            # Print BLAST instructions
            console.print("\n[bold yellow]📝 MANUAL BLAST INSTRUCTIONS:[/bold yellow]")
            console.print("1. Go to: https://blast.ncbi.nlm.nih.gov/Blast.cgi")
            console.print("2. Upload: FILTERED_probes_for_BLAST.fasta")
            console.print("3. Database: Mouse genomic + transcript")
            console.print("4. Look for probes with exactly 1 significant hit")
            console.print("5. Use FILTERED_probes_with_RTBC.fasta for synthesis")

            return generated_files

        except Exception as e:
            console.print(f"[red]Error generating BLAST FASTA files: {str(e)}[/red]")
            return []

    def _write_blast_fasta(self, df, output_file, description):
        """Write FASTA file optimized for BLAST (probe sequences only)"""

        records = []

        for idx, row in df.iterrows():
            # Use probe sequence only (no RTBC) for BLAST
            sequence = row["Seq"]
            gene_name = row["GeneName"]
            probe_id = f"probe_{idx}_{gene_name}"

            if self.fasta_description_format == "detailed":
                # Detailed description for tracking
                description = (
                    f"{gene_name} | {row['RegionType']} | "
                    f"chr{row['Chromosome']}:{row['theStartPos']}-{row['theEndPos']} | "
                    f"size:{row['ProbeSize']}nt | SNPs:{row.get('SNPs_Covered_Count', 0)} | "
                    f"GC:{row['GCpc']:.2f} | dG:{row['dG37']:.1f}"
                )
            else:
                # Simple description
                description = f"{gene_name}_probe_{idx}"

            record = SeqRecord(Seq(sequence), id=probe_id, description=description)
            records.append(record)

        SeqIO.write(records, output_file, "fasta")
        console.print(f"✅ BLAST FASTA: {output_file} ({len(records)} sequences)")

    def _write_rtbc_fasta(self, df, output_file):
        """Write FASTA file with RTBC sequences (for synthesis)"""

        records = []

        for idx, row in df.iterrows():
            # Use RTBC sequence for synthesis
            sequence = row["RTBC_5Prime_Sequence"]
            gene_name = row["GeneName"]
            probe_id = f"{gene_name}_RTBC_{idx}"

            description = (
                f"{gene_name} RT primer with RTBC | "
                f"Probe:{row['Seq']} | RTBC:{self.rtbc_sequence} | "
                f"Total length:{len(sequence)}nt"
            )

            record = SeqRecord(Seq(sequence), id=probe_id, description=description)
            records.append(record)

        SeqIO.write(records, output_file, "fasta")
        console.print(f"✅ RTBC FASTA: {output_file} ({len(records)} sequences)")
