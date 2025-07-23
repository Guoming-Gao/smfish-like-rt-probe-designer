# output_generator.py - Generate Consolidated Outputs for FISH-RT Analysis

import pandas as pd
import os
from pathlib import Path
from rich.console import Console
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

console = Console()


class OutputGenerator:
    """Generate consolidated output files for FISH-RT probe analysis"""

    def __init__(self, config):
        self.config = config
        self.add_rtbc = config.get("add_rtbc_barcode", True)
        self.rtbc_sequence = config.get("rtbc_sequence", "/5Phos/TGACTTGAGGAT")
        self.require_unique_hits = config.get("require_unique_hits_only", True)

    def generate_outputs(self, probe_list, output_dir):
        """
        Generate all output files for FISH-RT analysis
        Returns list of generated file paths
        """
        output_dir = Path(output_dir)
        generated_files = []

        try:
            console.print("[cyan]Generating output files...[/cyan]")

            # Add RTBC barcodes to all probes
            probe_list_with_rtbc = self._add_rtbc_barcodes(probe_list)

            # Create comprehensive DataFrame
            df_all = pd.DataFrame(probe_list_with_rtbc)

            # Sort by gene name, then by SNP coverage (descending)
            df_all = df_all.sort_values(
                ["GeneName", "SNPs_Covered_Count", "NbOfPNAS"],
                ascending=[True, False, False],
            )

            # Generate ALL probes file
            all_file = output_dir / "FISH_RT_probes_ALL.csv"
            df_all.to_csv(all_file, index=False)
            generated_files.append(all_file)
            console.print(f"✅ Generated: {all_file}")

            # Generate FILTERED probes file
            filtered_df = self._apply_filters(df_all)
            filt_file = output_dir / "FISH_RT_probes_FILT.csv"
            filtered_df.to_csv(filt_file, index=False)
            generated_files.append(filt_file)
            console.print(f"✅ Generated: {filt_file}")

            # Generate summary statistics
            summary_file = self._generate_summary(df_all, filtered_df, output_dir)
            generated_files.append(summary_file)

            # Generate FASTA files
            fasta_files = self._generate_fasta_files(df_all, filtered_df, output_dir)
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

    def _apply_filters(self, df_all):
        """Apply filtering criteria to get final probe set"""

        # Start with all probes
        filtered = df_all.copy()

        console.print(f"Starting with {len(filtered)} probes")

        # Filter 1: Basic quality filters (GC + PNAS)
        basic_filter = (filtered["GCFilter"] == 1) & (filtered["PNASFilter"] == 1)
        filtered = filtered[basic_filter]
        console.print(f"After GC+PNAS filters: {len(filtered)} probes")

        # Filter 2: BLAST uniqueness (if required)
        if self.require_unique_hits and "NumberOfHits" in filtered.columns:
            blast_filter = filtered["NumberOfHits"] == 1
            filtered = filtered[blast_filter]
            console.print(f"After BLAST uniqueness filter: {len(filtered)} probes")

        # Filter 3: dustmasker (if enabled)
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
                f.write(f"  Probes after filtering: {len(df_filtered)}\n")
                f.write(
                    f"  Filtering efficiency: {len(df_filtered)/len(df_all)*100:.1f}%\n\n"
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
                        f.write(
                            f"  Average SNPs covered: {df_all['SNPs_Covered_Count'].mean():.1f}\n"
                        )

                    if "NumberOfHits" in df_all.columns:
                        unique_pct = (df_all["NumberOfHits"] == 1).mean() * 100
                        f.write(f"  Probes with unique BLAST hits: {unique_pct:.1f}%\n")

                    f.write("\n")

                # Filter breakdown
                f.write("FILTER BREAKDOWN:\n")
                f.write(
                    f"  GC content filter pass: {(df_all['GCFilter'] == 1).sum()}/{len(df_all)}\n"
                )
                f.write(
                    f"  PNAS filter pass: {(df_all['PNASFilter'] == 1).sum()}/{len(df_all)}\n"
                )

                if "NumberOfHits" in df_all.columns:
                    unique_count = (df_all["NumberOfHits"] == 1).sum()
                    f.write(f"  BLAST unique hits: {unique_count}/{len(df_all)}\n")

                if "MaskedFilter" in df_all.columns:
                    masked_pass = (df_all["MaskedFilter"] == 1).sum()
                    f.write(f"  dustmasker filter pass: {masked_pass}/{len(df_all)}\n")

            console.print(f"✅ Generated summary: {summary_file}")
            return summary_file

        except Exception as e:
            console.print(f"[red]Error generating summary: {str(e)}[/red]")
            return None

    def _generate_fasta_files(self, df_all, df_filtered, output_dir):
        """Generate FASTA files for probe sequences"""
        generated_files = []

        try:
            # FASTA file for ALL probes (with RTBC)
            all_fasta = output_dir / "FISH_RT_probes_ALL.fasta"
            self._write_fasta(df_all, all_fasta, include_rtbc=True)
            generated_files.append(all_fasta)

            # FASTA file for FILTERED probes (with RTBC)
            filt_fasta = output_dir / "FISH_RT_probes_FILT.fasta"
            self._write_fasta(df_filtered, filt_fasta, include_rtbc=True)
            generated_files.append(filt_fasta)

            # FASTA file for probe sequences only (no RTBC)
            probes_only_fasta = output_dir / "FISH_RT_probe_sequences_only.fasta"
            self._write_fasta(df_filtered, probes_only_fasta, include_rtbc=False)
            generated_files.append(probes_only_fasta)

            console.print(f"✅ Generated {len(generated_files)} FASTA files")
            return generated_files

        except Exception as e:
            console.print(f"[red]Error generating FASTA files: {str(e)}[/red]")
            return []

    def _write_fasta(self, df, output_file, include_rtbc=True):
        """Write DataFrame to FASTA file"""

        records = []

        for idx, row in df.iterrows():
            # Choose sequence based on RTBC preference
            if include_rtbc and self.add_rtbc:
                sequence = row["RTBC_5Prime_Sequence"]
                seq_type = "with_RTBC"
            else:
                sequence = row["Seq"]
                seq_type = "probe_only"

            # Create descriptive FASTA header
            gene_name = row["GeneName"]
            probe_size = row["ProbeSize"]
            snp_count = row.get("SNPs_Covered_Count", 0)
            region_type = row.get("RegionType", "unknown")

            description = f"{gene_name} | {seq_type} | size:{probe_size}nt | SNPs:{snp_count} | region:{region_type}"

            record = SeqRecord(
                Seq(sequence), id=f"{gene_name}_probe_{idx}", description=description
            )
            records.append(record)

        SeqIO.write(records, output_file, "fasta")
