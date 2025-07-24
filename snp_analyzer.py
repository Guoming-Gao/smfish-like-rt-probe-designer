# snp_analyzer.py - SNP Coverage Analysis (Local Files Only)

import pandas as pd
import numpy as np
from pathlib import Path
from rich.console import Console
from rich.progress import track

console = Console()


class SNPCoverageAnalyzer:
    """Analyze SNP coverage using local B6xCast SNP file"""

    def __init__(self, config):
        self.config = config
        self.rt_coverage_length = config.get("rt_coverage_downstream", 100)
        self.snp_file_path = config.get("snp_file_path", "")

        # Load SNP database
        self.snp_database = {}
        self.load_snp_database()

    def load_snp_database(self):
        """Load B6xCast SNP file"""
        if not self.snp_file_path or not Path(self.snp_file_path).exists():
            console.print(f"[red]SNP file not found: {self.snp_file_path}[/red]")
            return

        try:
            snps = self._parse_snp_file(self.snp_file_path)
            console.print(f"✅ Loaded {len(snps)} B6xCast SNPs")

            # Organize by chromosome
            for snp in snps:
                chr_name = snp["chromosome"]
                if chr_name not in self.snp_database:
                    self.snp_database[chr_name] = []
                self.snp_database[chr_name].append(snp)

            # Sort by position
            for chr_name in self.snp_database:
                self.snp_database[chr_name].sort(key=lambda x: x["start"])

            console.print(f"Organized SNPs across {len(self.snp_database)} chromosomes")

        except Exception as e:
            console.print(f"[red]Error loading SNP database: {str(e)}[/red]")
            self.snp_database = {}

    def _parse_snp_file(self, file_path):
        """Parse B6xCast SNP file"""
        snps = []

        try:
            with open(file_path, "r") as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue

                    try:
                        tokens = line.split("\t")
                        if len(tokens) >= 4:
                            snps.append(
                                {
                                    "chromosome": tokens[0],
                                    "start": int(tokens[1]),
                                    "end": int(tokens[2]),
                                    "genotype": tokens[3],
                                    "line_number": line_num,
                                }
                            )
                    except (ValueError, IndexError):
                        continue

        except Exception as e:
            console.print(f"[red]Error reading SNP file: {e}[/red]")

        return snps

    def analyze_probes(self, probe_list):
        """Analyze SNP coverage for all probes"""
        if not self.snp_database:
            console.print(
                "[yellow]No SNP database loaded, skipping SNP analysis[/yellow]"
            )
            return probe_list

        console.print(
            f"[cyan]Analyzing SNP coverage for {len(probe_list)} probes...[/cyan]"
        )

        for probe in track(probe_list):
            try:
                # Calculate RT coverage region in genomic coordinates
                rt_coverage = self._calculate_rt_coverage_coords(probe)

                # Find SNPs in RT region
                covered_snps = self._find_snps_in_region(
                    chromosome=rt_coverage["chromosome"],
                    start=rt_coverage["start"],
                    end=rt_coverage["end"],
                )

                # Update probe with SNP information
                probe.update(
                    {
                        "RT_Coverage_Start": rt_coverage["start"],
                        "RT_Coverage_End": rt_coverage["end"],
                        "RT_Coverage_Strand": rt_coverage["strand"],
                        "SNPs_Covered_Count": len(covered_snps),
                        "SNPs_Covered_Positions": ",".join(
                            [f"{s['chromosome']}:{s['start']}" for s in covered_snps]
                        ),
                        "SNPs_Covered_Types": ",".join(
                            [s["genotype"] for s in covered_snps]
                        ),
                    }
                )

            except Exception as e:
                console.print(f"[red]Error analyzing SNP coverage: {e}[/red]")
                probe.update(
                    {
                        "RT_Coverage_Start": 0,
                        "RT_Coverage_End": 0,
                        "RT_Coverage_Strand": ".",
                        "SNPs_Covered_Count": 0,
                        "SNPs_Covered_Positions": "",
                        "SNPs_Covered_Types": "",
                    }
                )

        return probe_list

    def _calculate_rt_coverage_coords(self, probe):
        """Calculate RT coverage region (strand-aware)"""
        target_strand = probe["Target_Strand"]
        probe_start = probe["theStartPos"]
        probe_end = probe["theEndPos"]
        chromosome = probe["Chromosome"]

        if target_strand == "+":
            rt_start = probe_end + 1
            rt_end = rt_start + self.rt_coverage_length - 1
            rt_strand = "+"
        else:
            rt_end = probe_start - 1
            rt_start = rt_end - self.rt_coverage_length + 1
            rt_strand = "-"

        return {
            "chromosome": chromosome,
            "start": rt_start,
            "end": rt_end,
            "strand": rt_strand,
        }

    def _find_snps_in_region(self, chromosome, start, end):
        """Find SNPs overlapping with region"""
        # Handle chromosome naming (with or without 'chr' prefix)
        chr_variants = [chromosome, f"chr{chromosome}", chromosome.replace("chr", "")]

        overlapping_snps = []
        for chr_name in chr_variants:
            if chr_name in self.snp_database:
                chr_snps = self.snp_database[chr_name]
                for snp in chr_snps:
                    if snp["start"] <= end and snp["end"] >= start:
                        overlapping_snps.append(snp)
                break  # Found the right chromosome format

        return overlapping_snps
