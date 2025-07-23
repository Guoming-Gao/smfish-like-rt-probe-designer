# snp_analyzer.py - SNP Coverage Analysis with Coordinate Conversion

import pandas as pd
import numpy as np
from pathlib import Path
from rich.console import Console

console = Console()


class SNPCoverageAnalyzer:
    """Analyze SNP coverage in RT region with coordinate conversion (gene-focused)"""

    def __init__(self, config):
        self.config = config
        self.rt_coverage_length = config.get("rt_coverage_downstream", 100)
        self.snp_file_path = config.get("snp_file_path", "")

        # Load SNP database (single file)
        self.snp_database = {}
        self.load_snp_database()

    def load_snp_database(self):
        """Load SNP database from single B6xCast file"""
        if not self.snp_file_path or not Path(self.snp_file_path).exists():
            console.print(f"[red]SNP file not found: {self.snp_file_path}[/red]")
            return

        try:
            snps = self._parse_snp_file(self.snp_file_path)
            console.print(f"✅ Loaded {len(snps)} B6xCast SNPs from mm10 coordinates")

            # Organize by chromosome for efficient lookup
            for snp in snps:
                chr_name = snp["chromosome"]
                if chr_name not in self.snp_database:
                    self.snp_database[chr_name] = []
                self.snp_database[chr_name].append(snp)

            # Sort SNPs by position for efficient searching
            for chr_name in self.snp_database:
                self.snp_database[chr_name].sort(key=lambda x: x["start"])

            console.print(f"Organized SNPs across {len(self.snp_database)} chromosomes")

        except Exception as e:
            console.print(f"[red]Error loading SNP database: {str(e)}[/red]")
            self.snp_database = {}

    def _parse_snp_file(self, file_path):
        """Parse SNP file with format: chr1 \t 3001490 \t 3001490 \t C/C,A/A"""
        snps = []

        try:
            with open(file_path, "r") as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue

                    try:
                        # Format: chromosome \t start \t end \t genotype_info
                        tokens = line.split("\t")
                        if len(tokens) >= 4:
                            snps.append(
                                {
                                    "chromosome": tokens[0],  # chr1
                                    "start": int(tokens[1]),  # 3001490
                                    "end": int(tokens[2]),  # 3001490
                                    "genotype": tokens[3],  # C/C,A/A (ignore for now)
                                    "line_number": line_num,
                                }
                            )
                    except (ValueError, IndexError) as e:
                        console.print(
                            f"[yellow]Skipping malformed line {line_num}: {line}[/yellow]"
                        )
                        continue

        except Exception as e:
            console.print(f"[red]Error reading SNP file {file_path}: {str(e)}[/red]")

        return snps

    def analyze_probes(self, probe_list):
        """
        Analyze SNP coverage for all probes
        KEY: Only analyze SNPs within genes of interest
        """
        if not self.snp_database:
            console.print(
                "[yellow]No SNP database loaded, skipping SNP analysis[/yellow]"
            )
            return probe_list

        console.print(
            f"[cyan]Analyzing SNP coverage for {len(probe_list)} probes...[/cyan]"
        )

        for probe in probe_list:
            try:
                # Calculate RT coverage region in GENOMIC coordinates
                rt_coverage_genomic = self._calculate_rt_coverage_genomic_coords(probe)

                # Find SNPs in RT coverage region (genomic coordinates)
                covered_snps = self._find_snps_in_genomic_region(
                    chromosome=rt_coverage_genomic["chromosome"],
                    start=rt_coverage_genomic["start"],
                    end=rt_coverage_genomic["end"],
                )

                # Update probe with SNP information
                probe.update(
                    {
                        "RT_Coverage_Start": rt_coverage_genomic["start"],
                        "RT_Coverage_End": rt_coverage_genomic["end"],
                        "RT_Coverage_Strand": rt_coverage_genomic["strand"],
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
                console.print(
                    f"[red]Error analyzing SNP coverage for probe {probe.get('Seq', 'unknown')}: {str(e)}[/red]"
                )
                # Set default values on error
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

    def _calculate_rt_coverage_genomic_coords(self, probe):
        """
        Calculate RT coverage region in GENOMIC coordinates (mm10)
        CRITICAL: Strand-aware and converts relative → genomic coordinates
        """
        # Probe positions are already in genomic coordinates from oligostan_core
        target_strand = probe["Target_Strand"]  # Gene's RNA strand
        probe_start_genomic = probe["theStartPos"]  # Already genomic coords
        probe_end_genomic = probe["theEndPos"]  # Already genomic coords
        chromosome = probe["Chromosome"]

        if target_strand == "+":
            # RNA goes 5' → 3' in + direction
            # RT coverage is downstream (towards higher coordinates)
            rt_start = probe_end_genomic + 1  # Start after probe ends
            rt_end = rt_start + self.rt_coverage_length - 1
            rt_strand = "+"

        else:  # target_strand == '-'
            # RNA goes 5' → 3' in - direction
            # RT coverage is downstream (towards lower coordinates)
            rt_end = probe_start_genomic - 1  # End before probe starts
            rt_start = rt_end - self.rt_coverage_length + 1
            rt_strand = "-"

        return {
            "chromosome": chromosome,
            "start": rt_start,
            "end": rt_end,
            "strand": rt_strand,
            "length": self.rt_coverage_length,
        }

    def _find_snps_in_genomic_region(self, chromosome, start, end):
        """Find all SNPs that overlap with the genomic region"""
        if chromosome not in self.snp_database:
            return []

        chr_snps = self.snp_database[chromosome]
        overlapping_snps = []

        for snp in chr_snps:
            # Check if SNP overlaps with RT coverage region
            if snp["start"] <= end and snp["end"] >= start:
                overlapping_snps.append(snp)

        return overlapping_snps
