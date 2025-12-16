# snp_analyzer.py - SNP Coverage Analysis (VCF format with pysam/tabix)

import pandas as pd
import numpy as np
from pathlib import Path
from rich.console import Console
from rich.progress import track

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False

console = Console()


class SNPCoverageAnalyzer:
    """Analyze SNP coverage using VCF file with tabix index (B6 vs Cast differences)"""

    def __init__(self, config):
        self.config = config
        self.rt_coverage_length = config.get("rt_coverage_downstream", 100)
        self.snp_file_path = config.get("snp_file_path", "")
        self.b6_sample = config.get("vcf_b6_sample", "C57BL_6NJ")
        self.cast_sample = config.get("vcf_cast_sample", "CAST_EiJ")

        # VCF file handle and sample indices
        self.vcf_file = None
        self.b6_idx = None
        self.cast_idx = None
        self.samples = []

        # Initialize VCF connection
        self._init_vcf()

    def _init_vcf(self):
        """Initialize VCF file with tabix"""
        if not PYSAM_AVAILABLE:
            console.print("[red]pysam not installed. Run: pip install pysam[/red]")
            return

        if not self.snp_file_path or not Path(self.snp_file_path).exists():
            console.print(f"[red]VCF file not found: {self.snp_file_path}[/red]")
            return

        # Check for tabix index
        tbi_path = f"{self.snp_file_path}.tbi"
        if not Path(tbi_path).exists():
            console.print(f"[red]Tabix index not found: {tbi_path}[/red]")
            console.print("[yellow]Run: tabix -p vcf {self.snp_file_path}[/yellow]")
            return

        try:
            self.vcf_file = pysam.TabixFile(self.snp_file_path)

            # Parse header to get sample indices
            self._parse_vcf_header()

            if self.b6_idx is not None and self.cast_idx is not None:
                console.print(f"[green]✅ VCF loaded: {self.snp_file_path}[/green]")
                console.print(f"[cyan]   B6 sample: {self.b6_sample} (column {self.b6_idx})[/cyan]")
                console.print(f"[cyan]   Cast sample: {self.cast_sample} (column {self.cast_idx})[/cyan]")
            else:
                console.print(f"[red]Could not find B6/Cast samples in VCF[/red]")

        except Exception as e:
            console.print(f"[red]Error loading VCF: {str(e)}[/red]")
            import traceback
            traceback.print_exc()
            self.vcf_file = None

    def _parse_vcf_header(self):
        """Parse VCF header to find sample column indices"""
        try:
            # Get header from VCF
            header_line = None
            for line in self.vcf_file.header:
                line = line.decode() if isinstance(line, bytes) else line
                if line.startswith("#CHROM"):
                    header_line = line
                    break

            if not header_line:
                console.print("[red]No #CHROM header line found in VCF[/red]")
                return

            # Parse header columns
            columns = header_line.strip().split("\t")
            # Standard VCF columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, then samples
            # Sample columns start at index 9
            self.samples = columns[9:]

            # Find B6 and Cast indices
            if self.b6_sample in self.samples:
                self.b6_idx = self.samples.index(self.b6_sample)
            else:
                console.print(f"[red]B6 sample '{self.b6_sample}' not found in VCF[/red]")
                console.print(f"[yellow]Available samples: {', '.join(self.samples[:10])}...[/yellow]")

            if self.cast_sample in self.samples:
                self.cast_idx = self.samples.index(self.cast_sample)
            else:
                console.print(f"[red]Cast sample '{self.cast_sample}' not found in VCF[/red]")

        except Exception as e:
            console.print(f"[red]Error parsing VCF header: {e}[/red]")

    def load_snp_database(self):
        """Legacy method - VCF uses on-demand queries, no preloading needed"""
        pass  # VCF queries are done on-demand via tabix

    def analyze_probes(self, probe_list):
        """Analyze SNP coverage for all probes"""
        if not self.vcf_file or self.b6_idx is None or self.cast_idx is None:
            console.print(
                "[yellow]VCF not properly loaded, skipping SNP analysis[/yellow]"
            )
            # Return probes with empty SNP fields
            for probe in probe_list:
                probe.update({
                    "RT_Coverage_Start": 0,
                    "RT_Coverage_End": 0,
                    "RT_Coverage_Strand": ".",
                    "SNPs_Covered_Count": 0,
                    "SNPs_Covered_Positions": "",
                    "SNPs_Covered_Types": "",
                })
            return probe_list

        console.print(
            f"[cyan]Analyzing SNP coverage for {len(probe_list)} probes...[/cyan]"
        )

        for probe in track(probe_list):
            try:
                # Calculate RT coverage region in genomic coordinates
                rt_coverage = self._calculate_rt_coverage_coords(probe)

                # Find SNPs in RT region (B6 vs Cast differences)
                covered_snps = self._find_snps_in_region_vcf(
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
                            [f"{s['chromosome']}:{s['pos']}" for s in covered_snps]
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

    def _find_snps_in_region_vcf(self, chromosome, start, end):
        """Find SNPs in region where B6 and Cast have different genotypes"""
        overlapping_snps = []

        # VCF uses chromosome names without 'chr' prefix
        vcf_chrom = chromosome.replace("chr", "") if chromosome.startswith("chr") else chromosome

        try:
            # Query VCF using tabix
            # pysam.TabixFile.fetch returns iterator of lines
            for line in self.vcf_file.fetch(vcf_chrom, start - 1, end):  # tabix uses 0-based
                line = line.decode() if isinstance(line, bytes) else line
                fields = line.strip().split("\t")

                if len(fields) < 10:
                    continue

                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]

                # Get sample genotypes (fields[9:] are samples)
                sample_fields = fields[9:]
                if len(sample_fields) <= max(self.b6_idx, self.cast_idx):
                    continue

                b6_data = sample_fields[self.b6_idx]
                cast_data = sample_fields[self.cast_idx]

                # Parse genotypes
                b6_gt = self._parse_genotype(b6_data, ref, alt)
                cast_gt = self._parse_genotype(cast_data, ref, alt)

                # Skip if either genotype is missing
                if b6_gt is None or cast_gt is None:
                    continue

                # Only include SNPs where B6 and Cast are DIFFERENT
                if b6_gt != cast_gt:
                    overlapping_snps.append({
                        "chromosome": f"chr{chrom}",  # Add chr prefix for consistency
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "b6_gt": b6_gt,
                        "cast_gt": cast_gt,
                        "genotype": f"{b6_gt},{cast_gt}",  # B6/Cast alleles
                    })

        except ValueError as e:
            # Chromosome not in VCF (e.g., chrM, chrY)
            pass
        except Exception as e:
            console.print(f"[yellow]Warning: Error querying {vcf_chrom}:{start}-{end}: {e}[/yellow]")

        return overlapping_snps

    def _parse_genotype(self, sample_data, ref, alt):
        """Parse VCF genotype field to actual alleles"""
        # Format: GT:GQ:DP:... (GT is first field)
        gt_field = sample_data.split(":")[0]

        # Missing data
        if "." in gt_field:
            return None

        # Parse allele indices
        try:
            alleles = [ref] + alt.split(",")  # Handle multi-allelic
            indices = [int(i) for i in gt_field.replace("|", "/").split("/")]
            return "/".join(alleles[i] for i in indices)
        except (ValueError, IndexError):
            return None

    def __del__(self):
        """Clean up VCF file handle"""
        if self.vcf_file:
            self.vcf_file.close()
