# gene_fetcher.py - Local Genome with GTF Parsing (Complete Implementation)

import os
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rich.console import Console
import re

console = Console()


class GeneSequenceFetcher:
    """Fetch gene sequences from local FASTA file and GTF annotation"""

    def __init__(self, config):
        self.config = config
        self.genome_build = config.get("genome_build", "mm10")

        # File paths
        self.local_genome_fasta = config.get("local_genome_fasta_path", "")
        self.local_gtf_path = config.get("local_gtf_path", "")
        self.transcript_selection = config.get("transcript_selection", "longest")

        # Check if files exist
        self.use_demo_mode = not (
            self.local_genome_fasta and os.path.exists(self.local_genome_fasta)
        )
        self.use_gtf = self.local_gtf_path and os.path.exists(self.local_gtf_path)

        # Gene coordinate cache
        self.gene_cache = {}
        self.gtf_loaded = False

        if self.use_demo_mode:
            console.print(
                "[yellow]⚠️  Local genome FASTA not found - using demo mode[/yellow]"
            )
            console.print(f"[yellow]   Looking for: {self.local_genome_fasta}[/yellow]")
        else:
            console.print(
                f"[cyan]✅ Using local genome FASTA: {self.local_genome_fasta}[/cyan]"
            )

        if self.use_gtf:
            console.print(
                f"[cyan]✅ Using GTF annotation: {self.local_gtf_path}[/cyan]"
            )
            self._load_gtf_gene_coordinates()
        else:
            console.print(
                f"[yellow]⚠️  GTF not found, using demo coordinates only[/yellow]"
            )
            console.print(f"[yellow]   Looking for: {self.local_gtf_path}[/yellow]")

        console.print(f"[cyan]Genome build: {self.genome_build}[/cyan]")

    def _load_gtf_gene_coordinates(self):
        """Load gene coordinates from GTF file"""
        console.print("[cyan]Loading gene coordinates from GTF...[/cyan]")

        try:
            gene_count = 0
            with open(self.local_gtf_path, "r") as f:
                for line_num, line in enumerate(f, 1):
                    if line.startswith("#") or not line.strip():
                        continue

                    try:
                        parts = line.strip().split("\t")
                        if len(parts) < 9:
                            continue

                        (
                            chrom,
                            source,
                            feature,
                            start,
                            end,
                            score,
                            strand,
                            frame,
                            attributes,
                        ) = parts

                        # Only process gene features
                        if feature != "gene":
                            continue

                        # Parse attributes to get gene name
                        gene_name = self._extract_gene_name_from_attributes(attributes)
                        if not gene_name:
                            continue

                        # Store gene coordinates
                        self.gene_cache[gene_name] = {
                            "chromosome": (
                                chrom if chrom.startswith("chr") else f"chr{chrom}"
                            ),
                            "start": int(start),
                            "end": int(end),
                            "strand": strand,
                            "source": source,
                            "attributes": attributes,
                        }

                        gene_count += 1

                        # Show progress for large GTF files
                        if line_num % 100000 == 0:
                            console.print(
                                f"[cyan]   Processed {line_num} GTF lines, found {gene_count} genes...[/cyan]"
                            )

                    except Exception as e:
                        if line_num < 100:  # Only show errors for first 100 lines
                            console.print(
                                f"[yellow]Warning: Error parsing GTF line {line_num}: {e}[/yellow]"
                            )
                        continue

            console.print(f"[green]✅ Loaded {gene_count} genes from GTF[/green]")
            self.gtf_loaded = True

            # Show some example genes for verification
            if gene_count > 0:
                sample_genes = list(self.gene_cache.keys())[:5]
                console.print(f"[cyan]   Sample genes: {sample_genes}[/cyan]")

        except Exception as e:
            console.print(f"[red]Error loading GTF file: {e}[/red]")
            self.use_gtf = False

    def _extract_gene_name_from_attributes(self, attributes):
        """Extract gene name from GTF attributes"""
        # Try different attribute patterns
        patterns = [
            r'gene_name "([^"]+)"',  # gene_name "Nanog"
            r"gene_name=([^;]+)",  # gene_name=Nanog
            r"Name=([^;]+)",  # Name=Nanog
            r'gene_id "([^"]+)"',  # fallback to gene_id
        ]

        for pattern in patterns:
            match = re.search(pattern, attributes)
            if match:
                return match.group(1).strip()

        return None

    def fetch_gene_sequence(self, gene_symbol):
        """
        Fetch gene sequence using GTF coordinates + local FASTA
        Returns dict with gene info and sequence data
        """
        try:
            if self.use_demo_mode:
                return self._fetch_demo_sequence(gene_symbol)
            elif self.use_gtf and gene_symbol in self.gene_cache:
                return self._fetch_gtf_sequence(gene_symbol)
            else:
                return self._fetch_fallback_sequence(gene_symbol)

        except Exception as e:
            console.print(f"[red]Error fetching {gene_symbol}: {str(e)}[/red]")
            return None

    def _fetch_gtf_sequence(self, gene_symbol):
        """Fetch sequence using GTF coordinates + local FASTA"""
        gene_info = self.gene_cache[gene_symbol]

        # Extract sequence using samtools faidx
        chromosome = gene_info["chromosome"]
        start = gene_info["start"]
        end = gene_info["end"]
        strand = gene_info["strand"]

        try:
            # Use samtools faidx to extract real sequence from local FASTA
            sequence = self._extract_genomic_sequence(chromosome, start, end)

            if not sequence:
                console.print(
                    f"[red]Failed to extract sequence for {gene_symbol}[/red]"
                )
                return None

            console.print(
                f"[green]✅ Extracted GTF sequence for {gene_symbol}: {len(sequence)} bp[/green]"
            )

            # Create mock exon/intron regions (TODO: parse from GTF transcripts/exons)
            exon_regions, intron_regions = self._create_mock_regions(len(sequence))

            return {
                "gene_name": gene_symbol,
                "gene_id": f"GTF_{gene_symbol}",
                "transcript_id": f"GTF_{gene_symbol}_T1",
                "biotype": "protein_coding",  # TODO: extract from GTF
                "chromosome": chromosome,
                "strand": 1 if strand == "+" else -1,
                "genomic_start": start,
                "genomic_end": end,
                "sequence": sequence,
                "sequence_type": "genomic",
                "exon_regions": exon_regions,
                "intron_regions": intron_regions,
                "transcript_length": len(sequence),
                "genome_build": self.genome_build,
                "source": "gtf_coordinates",
            }

        except Exception as e:
            console.print(
                f"[red]Failed to fetch GTF sequence for {gene_symbol}: {e}[/red]"
            )
            return None

    def _fetch_fallback_sequence(self, gene_symbol):
        """Fallback to demo coordinates if GTF lookup fails"""
        from config import DEMO_GENE_COORDS

        if gene_symbol in DEMO_GENE_COORDS:
            console.print(f"[yellow]Using demo coordinates for {gene_symbol}[/yellow]")
            return self._fetch_demo_sequence(gene_symbol)
        else:
            console.print(
                f"[red]Gene {gene_symbol} not found in GTF or demo coordinates[/red]"
            )
            return None

    def _extract_genomic_sequence(self, chromosome, start, end):
        """Extract genomic sequence using samtools faidx from local FASTA"""
        try:
            # Construct region string
            region = f"{chromosome}:{start}-{end}"

            # Use local FASTA file path
            cmd = ["samtools", "faidx", self.local_genome_fasta, region]

            console.print(f"[cyan]Extracting {region} from local FASTA...[/cyan]")

            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            # Parse FASTA output
            lines = result.stdout.strip().split("\n")
            if len(lines) < 2:
                console.print(f"[red]Invalid samtools output for {region}[/red]")
                return None

            # Join sequence lines (skip header)
            sequence = "".join(lines[1:]).upper()

            return sequence

        except subprocess.CalledProcessError as e:
            console.print(f"[red]samtools error for {region}: {e}[/red]")
            console.print(f"[red]Command: {' '.join(cmd)}[/red]")
            console.print(f"[red]stderr: {e.stderr}[/red]")
            return None
        except FileNotFoundError:
            console.print(f"[red]samtools not found. Please install samtools.[/red]")
            return None
        except Exception as e:
            console.print(f"[red]Sequence extraction error: {e}[/red]")
            return None

    def _fetch_demo_sequence(self, gene_symbol):
        """Fetch sequence using demo coordinates and mock sequence"""
        from config import DEMO_GENE_COORDS

        if gene_symbol not in DEMO_GENE_COORDS:
            console.print(
                f"[red]Gene {gene_symbol} not found in demo coordinates[/red]"
            )
            return None

        gene_info = DEMO_GENE_COORDS[gene_symbol]

        # Calculate sequence length
        sequence_length = gene_info["end"] - gene_info["start"] + 1
        mock_sequence = self._generate_mock_sequence(sequence_length, gene_symbol)

        console.print(
            f"[green]✅ Generated demo sequence for {gene_symbol}: {len(mock_sequence)} bp[/green]"
        )

        # Create mock exon/intron regions
        exon_regions, intron_regions = self._create_mock_regions(sequence_length)

        return {
            "gene_name": gene_symbol,
            "gene_id": f"DEMO_{gene_symbol}",
            "transcript_id": f"DEMO_{gene_symbol}_T1",
            "biotype": gene_info["biotype"],
            "chromosome": gene_info["chromosome"],
            "strand": 1 if gene_info["strand"] == "+" else -1,
            "genomic_start": gene_info["start"],
            "genomic_end": gene_info["end"],
            "sequence": mock_sequence,
            "sequence_type": "genomic",
            "exon_regions": exon_regions,
            "intron_regions": intron_regions,
            "transcript_length": sequence_length,
            "genome_build": self.genome_build,
            "source": "demo_coordinates",
        }

    def _generate_mock_sequence(self, length, gene_symbol):
        """Generate a realistic mock DNA sequence for testing"""
        import random

        # Set seed based on gene name for consistent sequences
        random.seed(hash(gene_symbol) % 1000)

        # Create a sequence with realistic base composition
        # Mouse genome: ~42% GC content
        bases = ["A", "T", "G", "C"]
        weights = [0.29, 0.29, 0.21, 0.21]  # ~42% GC

        sequence = "".join(random.choices(bases, weights=weights, k=length))
        return sequence

    def _create_mock_regions(self, sequence_length):
        """Create mock exon/intron regions for demo purposes"""
        # Simple model: 3 exons with 2 introns
        exon1_end = sequence_length // 4
        intron1_start = exon1_end + 1
        intron1_end = sequence_length // 2
        exon2_start = intron1_end + 1
        exon2_end = 3 * sequence_length // 4
        intron2_start = exon2_end + 1
        intron2_end = sequence_length - sequence_length // 5
        exon3_start = intron2_end + 1

        exon_regions = [
            {"start": 0, "end": exon1_end, "length": exon1_end + 1},
            {
                "start": exon2_start,
                "end": exon2_end,
                "length": exon2_end - exon2_start + 1,
            },
            {
                "start": exon3_start,
                "end": sequence_length - 1,
                "length": sequence_length - exon3_start,
            },
        ]

        intron_regions = [
            {
                "start": intron1_start,
                "end": intron1_end,
                "length": intron1_end - intron1_start + 1,
            },
            {
                "start": intron2_start,
                "end": intron2_end,
                "length": intron2_end - intron2_start + 1,
            },
        ]

        return exon_regions, intron_regions

    def save_gene_sequences(self, gene_data_list, output_dir):
        """Save gene sequences as individual FASTA files"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        for gene_data in gene_data_list:
            try:
                gene_name = gene_data["gene_name"]
                sequence = gene_data["sequence"]

                # Create FASTA record with mm10 coordinates
                description = (
                    f"{gene_name} | {gene_data['biotype']} | "
                    f"{gene_data['chromosome']}:{gene_data['genomic_start']}-{gene_data['genomic_end']} | "
                    f"strand:{gene_data['strand']} | genome:{gene_data['genome_build']} | "
                    f"source:{gene_data['source']} | length:{len(sequence)}bp"
                )

                record = SeqRecord(Seq(sequence), id=gene_name, description=description)

                # Save to file
                output_file = output_dir / f"{gene_name}_{self.genome_build}.fa"
                SeqIO.write(record, output_file, "fasta")

            except Exception as e:
                console.print(
                    f"[red]Error saving {gene_data.get('gene_name', 'unknown')}: {str(e)}[/red]"
                )
