# gene_fetcher.py - Local File Gene Fetcher (NO ENSEMBL API)

import pandas as pd
import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rich.console import Console
import re
import os

console = Console()


class GeneSequenceFetcher:
    """Local gene fetcher using GTF + FASTA files or demo coordinates"""

    def __init__(self, config):
        self.config = config
        self.gtf_path = config.get("local_gtf_path", "")
        self.genome_fasta_path = config.get("local_genome_fasta_path", "")
        self.transcript_selection = config.get("transcript_selection", "longest")

        # Check if local files exist
        self.use_local_files = os.path.exists(self.gtf_path) and os.path.exists(
            self.genome_fasta_path
        )

        if self.use_local_files:
            console.print(f"[green]✅ Using local files:[/green]")
            console.print(f"[cyan]  GTF: {self.gtf_path}[/cyan]")
            console.print(f"[cyan]  FASTA: {self.genome_fasta_path}[/cyan]")
            # Load GTF data
            self.gene_cache = {}
            self._load_gtf_data()
        else:
            console.print(
                f"[yellow]⚠️  Local files not accessible, using demo coordinates[/yellow]"
            )
            console.print(f"[yellow]  GTF: {self.gtf_path}[/yellow]")
            console.print(f"[yellow]  FASTA: {self.genome_fasta_path}[/yellow]")

    def _load_gtf_data(self):
        """Load and parse GTF file"""
        console.print(f"[cyan]Loading GTF data...[/cyan]")

        try:
            genes = {}
            transcripts = {}
            exons = {}

            with open(self.gtf_path, "r") as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue

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

                    # Parse attributes
                    attr_dict = self._parse_gtf_attributes(attributes)
                    gene_name = attr_dict.get("gene_name", "").strip('"')

                    if not gene_name:
                        continue

                    # Store gene information
                    if feature == "gene":
                        genes[gene_name] = {
                            "chromosome": chrom,
                            "start": int(start),
                            "end": int(end),
                            "strand": 1 if strand == "+" else -1,
                            "gene_id": attr_dict.get("gene_id", "").strip('"'),
                            "biotype": attr_dict.get("gene_biotype", "").strip('"'),
                        }

                    # Store transcript information
                    elif feature == "transcript":
                        transcript_id = attr_dict.get("transcript_id", "").strip('"')
                        if transcript_id and gene_name:
                            if gene_name not in transcripts:
                                transcripts[gene_name] = []

                            transcripts[gene_name].append(
                                {
                                    "transcript_id": transcript_id,
                                    "start": int(start),
                                    "end": int(end),
                                    "length": int(end) - int(start) + 1,
                                    "is_canonical": "basic" in attributes,
                                }
                            )

                    # Store exon information
                    elif feature == "exon":
                        transcript_id = attr_dict.get("transcript_id", "").strip('"')
                        if transcript_id and gene_name:
                            if gene_name not in exons:
                                exons[gene_name] = {}
                            if transcript_id not in exons[gene_name]:
                                exons[gene_name][transcript_id] = []

                            exons[gene_name][transcript_id].append(
                                {"start": int(start), "end": int(end)}
                            )

            self.gene_cache = {
                "genes": genes,
                "transcripts": transcripts,
                "exons": exons,
            }

            console.print(f"[green]✅ Loaded {len(genes)} genes from GTF[/green]")

        except Exception as e:
            console.print(f"[red]Error loading GTF: {e}[/red]")
            self.use_local_files = False

    def _parse_gtf_attributes(self, attributes):
        """Parse GTF attributes string"""
        attr_dict = {}
        for attr in attributes.split(";"):
            attr = attr.strip()
            if " " in attr:
                parts = attr.split(" ", 1)
                if len(parts) == 2:
                    key, value = parts
                    attr_dict[key] = value.strip('"')
        return attr_dict

    def fetch_gene_sequence(self, gene_symbol):
        """Fetch gene sequence from local files or demo coordinates"""
        try:
            if self.use_local_files:
                return self._fetch_from_local_files(gene_symbol)
            else:
                return self._fetch_from_demo_coords(gene_symbol)
        except Exception as e:
            console.print(f"[red]Error fetching {gene_symbol}: {str(e)}[/red]")
            return None

    def _fetch_from_local_files(self, gene_symbol):
        """Fetch from local GTF + FASTA files"""
        genes = self.gene_cache["genes"]

        if gene_symbol not in genes:
            console.print(f"[red]Gene {gene_symbol} not found in GTF[/red]")
            return None

        gene_info = genes[gene_symbol]

        # Get transcript information
        transcripts = self.gene_cache["transcripts"].get(gene_symbol, [])
        if not transcripts:
            console.print(f"[yellow]No transcripts found for {gene_symbol}[/yellow]")
            return None

        # Select transcript
        selected_transcript = self._select_transcript(transcripts)

        # Extract genomic sequence using samtools
        genomic_sequence = self._extract_sequence_samtools(
            gene_info["chromosome"], gene_info["start"], gene_info["end"]
        )

        if not genomic_sequence:
            console.print(f"[red]Failed to extract sequence for {gene_symbol}[/red]")
            return None

        # Get exon/intron regions
        exon_intron_regions = self._get_exon_intron_regions(
            gene_symbol, selected_transcript, gene_info
        )

        return {
            "gene_name": gene_symbol,
            "gene_id": gene_info["gene_id"],
            "transcript_id": selected_transcript["transcript_id"],
            "biotype": gene_info["biotype"],
            "chromosome": gene_info["chromosome"],
            "strand": gene_info["strand"],
            "genomic_start": gene_info["start"],
            "genomic_end": gene_info["end"],
            "sequence": genomic_sequence,
            "sequence_type": "genomic",
            "exon_regions": exon_intron_regions["exons"],
            "intron_regions": exon_intron_regions["introns"],
            "transcript_length": selected_transcript["length"],
            "genome_build": "GRCm38",
            "source": "local_files",
        }

    def _fetch_from_demo_coords(self, gene_symbol):
        """Fetch using demo coordinates"""
        from config import DEMO_GENE_COORDS

        if gene_symbol not in DEMO_GENE_COORDS:
            console.print(f"[red]Gene {gene_symbol} not in demo coordinates[/red]")
            return None

        gene_info = DEMO_GENE_COORDS[gene_symbol]
        sequence_length = gene_info["end"] - gene_info["start"] + 1

        # Generate mock sequence
        mock_sequence = self._generate_mock_sequence(sequence_length, gene_symbol)

        # Create mock regions
        exon_regions, intron_regions = self._create_mock_regions(sequence_length)

        console.print(
            f"[yellow]Generated demo sequence for {gene_symbol}: {len(mock_sequence)} bp[/yellow]"
        )

        return {
            "gene_name": gene_symbol,
            "gene_id": f"DEMO_{gene_symbol}",
            "transcript_id": f"DEMO_{gene_symbol}_T1",
            "biotype": gene_info["biotype"],
            "chromosome": gene_info["chromosome"],
            "strand": gene_info["strand"],
            "genomic_start": gene_info["start"],
            "genomic_end": gene_info["end"],
            "sequence": mock_sequence,
            "sequence_type": "genomic",
            "exon_regions": exon_regions,
            "intron_regions": intron_regions,
            "transcript_length": sequence_length,
            "genome_build": "GRCm38",
            "source": "demo_coordinates",
        }

    def _extract_sequence_samtools(self, chromosome, start, end):
        """Extract sequence using samtools faidx"""
        try:
            region = f"{chromosome}:{start}-{end}"
            cmd = ["samtools", "faidx", self.genome_fasta_path, region]

            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            # Parse FASTA output
            lines = result.stdout.strip().split("\n")
            if len(lines) < 2:
                return None

            sequence = "".join(lines[1:]).upper()
            console.print(f"[green]Extracted {len(sequence)} bp for {region}[/green]")

            return sequence

        except subprocess.CalledProcessError as e:
            console.print(f"[red]samtools error: {e.stderr}[/red]")
            return None
        except Exception as e:
            console.print(f"[red]Sequence extraction error: {e}[/red]")
            return None

    def _select_transcript(self, transcripts):
        """Select transcript based on strategy"""
        if not transcripts:
            return None

        if self.transcript_selection == "longest":
            return max(transcripts, key=lambda t: t.get("length", 0))
        elif self.transcript_selection == "canonical":
            canonical = [t for t in transcripts if t.get("is_canonical")]
            return canonical[0] if canonical else transcripts[0]
        else:
            return transcripts[0]

    def _get_exon_intron_regions(self, gene_symbol, selected_transcript, gene_info):
        """Get exon/intron regions relative to gene start"""
        transcript_id = selected_transcript["transcript_id"]
        gene_start = gene_info["start"]

        exon_data = self.gene_cache["exons"].get(gene_symbol, {}).get(transcript_id, [])

        exon_regions = []
        for exon in exon_data:
            rel_start = exon["start"] - gene_start
            rel_end = exon["end"] - gene_start

            exon_regions.append(
                {
                    "start": rel_start,
                    "end": rel_end,
                    "length": rel_end - rel_start + 1,
                    "genomic_start": exon["start"],
                    "genomic_end": exon["end"],
                }
            )

        # Calculate introns
        intron_regions = []
        if len(exon_regions) > 1:
            exon_regions.sort(key=lambda x: x["start"])
            for i in range(len(exon_regions) - 1):
                intron_start = exon_regions[i]["end"] + 1
                intron_end = exon_regions[i + 1]["start"] - 1

                if intron_end > intron_start:
                    intron_regions.append(
                        {
                            "start": intron_start,
                            "end": intron_end,
                            "length": intron_end - intron_start + 1,
                        }
                    )

        return {"exons": exon_regions, "introns": intron_regions}

    def _generate_mock_sequence(self, length, gene_symbol):
        """Generate realistic mock sequence for demo"""
        import random

        random.seed(hash(gene_symbol) % 1000)

        bases = ["A", "T", "G", "C"]
        weights = [0.29, 0.29, 0.21, 0.21]  # ~42% GC

        return "".join(random.choices(bases, weights=weights, k=length))

    def _create_mock_regions(self, sequence_length):
        """Create mock exon/intron regions"""
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
        """Save gene sequences as FASTA files"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        for gene_data in gene_data_list:
            try:
                gene_name = gene_data["gene_name"]
                sequence = gene_data["sequence"]

                description = (
                    f"{gene_name} | {gene_data['biotype']} | "
                    f"chr{gene_data['chromosome']}:{gene_data['genomic_start']}-{gene_data['genomic_end']} | "
                    f"strand:{gene_data['strand']} | genome:{gene_data['genome_build']} | "
                    f"source:{gene_data['source']} | length:{len(sequence)}bp"
                )

                record = SeqRecord(Seq(sequence), id=gene_name, description=description)

                output_file = output_dir / f"{gene_name}_mm10.fa"
                SeqIO.write(record, output_file, "fasta")

            except Exception as e:
                console.print(
                    f"[red]Error saving {gene_data.get('gene_name', 'unknown')}: {str(e)}[/red]"
                )
