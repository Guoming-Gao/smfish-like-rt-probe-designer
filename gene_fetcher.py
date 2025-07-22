# gene_fetcher.py - Ensembl API Integration for Gene Sequence Fetching

import requests
import time
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rich.console import Console

console = Console()


class GeneSequenceFetcher:
    """Fetch gene sequences from Ensembl with exon/intron annotation"""

    def __init__(self, config):
        self.config = config
        self.species = config["species"]
        self.ensembl_release = config.get("ensembl_release", "110")
        self.genome_build = config.get("genome_build", "GRCm39")
        self.transcript_selection = config.get("transcript_selection", "longest")

        # Ensembl REST API base URL
        self.base_url = "https://rest.ensembl.org"
        self.session = requests.Session()
        self.session.headers.update({"Content-Type": "application/json"})

    def fetch_gene_sequence(self, gene_symbol):
        """
        Fetch gene sequence with exon/intron annotation
        Returns dict with gene info and sequence data
        """
        try:
            # Step 1: Get gene information
            gene_info = self._get_gene_info(gene_symbol)
            if not gene_info:
                console.print(f"[red]Could not find gene: {gene_symbol}[/red]")
                return None

            # Step 2: Get transcript information
            transcript_info = self._get_transcript_info(gene_info["id"])
            if not transcript_info:
                console.print(
                    f"[red]Could not find transcripts for: {gene_symbol}[/red]"
                )
                return None

            # Step 3: Select transcript based on strategy
            selected_transcript = self._select_transcript(transcript_info)

            # Step 4: Get genomic sequence (includes introns for lncRNAs)
            genomic_sequence = self._get_genomic_sequence(gene_info)

            # Step 5: Annotate exon/intron regions
            exon_intron_regions = self._annotate_regions(selected_transcript, gene_info)

            return {
                "gene_name": gene_symbol,
                "gene_id": gene_info["id"],
                "transcript_id": selected_transcript["id"],
                "biotype": gene_info["biotype"],
                "chromosome": gene_info["seq_region_name"],
                "strand": gene_info["strand"],
                "genomic_start": gene_info["start"],
                "genomic_end": gene_info["end"],
                "sequence": genomic_sequence,  # Full genomic sequence
                "sequence_type": "genomic",  # Always genomic for intron coverage
                "exon_regions": exon_intron_regions["exons"],
                "intron_regions": exon_intron_regions["introns"],
                "transcript_length": selected_transcript.get(
                    "length", len(genomic_sequence)
                ),
            }

        except Exception as e:
            console.print(f"[red]Error fetching {gene_symbol}: {str(e)}[/red]")
            return None

    def _get_gene_info(self, gene_symbol):
        """Get gene information from Ensembl"""
        species_name = "mus_musculus" if self.species == "mouse" else "homo_sapiens"
        url = f"{self.base_url}/lookup/symbol/{species_name}/{gene_symbol}"

        try:
            response = self.session.get(url, params={"expand": "1"})
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            console.print(
                f"[yellow]Ensembl API error for {gene_symbol}: {str(e)}[/yellow]"
            )
            return None

    def _get_transcript_info(self, gene_id):
        """Get all transcripts for a gene"""
        url = f"{self.base_url}/lookup/id/{gene_id}"

        try:
            response = self.session.get(url, params={"expand": "1"})
            response.raise_for_status()
            gene_data = response.json()
            return gene_data.get("Transcript", [])
        except requests.exceptions.RequestException as e:
            console.print(f"[yellow]Could not get transcripts: {str(e)}[/yellow]")
            return []

    def _select_transcript(self, transcripts):
        """Select transcript based on strategy"""
        if not transcripts:
            return None

        if self.transcript_selection == "longest":
            return max(transcripts, key=lambda t: t.get("length", 0))
        elif self.transcript_selection == "canonical":
            canonical = [t for t in transcripts if t.get("is_canonical")]
            return canonical[0] if canonical else transcripts[0]
        else:  # 'all' or default - return first
            return transcripts[0]

    def _get_genomic_sequence(self, gene_info):
        """Get genomic sequence (includes introns)"""
        species_name = "mus_musculus" if self.species == "mouse" else "homo_sapiens"
        chromosome = gene_info["seq_region_name"]
        start = gene_info["start"]
        end = gene_info["end"]
        strand = gene_info["strand"]

        url = f"{self.base_url}/sequence/region/{species_name}/{chromosome}:{start}..{end}"
        params = {"strand": strand}

        try:
            response = self.session.get(url, params=params)
            response.raise_for_status()
            sequence_data = response.json()
            return sequence_data["seq"].upper()
        except requests.exceptions.RequestException as e:
            console.print(f"[yellow]Could not get genomic sequence: {str(e)}[/yellow]")
            return None

    def _annotate_regions(self, transcript, gene_info):
        """Annotate exon and intron regions in genomic coordinates"""
        exons = transcript.get("Exon", [])

        # Convert to relative coordinates (0-based, relative to gene start)
        gene_start = gene_info["start"]
        strand = gene_info["strand"]

        exon_regions = []
        for exon in exons:
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

        # Calculate introns (regions between exons)
        intron_regions = []
        if len(exon_regions) > 1:
            exon_regions.sort(key=lambda x: x["start"])
            for i in range(len(exon_regions) - 1):
                intron_start = exon_regions[i]["end"] + 1
                intron_end = exon_regions[i + 1]["start"] - 1

                if intron_end > intron_start:  # Valid intron
                    intron_regions.append(
                        {
                            "start": intron_start,
                            "end": intron_end,
                            "length": intron_end - intron_start + 1,
                        }
                    )

        return {"exons": exon_regions, "introns": intron_regions}

    def save_gene_sequences(self, gene_data_list, output_dir):
        """Save gene sequences as individual FASTA files"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        for gene_data in gene_data_list:
            try:
                gene_name = gene_data["gene_name"]
                sequence = gene_data["sequence"]

                # Create FASTA record
                description = (
                    f"{gene_name} | {gene_data['biotype']} | "
                    f"{gene_data['chromosome']}:{gene_data['genomic_start']}-{gene_data['genomic_end']} | "
                    f"strand:{gene_data['strand']} | length:{len(sequence)}bp"
                )

                record = SeqRecord(Seq(sequence), id=gene_name, description=description)

                # Save to file
                output_file = output_dir / f"{gene_name}_full_genomic.fa"
                SeqIO.write(record, output_file, "fasta")

            except Exception as e:
                console.print(
                    f"[red]Error saving {gene_data.get('gene_name', 'unknown')}: {str(e)}[/red]"
                )
