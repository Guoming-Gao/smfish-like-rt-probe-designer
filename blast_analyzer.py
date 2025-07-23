# blast_analyzer.py - BLAST Specificity Analysis (Based on smFISH_blast_analysis.py)

import pandas as pd
import subprocess
import tempfile
import os
import re
from pathlib import Path
from rich.console import Console
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

console = Console()


class BLASTSpecificityAnalyzer:
    """BLAST specificity analysis using existing smFISH logic"""

    def __init__(self, config):
        self.config = config
        self.blast_db_path = config.get("blast_database_path", "")
        self.require_unique_only = config.get("require_unique_hits_only", True)

    def analyze_probes(self, probe_list, output_dir):
        """
        Run BLAST analysis on all probes and update with specificity information
        """
        if not self.blast_db_path:
            console.print(
                "[yellow]No BLAST database configured, skipping BLAST analysis[/yellow]"
            )
            return probe_list

        console.print(
            f"[cyan]Running BLAST analysis on {len(probe_list)} probes...[/cyan]"
        )

        try:
            # Step 1: Create FASTA file from probe sequences
            blast_input_file = self._create_probe_fasta(probe_list, output_dir)

            # Step 2: Run BLAST
            blast_output_file = self._run_blast(blast_input_file, output_dir)

            # Step 3: Parse BLAST results using smFISH logic
            blast_results_df = self._parse_blast_results(blast_output_file)

            # Step 4: Merge BLAST results with probe data
            updated_probes = self._merge_blast_results(probe_list, blast_results_df)

            console.print(f"✅ BLAST analysis completed")

            return updated_probes

        except Exception as e:
            console.print(f"[red]BLAST analysis failed: {str(e)}[/red]")
            # Return probes with default BLAST values
            for probe in probe_list:
                probe.update(
                    {
                        "NumberOfHits": 0,
                        "PercentAlignment": 0,
                        "UniqueHitName": "",
                        "BLAST_Start": 0,
                        "BLAST_End": 0,
                    }
                )
            return probe_list

    def _create_probe_fasta(self, probe_list, output_dir):
        """Create FASTA file from probe sequences"""
        fasta_file = Path(output_dir) / "probe_sequences.fasta"

        records = []
        for i, probe in enumerate(probe_list):
            # Use probe sequence as ID for matching later
            probe_id = f"probe_{i}_{probe['GeneName']}"  # Unique identifier
            sequence = probe["Seq"]

            record = SeqRecord(
                Seq(sequence),
                id=probe_id,
                description=f"Gene:{probe['GeneName']} Size:{probe['ProbeSize']}",
            )
            records.append(record)

        SeqIO.write(records, fasta_file, "fasta")
        console.print(f"Created FASTA file: {fasta_file}")

        return fasta_file

    def _run_blast(self, input_fasta, output_dir):
        """Run BLASTN against the configured database"""
        output_file = Path(output_dir) / "blast_results.txt"

        # BLAST command (text output format for parsing)
        blast_cmd = [
            "blastn",
            "-query",
            str(input_fasta),
            "-db",
            self.blast_db_path,
            "-out",
            str(output_file),
            "-outfmt",
            "0",  # Text format (same as smFISH script expects)
            "-max_target_seqs",
            "100",  # Allow up to 100 hits
            "-evalue",
            "0.01",
            "-word_size",
            "7",  # Sensitive for short sequences
        ]

        console.print(f"Running BLAST command: {' '.join(blast_cmd)}")

        try:
            result = subprocess.run(
                blast_cmd, capture_output=True, text=True, check=True
            )
            console.print(f"BLAST completed successfully")
            return output_file

        except subprocess.CalledProcessError as e:
            console.print(f"[red]BLAST failed with return code {e.returncode}[/red]")
            console.print(f"[red]STDERR: {e.stderr}[/red]")
            raise
        except FileNotFoundError:
            raise Exception(
                "BLAST not found. Please ensure BLAST+ is installed and in PATH"
            )

    def _parse_blast_results(self, blast_output_file):
        """
        Parse BLAST results using exact logic from smFISH_blast_analysis.py
        """
        try:
            with open(blast_output_file, "r") as f:
                blast_text = f.read()

            # Use exact parsing function from smFISH script
            blast_df = self.parse_blast_results_smfish(blast_text)

            console.print(f"Parsed {len(blast_df)} BLAST results")
            return blast_df

        except Exception as e:
            console.print(f"[red]Error parsing BLAST results: {str(e)}[/red]")
            return pd.DataFrame()

    def parse_blast_results_smfish(self, blast_text):
        """
        EXACT COPY of parse_blast_results from smFISH_blast_analysis.py
        """
        # Split by queries
        queries = re.split(r"Query #\d+: ", blast_text)[1:]  # skip first empty split
        records = []

        for query in queries:
            # Extract probe name (first word in the query)
            probe_name_search = re.search(r"^(\S+)", query)
            probe_name = probe_name_search.group(1) if probe_name_search else None

            # CORRECTED: Count alignment sections starting with ">"
            # This counts actual genomic hits, not table entries
            alignment_headers = re.findall(r"^>([^\n]+)", query, re.MULTILINE)
            num_hits = len(alignment_headers)

            # Extract unique hit name if only one hit
            unique_hit_name = None
            if num_hits == 1:
                unique_hit_name = alignment_headers[0].strip()

            # Extract start and end positions (from first alignment)
            start_end_search = re.search(r"Range 1: (\d+) to (\d+)", query)
            start = int(start_end_search.group(1)) if start_end_search else None
            end = int(start_end_search.group(2)) if start_end_search else None

            # Extract percentage identity
            perc_identity_search = re.search(r"Identities:\s*\d+/\d+\((\d+)%\)", query)
            perc_identity = (
                int(perc_identity_search.group(1)) if perc_identity_search else None
            )

            # Extract probe sequence (longest query sequence from alignment)
            query_seqs = re.findall(r"Query\s+\d+\s+([A-Z]+)\s+\d+", query)
            probe_seq = max(query_seqs, key=len) if query_seqs else None

            records.append(
                {
                    "ProbeName": probe_name,
                    "ProbeSequence": probe_seq,
                    "PercentAlignment": perc_identity,
                    "NumberOfHits": num_hits,
                    "UniqueHitName": unique_hit_name,
                    "Start": start,
                    "End": end,
                }
            )

        return pd.DataFrame(records)

    def _merge_blast_results(self, probe_list, blast_df):
        """Merge BLAST results back into probe list"""

        # Create mapping from probe index to BLAST results
        blast_lookup = {}
        for _, row in blast_df.iterrows():
            if row["ProbeName"]:
                # Extract probe index from probe name (probe_0_GeneName format)
                try:
                    probe_idx = int(row["ProbeName"].split("_")[1])
                    blast_lookup[probe_idx] = row
                except (IndexError, ValueError):
                    continue

        # Update probes with BLAST results
        for i, probe in enumerate(probe_list):
            if i in blast_lookup:
                blast_result = blast_lookup[i]
                probe.update(
                    {
                        "NumberOfHits": blast_result["NumberOfHits"] or 0,
                        "PercentAlignment": blast_result["PercentAlignment"] or 0,
                        "UniqueHitName": blast_result["UniqueHitName"] or "",
                        "BLAST_Start": blast_result["Start"] or 0,
                        "BLAST_End": blast_result["End"] or 0,
                    }
                )
            else:
                # No BLAST result found
                probe.update(
                    {
                        "NumberOfHits": 0,
                        "PercentAlignment": 0,
                        "UniqueHitName": "No_BLAST_result",
                        "BLAST_Start": 0,
                        "BLAST_End": 0,
                    }
                )

        return probe_list
