import os
import sys

# Ensure current directory is in path
sys.path.append(os.getcwd())

import pandas as pd
from gene_fetcher import GeneSequenceFetcher
from utils.sequence_utils import determine_region_type
from rich.console import Console

console = Console()

def test_comprehensive_fix():
    # Mock config
    config = {
        "local_gtf_path": "dummy.gtf",
        "local_genome_fasta_path": "dummy.fa",
        "transcript_selection": "longest",
        "aggregate_exons_for_labeling": True
    }

    fetcher = GeneSequenceFetcher(config)
    fetcher.use_local_files = True # Force mock logic

    # Mock sequence extraction to avoid samtools dependency
    fetcher._extract_sequence_samtools = lambda *args: "A" * 2000

    # -----------------------------------------------------------------
    # Test Case 1: Forward strand gene with splicing variants
    # -----------------------------------------------------------------
    # T1 (Short): Exon 1 [1000, 1100], Exon 2 [1500, 1600] -> Spliced: 202
    # T2 (Long):  Exon 1 [1000, 1100], Exon 2 [1200, 1300], Exon 3 [1500, 1600] -> Spliced: 303
    fetcher.gene_cache = {
        "genes": {
            "GeneF": {"chromosome": "chr1", "start": 1000, "end": 2000, "strand": 1, "gene_id": "F1", "biotype": "pc"}
        },
        "transcripts": {
            "GeneF": [
                {"transcript_id": "T1", "start": 1000, "end": 2000, "length": 1001, "spliced_length": 202},
                {"transcript_id": "T2", "start": 1000, "end": 2000, "length": 1001, "spliced_length": 303}
            ]
        },
        "exons": {
            "GeneF": {
                "T1": [{"start": 1000, "end": 1100}, {"start": 1500, "end": 1600}],
                "T2": [{"start": 1000, "end": 1100}, {"start": 1200, "end": 1300}, {"start": 1500, "end": 1600}]
            }
        }
    }

    console.print("\n[bold]Testing Transcript Selection (T2 should be picked as longest spliced)...[/bold]")
    gene_data = fetcher.fetch_gene_sequence("GeneF")
    if gene_data["transcript_id"] == "T2":
        console.print("[green]PASS: Selected T2 (spliced length 303)[/green]")
    else:
        console.print(f"[red]FAIL: Selected {gene_data['transcript_id']} (expected T2)[/red]")

    # rel_exons (T2): [1, 101], [201, 301], [501, 601]
    # rel_introns (T2): [102, 200], [302, 500]

    console.print("\n[bold]Testing 1-indexed Inclusive Boundaries (Forward)...[/bold]")
    cases = [
        (1, 30, "exon"),
        (101, 130, "exon-intron"), # Spans 101 (exon) and 102 (intron)
        (102, 200, "intron"),
        (201, 230, "exon"),
        (301, 302, "exon-intron"),
    ]

    for start, end, expected in cases:
        result = determine_region_type(start, end, gene_data["exon_regions"], gene_data["intron_regions"])
        if result == expected:
            console.print(f"[green]PASS: [{start}, {end}] -> {result}[/green]")
        else:
            console.print(f"[red]FAIL: [{start}, {end}] -> {result} (expected {expected})[/red]")

    # -----------------------------------------------------------------
    # Test Case 2: Reverse strand gene
    # -----------------------------------------------------------------
    # Genomic: chr1:5000-6000 (-)
    # Exons: [5000, 5100], [5900, 6000]
    # rel_exons: [1, 101] (6000-5900+1), [901, 1001] (6000-5000+1)
    fetcher.gene_cache = {
        "genes": {
            "GeneR": {"chromosome": "chr1", "start": 5000, "end": 6000, "strand": -1, "gene_id": "R1", "biotype": "pc"}
        },
        "transcripts": {
            "GeneR": [{"transcript_id": "TR1", "start": 5000, "end": 6000, "length": 1001, "spliced_length": 202}]
        },
        "exons": {
            "GeneR": {"TR1": [{"start": 5000, "end": 5100}, {"start": 5900, "end": 6000}]}
        }
    }

    gene_data_r = fetcher.fetch_gene_sequence("GeneR")
    console.print("\n[bold]Testing 1-indexed Inclusive Boundaries (Reverse)...[/bold]")

    # 101 is the end of first exon (RNA 5' end)
    cases_r = [
        (1, 101, "exon"),
        (101, 102, "exon-intron"),
        (102, 130, "intron"),
        (901, 1001, "exon")
    ]

    for start, end, expected in cases_r:
        result = determine_region_type(start, end, gene_data_r["exon_regions"], gene_data_r["intron_regions"])
        if result == expected:
            console.print(f"[green]PASS: [{start}, {end}] -> {result}[/green]")
        else:
            console.print(f"[red]FAIL: [{start}, {end}] -> {result} (expected {expected})[/red]")

    # -----------------------------------------------------------------
    # Test Case 3: Aggregate Exon Labeling
    # -----------------------------------------------------------------
    # GeneA has T1 [1, 100] exon and T2 [150, 250] exon
    # Selected T1 (longest spliced if we say T1 > T2, or just pick one)
    fetcher.gene_cache = {
        "genes": {
            "GeneA": {"chromosome": "chr1", "start": 1000, "end": 2000, "strand": 1, "gene_id": "A1", "biotype": "pc"}
        },
        "transcripts": {
            "GeneA": [
                {"transcript_id": "T1", "start": 1000, "end": 2000, "length": 1001, "spliced_length": 100},
                {"transcript_id": "T2", "start": 1000, "end": 2000, "length": 1001, "spliced_length": 100}
            ]
        },
        "exons": {
            "GeneA": {
                "T1": [{"start": 1000, "end": 1099}],
                "T2": [{"start": 1150, "end": 1249}]
            }
        }
    }

    gene_data_a = fetcher.fetch_gene_sequence("GeneA")
    # T1 is picked. T1's introns: [] (only 1 exon)
    # T1's exons: [[1, 100]]
    # all_known_exons: [[1, 100], [151, 250]]

    console.print("\n[bold]Testing Aggregate Exon Labeling...[/bold]")
    # Test with standard T1 regions
    res_std = determine_region_type(151, 180, gene_data_a["exon_regions"], gene_data_a["intron_regions"])
    console.print(f"Standard (T1) check for [151, 180]: {res_std} (expected intergenic/intron-like)")

    # Test with aggregate regions
    res_agg = determine_region_type(151, 180, gene_data_a["all_known_exons"], gene_data_a["intron_regions"])
    if res_agg == "exon":
        console.print("[green]PASS: Aggregate logic correctly labeled [151, 180] as exon[/green]")
    else:
        console.print(f"[red]FAIL: Aggregate logic labeled [151, 180] as {res_agg}[/red]")

if __name__ == "__main__":
    test_comprehensive_fix()
