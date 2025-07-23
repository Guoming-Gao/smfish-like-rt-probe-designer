#!/usr/bin/env python3
# debug_e102_api.py - Debug Ensembl Release 102 API thoroughly

import requests
import json
from rich.console import Console
from rich.table import Table

console = Console()


class TestGeneSequenceFetcher:
    """Exact copy of gene_fetcher logic for debugging"""

    def __init__(self):
        self.base_url = "https://e102.rest.ensembl.org"
        self.session = requests.Session()
        self.session.headers.update({"Content-Type": "application/json"})

        console.print(f"[cyan]Testing API endpoint: {self.base_url}[/cyan]")

    def test_gene_info(self, gene_symbol):
        """Test _get_gene_info() method exactly"""
        url = f"{self.base_url}/lookup/symbol/mus_musculus/{gene_symbol}"

        console.print(f"\n[bold]Testing Gene Info Fetch[/bold]")
        console.print(f"URL: {url}")

        try:
            params = {"expand": "1"}
            console.print(f"Params: {params}")
            console.print(f"Headers: {dict(self.session.headers)}")

            # Make request exactly like gene_fetcher.py
            response = self.session.get(url, params=params)

            console.print(f"Status Code: {response.status_code}")
            console.print(f"Response Headers: {dict(response.headers)}")
            console.print(
                f"Content-Type: {response.headers.get('content-type', 'Unknown')}"
            )
            console.print(f"Response Length: {len(response.text)}")

            # Show raw response first
            console.print(f"\n[yellow]Raw Response (first 500 chars):[/yellow]")
            console.print(response.text[:500])

            if response.status_code == 200:
                try:
                    gene_data = response.json()
                    console.print(f"\n[green]✅ JSON Parse Success![/green]")

                    # Show key fields
                    key_info = {
                        "Gene Name": gene_data.get("display_name", "Unknown"),
                        "Gene ID": gene_data.get("id", "Unknown"),
                        "Assembly": gene_data.get("assembly_name", "Unknown"),
                        "Chromosome": gene_data.get("seq_region_name", "Unknown"),
                        "Start": gene_data.get("start", "Unknown"),
                        "End": gene_data.get("end", "Unknown"),
                        "Strand": gene_data.get("strand", "Unknown"),
                    }

                    table = Table(title="Gene Information")
                    table.add_column("Field", style="cyan")
                    table.add_column("Value", style="green")

                    for key, value in key_info.items():
                        table.add_row(key, str(value))

                    console.print(table)
                    return True, gene_data

                except json.JSONDecodeError as e:
                    console.print(f"[red]❌ JSON Parse Error: {e}[/red]")
                    console.print(f"[red]Full Response: {response.text}[/red]")
                    return False, response.text
            else:
                console.print(f"[red]❌ HTTP Error: {response.status_code}[/red]")
                console.print(f"[red]Response: {response.text}[/red]")
                return False, response.text

        except requests.exceptions.RequestException as e:
            console.print(f"[red]❌ Request Exception: {e}[/red]")
            return False, str(e)

    def test_transcript_info(self, gene_id):
        """Test _get_transcript_info() method exactly"""
        url = f"{self.base_url}/lookup/id/{gene_id}"

        console.print(f"\n[bold]Testing Transcript Info Fetch[/bold]")
        console.print(f"URL: {url}")

        try:
            params = {"expand": "1"}
            response = self.session.get(url, params=params)

            console.print(f"Status Code: {response.status_code}")

            if response.status_code == 200:
                try:
                    gene_data = response.json()
                    transcripts = gene_data.get("Transcript", [])
                    console.print(
                        f"[green]✅ Found {len(transcripts)} transcripts[/green]"
                    )

                    if transcripts:
                        # Show first transcript
                        first_transcript = transcripts[0]
                        console.print(
                            f"First transcript ID: {first_transcript.get('id', 'Unknown')}"
                        )
                        console.print(
                            f"First transcript length: {first_transcript.get('length', 'Unknown')}"
                        )

                    return True, transcripts

                except json.JSONDecodeError as e:
                    console.print(f"[red]❌ Transcript JSON Parse Error: {e}[/red]")
                    return False, response.text
            else:
                console.print(
                    f"[red]❌ Transcript HTTP Error: {response.status_code}[/red]"
                )
                return False, response.text

        except requests.exceptions.RequestException as e:
            console.print(f"[red]❌ Transcript Request Exception: {e}[/red]")
            return False, str(e)

    def test_genomic_sequence(self, gene_info):
        """Test _get_genomic_sequence() method exactly"""
        chromosome = gene_info["seq_region_name"]
        start = gene_info["start"]
        end = gene_info["end"]
        strand = gene_info["strand"]

        url = (
            f"{self.base_url}/sequence/region/mus_musculus/{chromosome}:{start}..{end}"
        )

        console.print(f"\n[bold]Testing Genomic Sequence Fetch[/bold]")
        console.print(f"URL: {url}")

        try:
            params = {"strand": strand}
            response = self.session.get(url, params=params)

            console.print(f"Status Code: {response.status_code}")

            if response.status_code == 200:
                try:
                    sequence_data = response.json()
                    sequence = sequence_data.get("seq", "")
                    console.print(
                        f"[green]✅ Sequence fetched: {len(sequence)} bp[/green]"
                    )
                    console.print(f"First 100bp: {sequence[:100]}")
                    return True, sequence_data

                except json.JSONDecodeError as e:
                    console.print(f"[red]❌ Sequence JSON Parse Error: {e}[/red]")
                    return False, response.text
            else:
                console.print(
                    f"[red]❌ Sequence HTTP Error: {response.status_code}[/red]"
                )
                return False, response.text

        except requests.exceptions.RequestException as e:
            console.print(f"[red]❌ Sequence Request Exception: {e}[/red]")
            return False, str(e)


def test_alternative_headers():
    """Test different header combinations"""
    console.print(f"\n[bold yellow]Testing Different Header Combinations[/bold yellow]")

    base_url = "https://e102.rest.ensembl.org"
    gene_symbol = "Nanog"
    url = f"{base_url}/lookup/symbol/mus_musculus/{gene_symbol}"

    header_combinations = [
        {
            "name": "JSON headers (current)",
            "headers": {"Content-Type": "application/json"},
        },
        {"name": "No special headers", "headers": {}},
        {"name": "Accept JSON only", "headers": {"Accept": "application/json"}},
        {
            "name": "User-Agent + JSON",
            "headers": {
                "Content-Type": "application/json",
                "User-Agent": "smfish-rt-probe-designer/1.0",
            },
        },
    ]

    for combo in header_combinations:
        console.print(f"\n[cyan]Testing: {combo['name']}[/cyan]")
        console.print(f"Headers: {combo['headers']}")

        try:
            response = requests.get(
                url, params={"expand": "1"}, headers=combo["headers"], timeout=10
            )

            console.print(f"Status: {response.status_code}")
            console.print(
                f"Content-Type: {response.headers.get('content-type', 'Unknown')}"
            )
            console.print(f"Response length: {len(response.text)}")

            if response.status_code == 200:
                try:
                    data = response.json()
                    assembly = data.get("assembly_name", "Unknown")
                    console.print(f"[green]✅ Success! Assembly: {assembly}[/green]")
                except json.JSONDecodeError:
                    console.print(f"[red]❌ JSON parse failed[/red]")
                    console.print(f"Raw response: {response.text[:200]}")
            else:
                console.print(f"[red]❌ HTTP error: {response.text[:200]}[/red]")

        except Exception as e:
            console.print(f"[red]❌ Exception: {e}[/red]")


def test_network_connectivity():
    """Test basic network connectivity to different Ensembl endpoints"""
    console.print(f"\n[bold yellow]Testing Network Connectivity[/bold yellow]")

    test_urls = [
        "https://rest.ensembl.org",
        "https://e102.rest.ensembl.org",
        "https://e104.rest.ensembl.org",
        "https://www.ensembl.org",
    ]

    for url in test_urls:
        try:
            response = requests.get(url, timeout=5)
            console.print(f"[green]✅ {url}: {response.status_code}[/green]")
        except Exception as e:
            console.print(f"[red]❌ {url}: {e}[/red]")


def main():
    """Main debugging function"""
    console.print("[bold blue]🔍 Debug Ensembl Release 102 API[/bold blue]")
    console.print("=" * 60)

    # Test basic connectivity first
    test_network_connectivity()

    # Test different header combinations
    test_alternative_headers()

    # Test the exact gene_fetcher logic
    console.print(f"\n[bold yellow]Testing Exact Gene Fetcher Logic[/bold yellow]")
    console.print("-" * 50)

    fetcher = TestGeneSequenceFetcher()

    # Test gene info fetch
    success, gene_data = fetcher.test_gene_info("Nanog")

    if success:
        # Test transcript fetch
        gene_id = gene_data.get("id")
        if gene_id:
            transcript_success, transcripts = fetcher.test_transcript_info(gene_id)

        # Test sequence fetch
        sequence_success, sequence_data = fetcher.test_genomic_sequence(gene_data)

        console.print(f"\n[bold green]🎉 All tests passed for Nanog![/bold green]")
        console.print(f"Assembly: {gene_data.get('assembly_name')}")
        console.print(
            f"Coordinates: chr{gene_data.get('seq_region_name')}:{gene_data.get('start')}-{gene_data.get('end')}"
        )

    else:
        console.print(
            f"\n[bold red]❌ Gene fetch failed - this explains the main pipeline issue![/bold red]"
        )
        console.print(f"Error details: {gene_data}")


if __name__ == "__main__":
    main()
