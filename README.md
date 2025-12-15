# smfish-like-rt-probe-designer

Design focused RT primers for SWIFT-seq using smFISH probe design principles, with comprehensive SNP coverage analysis for allelic expression studies.

## Overview

This tool applies the thermodynamic rigor of the Oligostan algorithm (originally developed for smFISH probe design) to design **focused RT primers for SWIFT-seq workflows**. By leveraging proven oligonucleotide design principles, it generates high-quality primers optimized for reverse transcription.

Key features:
- **Gene-based input** (gene symbols instead of FASTA files)
- **Oligostan thermodynamic optimization** (dG37 targeting for optimal hybridization)
- **Strand-aware RT coverage analysis** (100 nt downstream of primer binding site)
- **SNP profiling** for allelic expression studies (B6xCast mouse strains)
- **BLAST specificity filtering** for unique genomic targets
- **RTBC barcode integration** for SWIFT-seq compatibility

## How It Works

The pipeline uses smFISH probe design principles to create RT primers:

1. **Thermodynamic scoring**: Calculates dG37 (Gibbs free energy at 37°C) using nearest-neighbor dinucleotide values, targeting optimal dG = -32.0
2. **Quality filtering**: Applies GC content (40-60%), PNAS biochemical rules, and dustmasker repeat filtering
3. **RT coverage analysis**: For each primer, calculates the 100 nt region downstream (in RNA 5'→3' direction) that will be reverse-transcribed
4. **SNP profiling**: Identifies SNPs within the RT coverage region for allelic expression detection

## Quick Start

Install dependencies:
```bash
pip install -r requirements.txt
```

Run with test gene set:
```bash
python main.py --test --output ./test_output
```

Run with specific genes:
```bash
python main.py --genes Nanog Pou5f1 Sox2 --output ./my_analysis
```

Run with custom gene list from config:
```bash
python main.py --output ./full_analysis
```

## Complete Workflow: From Design to Final Primer Selection

### Step 1: Initial Primer Design
Run the main pipeline to generate high-quality RT primers with SNP coverage analysis:

```bash
python main.py --genes Nanog Pou5f1 Sox2 --output ./primer_design
```

This generates stringently filtered primers (GC + PNAS rules + dustmasker + SNP coverage ≥2) in:
- `FISH_RT_probes_FILTERED.csv` - All high-quality primers
- `FISH_RT_probes_HIGH_SNP_2plus.csv` - Primers with high SNP coverage

### Step 2: Primer Selection (Interactive GUI)
Use the primer picker to select top primers per gene (default: 10 per gene):

```bash
python probe_picker_gui.py
```

1. **Load your CSV file** (e.g., `FISH_RT_probes_HIGH_SNP_2plus.csv`)
2. **Adjust parameters** (primers per gene, filters)
3. **Select top primers** - sorted by SNP coverage and PNAS score
4. **Export results** as `FISH_RT_probes_TOP10.csv` and `FISH_RT_probes_TOP10.fasta`

### Step 3: BLAST Specificity Check
Submit the selected primers to NCBI BLAST for specificity validation:

1. **Go to NCBI BLAST**: https://blast.ncbi.nlm.nih.gov/Blast.cgi
2. **Upload FASTA**: Use `FISH_RT_probes_TOP10.fasta` (primer sequences only, no RTBC)
3. **Database**: Select "Mouse genomic + transcript"
4. **Parameters**: Use default settings or adjust e-value as needed
5. **Download results**: Save as plain text format (`.txt` file)

### Step 4: BLAST Analysis & Final Selection
Analyze BLAST results to identify primers with unique genomic targets:

```bash
python top_probes_blast_analysis.py
```

1. **Select BLAST text file** from Step 3
2. **Select CSV file** (`FISH_RT_probes_TOP10.csv`) from Step 2
3. **Review analysis** - identifies primers with exactly 1 genomic hit
4. **Use final output**: `*_UNIQUE_HITS.csv` contains primers ready for synthesis

### Final Output
- **Experiment-ready primers**: Use `*_UNIQUE_HITS.csv` for primer synthesis
- **Quality metrics**: 100% identity, unique genomic targets, high SNP coverage
- **Synthesis sequences**: RTBC-containing sequences ready for ordering

## Pipeline Architecture

```
Gene Symbols (e.g., Nanog, Pou5f1, Sox2)
    ↓
[Gene Fetcher] → Parse GTF, extract sequences from FASTA
    ↓
[Oligostan Algorithm] → Design primers (26-32 nt) with dG37 optimization
    ↓
[Stringent Filtering] → GC content, PNAS rules, dustmasker
    ↓
[SNP Analysis] → Calculate RT coverage region, find overlapping SNPs
    ↓
[Output] → CSV files, FASTA files (with/without RTBC barcodes)
```

## Features

- **Gene Symbol Input**: Automatic sequence fetching via local GTF/FASTA files
- **Thermodynamic Optimization**: Oligostan algorithm targeting dG37 = -32.0
- **Strand-Aware RT Coverage**: Calculates downstream region in RNA 5'→3' direction
- **SNP Analysis**: Coverage profiling for B6xCast allelic expression studies
- **Stringent Filtering**: GC content + PNAS biochemical rules + dustmasker
- **Interactive Selection**: GUI-based primer picker for optimal subset selection
- **RTBC Barcode Integration**: `/5Phos/TGACTTGAGGAT` prefix for SWIFT-seq
- **Batch Processing**: Process multiple genes efficiently

### PNAS Filter Rules

Biochemical quality rules to ensure primer synthesis and hybridization:
1. **Rule 1**: Adenine content < 28%
2. **Rule 2**: No AAAA homopolymer runs
3. **Rule 3**: Cytosine content 22-28%
4. **Rule 4**: No CCCC homopolymer runs
5. **Rule 5**: No 6-nt window with >50% cytosine

## Configuration

Edit `config.py` to customize:

```python
FISH_RT_CONFIG = {
    'gene_list': ['Nanog', 'Pou5f1', 'Sox2'],  # Your genes
    'rt_coverage_downstream': 100,  # RT coverage length (nt)
    'rtbc_sequence': '/5Phos/TGACTTGAGGAT',  # RTBC barcode for SWIFT-seq
    'local_gtf_path': '/path/to/annotations.gtf',  # Local GTF file
    'local_genome_fasta_path': '/path/to/genome.fa',  # Local FASTA file
    'snp_file_path': '/path/to/snps.txt',  # SNP file (B6xCast)
    'min_snp_coverage_for_final': 2,  # Minimum SNPs for final output
    'fixed_dg37_value': -32.0,  # Target Gibbs free energy
    'probe_length_min': 26,  # Minimum primer length
    'probe_length_max': 32,  # Maximum primer length
}
```

## Requirements

- Python 3.7+
- Local genome files (GTF + FASTA + SNP files for mm10)
- NCBI BLAST web interface (for specificity analysis)

### Python Dependencies
- biopython
- pandas
- numpy
- rich (for progress display)
- tkinter (for GUI)

## File Descriptions

| File | Purpose |
|------|---------|
| `main.py` | Main pipeline orchestration |
| `config.py` | Configuration parameters and gene lists |
| `gene_fetcher.py` | Extracts gene sequences from GTF/FASTA |
| `snp_analyzer.py` | SNP coverage analysis for RT regions |
| `output_generator.py` | CSV/FASTA output generation |
| `probe_picker_gui.py` | Interactive GUI for primer selection |
| `top_probes_blast_analysis.py` | BLAST result parsing and filtering |
| `extract_rt_sequences.py` | RT region sequence extraction utility |
| `utils/oligostan_core.py` | Core Oligostan algorithm |
| `utils/thermodynamics.py` | dG37 calculations |
| `utils/filters.py` | GC, PNAS, dustmasker filters |

## Citation

The thermodynamic design principles are based on the Oligostan algorithm from:
> Tsanov, Nikolay, et al. "smiFISH and FISH-Quant – a Flexible Single RNA Detection Approach with Super-Resolution Capability." *Nucleic Acids Research* 44, no. 22 (2016): e165. https://doi.org/10.1093/NAR/GKW784

## License

MIT License - see LICENSE file for details.
