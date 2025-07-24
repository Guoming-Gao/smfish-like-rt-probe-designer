# smfish-like-rt-probe-designer

Design FISH probes optimized for RT-PCR applications with comprehensive SNP coverage analysis for allelic expression studies.

## Overview

This tool combines the thermodynamic rigor of Oligostan probe design with:
- **Gene-based input** (gene symbols instead of FASTA files)
- **Strand-aware RT coverage analysis** (100nt downstream)
- **SNP profiling** for allelic expression studies
- **BLAST specificity filtering**
- **RTBC barcode integration** for SWIFT-seq workflows

## Quick Start

Install dependencies
```bash
pip install -r requirements.txt
```

Run with test gene set
```bash
python main.py --test --output ./test_output
```

Run with specific genes
```bash
python main.py --genes Nanog Pou5f1 Sox2 --output ./my_analysis
```

Run with custom gene list from config
```bash
python main.py --output ./full_analysis
```

## Complete Workflow: From Design to Final Probe Selection

### Step 1: Initial Probe Design
Run the main pipeline to generate high-quality probes with SNP coverage analysis:

```bash
python main.py --genes Nanog Pou5f1 Sox2 --output ./probe_design
```

This generates stringently filtered probes (GC + all PNAS rules + dustmasker + SNP coverage ≥2) in:
- `FISH_RT_probes_FILTERED.csv` - All high-quality probes
- `FISH_RT_probes_HIGH_SNP_2plus.csv` - Probes with high SNP coverage

### Step 2: Probe Selection (Interactive GUI)
Use the probe picker to select top probes per gene (default: 10 per gene):

```bash
python probe_picker_gui.py
```

1. **Load your CSV file** (e.g., `FISH_RT_probes_HIGH_SNP_2plus.csv`)
2. **Adjust parameters** (probes per gene, filters)
3. **Select top probes** - sorted by SNP coverage and PNAS score
4. **Export results** as `FISH_RT_probes_TOP10.csv` and `FISH_RT_probes_TOP10.fasta`

### Step 3: BLAST Specificity Check
Submit the selected probes to NCBI BLAST for specificity validation:

1. **Go to NCBI BLAST**: https://blast.ncbi.nlm.nih.gov/Blast.cgi
2. **Upload FASTA**: Use `FISH_RT_probes_TOP10.fasta` (probe sequences only, no RTBC)
3. **Database**: Select "Mouse genomic + transcript"
4. **Parameters**: Use default settings or adjust e-value as needed
5. **Download results**: Save as plain text format (`.txt` file)

### Step 4: BLAST Analysis & Final Selection
Analyze BLAST results to identify probes with unique genomic targets:

```bash
python top_probes_blast_analysis.py
```

1. **Select BLAST text file** from Step 3
2. **Select CSV file** (`FISH_RT_probes_TOP10.csv`) from Step 2
3. **Review analysis** - identifies probes with exactly 1 genomic hit
4. **Use final output**: `*_UNIQUE_HITS.csv` contains probes ready for synthesis

### Final Output
- **Experiment-ready probes**: Use `*_UNIQUE_HITS.csv` for probe synthesis
- **Quality metrics**: 100% identity, unique genomic targets, high SNP coverage
- **Synthesis sequences**: RTBC-containing sequences ready for ordering

## Features

- **Gene Symbol Input**: Automatic sequence fetching via local GTF/FASTA files
- **Thermodynamic Optimization**: Proven Oligostan algorithm (dG37 = -32.0)
- **Strand-Aware Coverage**: Precise 100nt downstream RT coverage calculation
- **SNP Analysis**: Coverage profiling for allelic expression studies
- **Stringent PNAS+SNP Filtering**: Keep only high-quality probes for BLAST verification
- **Interactive Probe Selection**: GUI-based probe picker for optimal subset selection
- **Batch Processing**: Process multiple genes efficiently
- **Comprehensive Output**: Consolidated CSV and FASTA files with all metrics

### PNAS Filter Rules

1. **Rule 1**: Adenine content 50% cytosine

## Configuration

Edit `config.py` to customize:

```python
FISH_RT_CONFIG = {
    'gene_list': ['Nanog', 'Pou5f1', 'Sox2'],  # Your genes
    'rt_coverage_downstream': 100,  # RT coverage length
    'rtbc_sequence': '/5Phos/TGACTTGAGGAT',  # RTBC barcode
    'local_gtf_path': '/path/to/annotations.gtf',  # Local GTF file
    'local_genome_fasta_path': '/path/to/genome.fa',  # Local FASTA file
    'snp_file_path': '/path/to/snps.txt',  # Your SNP file
    'min_snp_coverage_for_final': 2,  # Minimum SNPs for final output
}
```

## Requirements

- Python 3.7+
- samtools (for local genome sequence extraction)
- Local genome files (GTF + FASTA + SNP files)
- NCBI BLAST web interface (for specificity analysis)

## Citation

Based on the original Oligostan algorithm from:
> Tsanov, Nikolay, et al. "smiFISH and FISH-Quant – a Flexible Single RNA Detection Approach with Super-Resolution Capability." *Nucleic Acids Research* 44, no. 22 (2016): e165. https://doi.org/10.1093/NAR/GKW784

## License

MIT License - see LICENSE file for details.