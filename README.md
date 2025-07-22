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

Install dependencies\
```pip install -r requirements.txt```

Run with test gene set\
```python main.py --test --output ./test_output```

Run with specific genes\
```python main.py --genes Nanog Pou5f1 Sox2 --output ./my_analysis```

Run with custom gene list from config\
```python main.py --output ./full_analysis```



## Features

- **Gene Symbol Input**: Automatic sequence fetching via Ensembl API
- **Thermodynamic Optimization**: Proven Oligostan algorithm (dG37 = -32.0)
- **Strand-Aware Coverage**: Precise 100nt downstream RT coverage calculation
- **SNP Analysis**: Coverage profiling for allelic expression studies
- **BLAST Filtering**: Remove non-specific probes automatically
- **Batch Processing**: Process multiple genes efficiently
- **Comprehensive Output**: Consolidated CSV with all metrics

## Configuration

Edit `config.py` to customize:

```
FISH_RT_CONFIG = {
'gene_list': ['Nanog', 'Pou5f1', 'Sox2'], # Your genes
'species': 'mouse', # Currently mouse only
'rt_coverage_downstream': 100, # RT coverage length
'rtbc_sequence': '/5Phos/TGACTTGAGGAT', # RTBC barcode
'local_snp_file1': '/path/to/snps1.txt', # Your SNP files
'blast_database_path': '/path/to/blast_db' # BLAST database
}
```


## Requirements

- Python 3.7+
- BLAST+ suite (for specificity analysis)
- Internet connection (for Ensembl API)
- Custom SNP files (for allelic analysis)

## Citation

Based on the original Oligostan algorithm from:
> Tsanov, Nikolay, et al. "smiFISH and FISH-Quant – a Flexible Single RNA Detection Approach with Super-Resolution Capability." *Nucleic Acids Research* 44, no. 22 (2016): e165. https://doi.org/10.1093/NAR/GKW784

## License

MIT License - see LICENSE file for details.