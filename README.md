# smfish-like-rt-probe-designer

Design focused RT primers for SWIFT-seq using smFISH probe design principles (Oligostan), with SNP coverage analysis for B6xCast allelic expression studies.

## Quick Start

```bash
pip install -r requirements.txt
python main.py --genes Nanog Pou5f1 Xist --output ./my_output
```

## Pipeline Overview

```
Gene Symbols → [GTF/FASTA] → [Oligostan dG37] → [Filters] → [VCF SNP Analysis] → Output
```

1. **Fetch sequences** from GTF + FASTA (supports gzipped refGene format)
2. **Design primers** (26-32 nt) using Oligostan thermodynamic optimization (dG37 = -32.0)
3. **Filter** by GC content (40-60%), PNAS rules, dustmasker
4. **Analyze SNPs** in RT coverage region (100 nt downstream) using VCF with tabix
5. **Output** CSV and FASTA files with RTBC barcodes for SWIFT-seq

## Output Files

| Step | File | Description |
|------|------|-------------|
| 1 | `FISH_RT_probes_FILTERED.csv` | All primers passing quality filters |
| 1 | `FISH_RT_probes_HIGH_SNP_5plus.csv` | Primers with ≥5 B6/Cast SNP differences |
| 1 | `FISH_RT_probe_sequences_for_BLAST.fasta` | For NCBI BLAST validation |
| 2 | `FISH_RT_probes_TOP10.csv` | Top N per gene (via `probe_picker_gui.py`) |
| 4 | `*_UNIQUE_HITS.csv` | Final synthesis-ready primers (via `top_probes_blast_analysis.py`) |

## Workflow

```bash
# Step 1: Design primers
python main.py --test --output ./results

# Step 2: Select top primers per gene (GUI)
python probe_picker_gui.py

# Step 3: BLAST validation (manual at NCBI)

# Step 4: Filter for unique hits
python top_probes_blast_analysis.py
```

## Configuration

Edit `config.py`:

```python
FISH_RT_CONFIG = {
    # Input files (mm10)
    "local_gtf_path": "/Volumes/guttman/genomes/mm10/annotation/mm10.refGene.gtf.gz",
    "local_genome_fasta_path": "/Volumes/guttman/genomes/mm10/fasta/mm10.fa",
    "snp_file_path": "/Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz",

    # VCF samples for allelic SNP detection
    "vcf_b6_sample": "C57BL_6NJ",
    "vcf_cast_sample": "CAST_EiJ",

    # Parameters
    "rt_coverage_downstream": 100,      # nt downstream for SNP detection
    "min_snp_coverage_for_final": 5,    # Minimum B6/Cast SNP differences
    "fixed_dg37_value": -32.0,          # Target Gibbs free energy
    "probe_length_min": 26,
    "probe_length_max": 32,
}
```

## Requirements

- Python 3.7+
- **pysam** (for VCF/tabix queries)
- biopython, pandas, numpy, rich
- samtools (for FASTA indexing)

**One-time setup** (if FASTA not indexed):
```bash
cd /Volumes/guttman/genomes/mm10/fasta/
gunzip -k mm10.fa.gz && samtools faidx mm10.fa
```

## File Formats Supported

| File | Format | Notes |
|------|--------|-------|
| GTF | Gzipped refGene or Ensembl | Auto-infers genes from transcripts if no "gene" features |
| FASTA | Indexed (`.fai` required) | Use `samtools faidx` to create index |
| SNP | VCF with tabix index (`.tbi`) | Filters for B6 vs Cast **different** genotypes only |

## Test Results (21 genes, mm10)

| Metric | Value |
|--------|-------|
| Total probes designed | 41,118 |
| After quality filtering | 12,603 (30.7%) |
| Probes with ≥5 SNPs | 129 |
| SNP analysis time | ~2 min (tabix) |

### Validation

Probes at identical genomic positions produce identical results:
- Same sequences, dG37 scores, and filter outcomes
- Example: Mecp2 @ chr74027166 → `ATCTTTCTCTCTGCCTTGTCTGCCTGCTCCC`, dG37=-42.978

### SNP Detection Comparison (vs old custom .snps format)

| SNP Count | Old Format | New VCF |
|-----------|------------|---------|
| ≥5 SNPs | 78 probes | **129 probes** (+65%) |

The VCF-based approach finds more high-SNP probes by properly filtering for B6 vs Cast genotype differences.

## Key Files

| File | Purpose |
|------|---------|
| `main.py` | Pipeline orchestration |
| `config.py` | All configuration parameters |
| `gene_fetcher.py` | GTF/FASTA parsing (gzip support) |
| `snp_analyzer.py` | VCF-based SNP analysis with pysam/tabix |
| `probe_picker_gui.py` | Interactive primer selection |
| `top_probes_blast_analysis.py` | BLAST result filtering |

## Citation

Thermodynamic design based on:
> Tsanov et al. "smiFISH and FISH-Quant" *Nucleic Acids Research* 44(22):e165, 2016. https://doi.org/10.1093/NAR/GKW784

## License

MIT License
