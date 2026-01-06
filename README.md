# smfish-like-rt-probe-designer

Design focused RT primers for SWIFT-seq using smFISH probe design principles (Oligostan), with SNP coverage analysis for B6xCast allelic expression studies.

## Quick Start

```bash
conda activate blast
python main.py --genes Nanog Pou5f1 Xist --output ./my_output
```

## Pipeline Overview

```
Gene Symbols → [GTF/FASTA] → [Oligostan dG37] → [Filters] → [SNP Analysis] → [Local BLAST] → Output
```

1. **Fetch sequences** from GTF + FASTA (supports gzipped refGene format)
2. **Design primers** (26-32 nt) using Oligostan thermodynamic optimization (dG37 = -32.0)
3. **Filter** by GC content (40-60%), PNAS rules, dustmasker
4. **Analyze SNPs** in RT coverage region (200 nt upstream of probe) using VCF with tabix
5. **Local BLAST** for specificity validation (requires mm10 BLAST database)
6. **Output** CSV and FASTA files with RTBC barcodes for SWIFT-seq

---

## RT Coverage and SNP Detection Logic

The key to understanding SNP coverage is the directionality of RT extension:

### For **+ strand** genes (e.g., Nanog on chr6):

```
                        gene direction →
                        mRNA 5'→3'
                        ─────────────────────────────────────────▶
Genomic DNA:  5'═══════════════════════════════════════════════════3'
                        ▲           ▲
                        │           │
              RT region │  probe    │
              (200 bp)  │  binds    │
              ◀─────────┤  here     │
                        │           │
RT extends this way ◀───┘           │
(3'→5' on mRNA)                     │
                                    ▼
                        ┌───────────┐
Probe (DNA):            │ 5'────3'  │  (complementary to mRNA)
                        └───────────┘
                              │
                              ▼ RT extends
                        ┌─────────────────────────┐
RT Product (cDNA):      │ 5'────────────────────3'│  (extends ~200bp UPSTREAM)
                        └─────────────────────────┘
                        ◀─────────────────────────▶
                        RT_Region_Start    Probe_End
                        (lower coords)     (higher coords)
```

**+ strand summary:**
- Probe binds at `Probe_Start..Probe_End` (higher genomic coordinates)
- RT extends **UPSTREAM** (toward lower genomic coordinates)
- `RT_Region = [Probe_Start - 200, Probe_Start - 1]`
- SNPs checked in RT region for allelic distinction

---

### For **- strand** genes (e.g., Xist on chrX):

```
                        ◀ gene direction
                          mRNA 5'→3'
                        ◀─────────────────────────────────────────
Genomic DNA:  5'═══════════════════════════════════════════════════3'
                                    ▲           ▲
                                    │           │
                                 probe      RT region
                                 binds      (200 bp)
                                 here       ─────────▶
                                    │
                                    └───▶ RT extends this way
                                          (3'→5' on mRNA)

                        ┌───────────┐
Probe (DNA):            │ 5'────3'  │  (complementary to mRNA)
                        └───────────┘
                              │
                              ▼ RT extends
                        ┌─────────────────────────┐
RT Product (cDNA):      │ 5'────────────────────3'│  (extends ~200bp DOWNSTREAM)
                        └─────────────────────────┘
                        ◀─────────────────────────▶
                        Probe_Start       RT_Region_End
                        (lower coords)    (higher coords)
```

**- strand summary:**
- Probe binds at `Probe_Start..Probe_End` (lower genomic coordinates)
- RT extends **DOWNSTREAM** (toward higher genomic coordinates)
- `RT_Region = [Probe_End + 1, Probe_End + 200]`
- SNPs checked in RT region for allelic distinction

---

### Key Insight

**RT always extends 3'→5' on the mRNA template**, which means:
- `+ strand gene`: RT goes toward **lower** genomic coordinates (upstream)
- `- strand gene`: RT goes toward **higher** genomic coordinates (downstream)

This is critical for correct SNP coverage calculation!

---

## Output Files

| Step | File | Description |
|------|------|-------------|
| 1 | `FISH_RT_probes_PRE_BLAST_CANDIDATES.csv` | Intermediate: All high-quality candidates before BLAST |
| 1 | `FISH_RT_probes_FINAL_SELECTION.csv` | **Main Output**: Final probe set after all filters + BLAST |
| 1 | `FISH_RT_probes_FINAL_SELECTION.fasta` | FASTA format for synthesis (with RTBC if enabled) |
| 1 | `FISH_RT_probes_FINAL_SELECTION_summary.txt` | Design statistics and quality summary |

## Workflow

```bash
# Single command - runs entire pipeline including local BLAST
conda activate blast
python main.py --genes Nanog Xist Mecp2 --output ./results

# Or run with all 21 test genes
python main.py --test --output ./results
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
    "max_probes_per_gene": 200,         # Maximum probes per gene
    "min_snps_for_selection": 3,        # Min SNPs to include probe in selection
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
