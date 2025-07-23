# config.py - SIMPLIFIED (No BLAST, No Redundant Folders)
import os

FISH_RT_CONFIG = {
    # INPUT SETTINGS
    "gene_list": [
        "Atrx",
        "Diaph2",
        "Gpc4",
        "Hdac3",
        "Hnrnpu",
        "Kdm5c",
        "Kdm6a",
        "Mecp2",
        "Mid1",
        "Nanog",
        "Pir",
        "Pou5f1",
        "Rbmx",
        "Rlim",
        "Rps6ka3",
        "Rps6ka6",
        "Smc1a",
        "Spen",
        "Tsix",
        "Xist",
        "Zfp42",
    ],
    "transcript_selection": "longest",
    "output_directory": "/Users/gmgao/Dropbox/Caltech_PostDoc_GuttmanLab/constructs_and_smiFISH/smFISH_like_focusedRT-XCI",
    # LOCAL FILE PATHS (Simplified - direct paths)
    "genome_build": "mm10",
    "local_gtf_path": "/Volumes/guttman-1/annotations/mm10/Mus_musculus.GRCm38.96.gtf",
    "local_genome_fasta_path": "/Volumes/guttman-1/genomes/mm10/mm10withchr.fa",
    "snp_file_path": "/Volumes/guttman-1/data/snps/Bl6xCast.mm10.snps",
    # RT COVERAGE SETTINGS
    "rt_coverage_downstream": 100,  # nt downstream in RNA 5'→3' direction
    "include_probe_in_coverage": False,  # Strictly downstream
    # RTBC BARCODE SETTINGS
    "add_rtbc_barcode": True,
    "rtbc_sequence": "/5Phos/TGACTTGAGGAT",
    # OLIGOSTAN PARAMETERS
    "fixed_dg37_value": -32.0,
    "score_min": 0.9,
    "probe_length_max": 32,
    "probe_length_min": 26,
    "probe_spacing_min": 2,
    "gc_content_min": 0.4,
    "gc_content_max": 0.6,
    "pnas_filter_rules": [1, 2, 4],
    "salt_concentration": 0.115,
    # OPTIONAL FILTERS
    "use_dustmasker": False,
    "max_masked_percent": 0.1,
    # FASTA OUTPUT SETTINGS (NEW)
    "generate_blast_fasta": True,  # Generate FASTA for manual BLAST
    "fasta_description_format": "detailed",  # 'simple' or 'detailed'
}

# =============================================================================
# THERMODYNAMIC PARAMETERS (unchanged)
# =============================================================================

DG37_VALUES = {
    "AA": -0.2,
    "AC": -1.5,
    "AG": -0.9,
    "AT": -1.0,
    "CA": -1.0,
    "CC": -2.2,
    "CG": -1.2,
    "CT": -1.4,
    "GA": -0.8,
    "GC": -2.4,
    "GG": -1.5,
    "GT": -1.0,
    "TA": -0.3,
    "TC": -1.4,
    "TG": -1.0,
    "TT": -0.4,
}

# =============================================================================
# TEST GENE SETS (unchanged)
# =============================================================================

TEST_GENES_21 = [
    "Atrx",
    "Diaph2",
    "Gpc4",
    "Hdac3",
    "Hnrnpu",
    "Kdm5c",
    "Kdm6a",
    "Mecp2",
    "Mid1",
    "Nanog",
    "Pir",
    "Pou5f1",
    "Rbmx",
    "Rlim",
    "Rps6ka3",
    "Rps6ka6",
    "Smc1a",
    "Spen",
    "Tsix",
    "Xist",
    "Zfp42",
]

TEST_GENES_SMALL = ["Nanog", "Pou5f1", "Sox2"]

# =============================================================================
# DEMO GENE COORDINATES (unchanged)
# =============================================================================

DEMO_GENE_COORDS = {
    "Nanog": {
        "chromosome": "chr6",
        "start": 122707489,
        "end": 122714633,
        "strand": "+",
        "biotype": "protein_coding",
    },
    "Pou5f1": {
        "chromosome": "chr17",
        "start": 35509945,
        "end": 35514747,
        "strand": "+",
        "biotype": "protein_coding",
    },
    "Sox2": {
        "chromosome": "chr3",
        "start": 34650840,
        "end": 34652882,
        "strand": "-",
        "biotype": "protein_coding",
    },
}
