# config.py - smfish-like-rt-probe-designer Configuration (Local Genome Approach)
import os

# =============================================================================
# MAIN CONFIGURATION (LOCAL GENOME APPROACH)
# =============================================================================

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
    "transcript_selection": "longest",  # 'longest', 'canonical', 'all'
    "output_directory": "/Users/gmgao/Dropbox/Caltech_PostDoc_GuttmanLab/constructs_and_smiFISH/smFISH_like_focusedRT-XCI",
    # LOCAL GENOME & ANNOTATION SOURCES (NO MORE ENSEMBL API)
    "genome_build": "mm10",  # Simplified identifier
    "local_genome_directory": "",  # TO BE FILLED - path to local genome FASTA files
    "gene_annotation_file": "",  # TO BE FILLED - path to gene annotation (GTF/GFF)
    # SNP ANALYSIS (Single file approach)
    "snp_file_path": "/Volumes/guttman-1/data/snps/Bl6xCast.mm10.snps",
    # RT COVERAGE SETTINGS (STRAND-AWARE)
    "rt_coverage_downstream": 100,  # nt downstream in RNA 5'→3' direction
    "include_probe_in_coverage": False,  # Strictly downstream
    # RTBC BARCODE SETTINGS
    "add_rtbc_barcode": True,
    "rtbc_sequence": "/5Phos/TGACTTGAGGAT",  # Configurable RTBC with phosphorylation
    # BLAST SPECIFICITY
    "blast_database_path": "",  # TO BE FILLED - path to BLAST database
    "require_unique_hits_only": True,  # Remove multi-hit probes
    "run_blast_analysis": False,  # Disabled by default until DB configured
    # OLIGOSTAN PARAMETERS
    "fixed_dg37_value": -32.0,  # Fixed dG37 (no optimization)
    "score_min": 0.9,  # Minimum probe score threshold
    "probe_length_max": 32,  # Maximum probe length
    "probe_length_min": 26,  # Minimum probe length
    "probe_spacing_min": 2,  # Minimum spacing between probes
    "gc_content_min": 0.4,  # Minimum GC content
    "gc_content_max": 0.6,  # Maximum GC content
    "pnas_filter_rules": [1, 2, 4],  # PNAS composition rules to apply
    "salt_concentration": 0.115,  # Salt concentration for dG calculations
    # OPTIONAL FILTERS
    "use_dustmasker": False,  # Enable dustmasker repeat filtering
    "max_masked_percent": 0.1,  # Maximum allowed masked percentage
}

# =============================================================================
# THERMODYNAMIC PARAMETERS (from original Oligostan)
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
# TEST GENE SETS
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

TEST_GENES_SMALL = ["Nanog", "Pou5f1", "Sox2"]  # For quick testing

# =============================================================================
# DEMO GENE COORDINATES (mm10) - For testing when local files not available
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
