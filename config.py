# config.py - smfish-like-rt-probe-designer Configuration
import os

# =============================================================================
# MAIN CONFIGURATION
# =============================================================================

FISH_RT_CONFIG = {
    # INPUT SETTINGS
    "species": "mouse",  # Currently only 'mouse' supported
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
    "output_directory": "/path/to/output",  # Will be set by user
    # SEQUENCE & ANNOTATION SOURCES (TO BE FILLED FROM JAVA PATHS)
    "ensembl_release": "110",  # Match Java script if possible
    "genome_build": "GRCm39",
    "local_genome_directory": "",  # [your_args[3]_path] - TO BE FILLED
    "local_snp_file1": "",  # [your_args[1]_path] - TO BE FILLED
    "local_snp_file2": "",  # [your_args[2]_path] - TO BE FILLED
    # RT COVERAGE SETTINGS (STRAND-AWARE!)
    "rt_coverage_downstream": 100,  # nt downstream in RNA 5'→3' direction
    "include_probe_in_coverage": False,  # Strictly downstream
    # RTBC BARCODE SETTINGS
    "add_rtbc_barcode": True,
    "rtbc_sequence": "/5Phos/TGACTTGAGGAT",  # Configurable RTBC with phosphorylation
    # BLAST SPECIFICITY
    "blast_database_path": "",  # TO BE FILLED - path to BLAST database
    "require_unique_hits_only": True,  # Remove multi-hit probes
    "run_blast_analysis": True,  # Enable BLAST specificity check
    # OLIGOSTAN PARAMETERS (inherited from original)
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
