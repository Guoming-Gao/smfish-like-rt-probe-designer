# config.py - smfish-like-rt-probe-designer Configuration (Local Files + Stringent Filtering)
import os

# =============================================================================
# MAIN CONFIGURATION (LOCAL FILES + STRINGENT FILTERING)
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
    "transcript_selection": "longest",
    "output_directory": "/Users/gmgao/Dropbox/Caltech_PostDoc_GuttmanLab/constructs_and_smiFISH/smFISH_like_focusedRT-XCI",
    # LOCAL FILE PATHS (YOUR SERVER FILES)
    "local_gtf_path": "/Volumes/guttman-1/annotations/mm10/Mus_musculus.GRCm38.96.sorted.gtf",
    "local_genome_fasta_path": "/Volumes/guttman-1/genomes/mm10/GRCm38_68.fa",
    "snp_file_path": "/Volumes/guttman-1/data/snps/Bl6xCast.mm10.snps",
    # RT COVERAGE SETTINGS (UNCHANGED as requested)
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
    "salt_concentration": 0.115,
    # STRINGENT FILTERING (FIXED - all features restored)
    "pnas_filter_rules": [1, 2, 3, 4, 5],  # ALL 5 PNAS rules (FIXED)
    "use_dustmasker": True,  # ENABLED (FIXED)
    "max_masked_percent": 0.1,
    # SNP COVERAGE FILTERING (FIXED - missing parameter added)
    "min_snp_coverage_for_final": 2,  # FIXED: Minimum SNPs for final output
    "generate_all_probes_file": False,  # FIXED: Skip ALL.csv generation
    "focus_on_high_quality_only": True,  # FIXED: Only process filtered probes
    # FASTA OUTPUT SETTINGS
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
# TEST GENE SETS (unchanged as requested)
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
# DEMO GENE COORDINATES (mm10) - For when local files aren't accessible
# =============================================================================

DEMO_GENE_COORDS = {
    "Nanog": {
        "chromosome": "6",
        "start": 122707489,
        "end": 122714633,
        "strand": 1,
        "biotype": "protein_coding",
    },
    "Pou5f1": {
        "chromosome": "17",
        "start": 35506018,
        "end": 35510772,
        "strand": 1,
        "biotype": "protein_coding",
    },
    "Sox2": {
        "chromosome": "3",
        "start": 34650840,
        "end": 34652882,
        "strand": -1,
        "biotype": "protein_coding",
    },
}
