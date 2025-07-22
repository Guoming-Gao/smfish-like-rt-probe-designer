# utils/oligostan_core.py - Modified Oligostan Core (No FLAP sequences)

import numpy as np
from .thermodynamics import dg_calc_rna_37, dg37_score_calc
from .filters import is_ok_4_pnas_filter, is_ok_4_gc_filter
from .sequence_utils import determine_region_type


def which_max_r(x):
    """Exact translation of R's WhichMax function"""
    max_val = np.max(x)
    max_indices = np.where(x == max_val)[0]

    if len(max_indices) >= 2:
        return [0, max_val]  # Return 0 for ties (R behavior)
    else:
        return [max_indices[0] + 1, max_val]  # Return 1-based index


def get_probes_from_sequence(
    seq,
    min_size_probe=26,
    max_size_probe=32,
    desired_dg=-32,
    min_score_value=0.9,
    inc_betw_prob=2,
):
    """Modified probe generation from Oligostan (same core algorithm)"""
    if isinstance(seq, list):
        seq = "".join(seq).upper()

    diff_size = max_size_probe - min_size_probe

    # Build thermodynamic matrix for all probe sizes
    the_tms_tmp = dg_calc_rna_37(seq, probe_length=max_size_probe)
    nb_of_probes = len(the_tms_tmp)

    if diff_size > 0:
        all_columns = [the_tms_tmp]  # Start with max size

        for i in range(diff_size - 1, -1, -1):
            probe_length = min_size_probe + i
            dg_values = dg_calc_rna_37(seq, probe_length=probe_length)
            min_len = min(len(dg_values), nb_of_probes)
            all_columns.insert(0, dg_values[:min_len])

        # Build matrix: rows = positions, columns = probe sizes
        max_len = max(len(col) for col in all_columns)
        the_tms_matrix = np.full((max_len, diff_size + 1), np.nan)

        for col_idx, column in enumerate(all_columns):
            the_tms_matrix[: len(column), col_idx] = column
    else:
        the_tms_matrix = np.array(the_tms_tmp).reshape(-1, 1)

    # Calculate scores
    tm_scores = np.zeros_like(the_tms_matrix)
    for col in range(the_tms_matrix.shape[1]):
        valid_mask = ~np.isnan(the_tms_matrix[:, col])
        if np.any(valid_mask):
            tm_scores[valid_mask, col] = dg37_score_calc(
                the_tms_matrix[valid_mask, col], desired_dg
            )

    # Find best probe size for each position
    best_scores = []
    for row in range(tm_scores.shape[0]):
        row_data = tm_scores[row, :]
        if not np.all(np.isnan(row_data)) and not np.all(row_data == 0):
            valid_data = row_data[~np.isnan(row_data)]
            if len(valid_data) > 0:
                result = which_max_r(valid_data)
                best_scores.append(result)
            else:
                best_scores.append([0, 0])
        else:
            best_scores.append([0, 0])

    best_scores = np.array(best_scores)

    # Convert to probe sizes
    for i in range(len(best_scores)):
        if best_scores[i, 0] != 0:  # Not a tie
            best_scores[i, 0] = best_scores[i, 0] + (min_size_probe - 1)
        else:  # Tie case - use minimum size
            best_scores[i, 0] = min_size_probe

    # Add positions
    positions = np.arange(1, len(best_scores) + 1).reshape(-1, 1)
    best_scores_with_pos = np.column_stack([best_scores, positions])

    # Filter by score threshold
    validated_scores = best_scores_with_pos[
        best_scores_with_pos[:, 1] >= min_score_value
    ]

    if len(validated_scores) == 0:
        return None

    # Apply spacing constraint
    return apply_spacing_constraint(validated_scores, seq, inc_betw_prob)


def apply_spacing_constraint(validated_scores, seq, min_spacing):
    """Apply minimum spacing between probes"""
    if len(validated_scores.shape) == 1:
        validated_scores = validated_scores.reshape(1, -1)

    validated_scores = validated_scores[np.argsort(validated_scores[:, 2])]

    the_probes = []
    pointeur = 0

    while pointeur < len(seq):
        valid_tmp = validated_scores[validated_scores[:, 2] >= pointeur]

        if len(valid_tmp) > 0:
            if len(valid_tmp.shape) == 1:
                valid_tmp = valid_tmp.reshape(1, -1)

            probe_size = int(valid_tmp[0, 0])
            score = valid_tmp[0, 1]
            position = int(valid_tmp[0, 2])

            # Extract probe sequence
            probe_seq = seq[position - 1 : position - 1 + probe_size].upper()

            the_probes.append(
                {
                    "size": probe_size,
                    "score": score,
                    "position": position,
                    "sequence": probe_seq,
                    "genomic_start": position,  # Will be adjusted later
                    "genomic_end": position + probe_size - 1,
                }
            )

            pointeur = position + probe_size + min_spacing
        else:
            break

    return the_probes


def design_fish_probes(gene_data, config):
    """
    Main probe design function for FISH-RT probes
    Returns list of probe dictionaries with all required information
    """
    try:
        # Extract configuration
        min_size = config["probe_length_min"]
        max_size = config["probe_length_max"]
        min_score = config["score_min"]
        spacing = config["probe_spacing_min"]
        desired_dg = config["fixed_dg37_value"]

        # Design probes using Oligostan algorithm
        raw_probes = get_probes_from_sequence(
            gene_data["sequence"],
            min_size_probe=min_size,
            max_size_probe=max_size,
            desired_dg=desired_dg,
            min_score_value=min_score,
            inc_betw_prob=spacing,
        )

        if not raw_probes:
            return []

        # Process probes and add gene context
        processed_probes = []

        for probe in raw_probes:
            # Calculate genomic coordinates
            if gene_data["strand"] == 1:  # Forward strand
                genomic_start = gene_data["genomic_start"] + probe["position"] - 1
                genomic_end = genomic_start + probe["size"] - 1
                target_strand = "+"
                probe_strand = "-"  # Probe is complementary
            else:  # Reverse strand
                genomic_end = gene_data["genomic_end"] - probe["position"] + 1
                genomic_start = genomic_end - probe["size"] + 1
                target_strand = "-"
                probe_strand = "+"  # Probe is complementary

            # Determine if probe is in exon or intron
            region_type = determine_region_type(
                probe["position"],
                probe["position"] + probe["size"] - 1,
                gene_data["exon_regions"],
                gene_data["intron_regions"],
            )

            # Calculate quality metrics
            gc_count = probe["sequence"].count("G") + probe["sequence"].count("C")
            gc_percentage = gc_count / len(probe["sequence"])

            # Apply filters
            gc_filter = is_ok_4_gc_filter(
                probe["sequence"], config["gc_content_min"], config["gc_content_max"]
            )
            pnas_filter = is_ok_4_pnas_filter(
                probe["sequence"], config["pnas_filter_rules"]
            )

            # Create probe record
            probe_record = {
                # Gene information
                "GeneName": gene_data["gene_name"],
                "Species": config["species"],
                "Chromosome": gene_data["chromosome"],
                "RegionType": region_type,  # 'exon', 'intron', or 'exon-intron'
                # Probe sequence and position
                "Seq": probe["sequence"],
                "ProbeSize": probe["size"],
                "theStartPos": genomic_start,
                "theEndPos": genomic_end,
                # Strand information (CRITICAL for RT coverage)
                "Target_Strand": target_strand,  # Gene's RNA strand
                "Probe_Strand": probe_strand,  # Probe's strand (opposite)
                # Thermodynamic properties
                "dGOpt": desired_dg,
                "dGScore": probe["score"],
                "dG37": dg_calc_rna_37(
                    probe["sequence"], probe_length=len(probe["sequence"])
                )[0],
                # Quality metrics
                "GCpc": gc_percentage,
                "GCFilter": 1 if gc_filter else 0,
                "PNASFilter": 1 if pnas_filter else 0,
                # Placeholders for later analysis
                "SNPs_Covered_Count": 0,
                "SNPs_Covered_Positions": "",
                "RT_Coverage_Start": 0,
                "RT_Coverage_End": 0,
                "RTBC_5Prime_Sequence": "",  # Will be added later
                "NumberOfHits": 0,  # BLAST results
                "PercentAlignment": 0,  # BLAST results
                "UniqueHitName": "",  # BLAST results
            }

            processed_probes.append(probe_record)

        return processed_probes

    except Exception as e:
        raise Exception(
            f"Error designing probes for {gene_data.get('gene_name', 'unknown')}: {str(e)}"
        )
