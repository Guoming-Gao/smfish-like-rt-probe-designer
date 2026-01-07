# Usage Guide: smFISH-like RT Probe Designer

This pipeline is designed to generate highly specific probes for smFISH-like applications, including forward primer design for RT-PCR validation.

## 4-Phase Workflow

### Phase 1: Candidate Design
**Script**: `design_candidate_probes.py`
- Generates raw probe candidates based on melting temperature (Tm), GC content, and sequence complexity.
- Filters for homopolymers and other unfavorable motifs.

### Phase 2: Specificity Validation (BLAST)
**Script**: `validate_probe_specificity.py`
- Validates each candidate probe against the local genome using BLAST.
- **Parameters**: `word_size=11`, `evalue=0.1`.
- Probes are marked as "Unique" if they have exactly one genomic hit with significant homology.

### Phase 3: Forward Primer Design
**Script**: `design_forward_primers.py`
- Designs forward primers 200-750bp upstream of the probe.
- **Complexity Filter**: Pre-screens for repetitive motifs.
- **BLAST Specificity**: Ensures primers are uniquely targetable.
- *Note*: If primers are not found in the 750bp window (common in repetitive regions like Xist), the `--max-distance` can be extended.

### Phase 4: Synthesis Preparation
**Notebook Cell**: Synthesis Prep
- Consolidates results into a single synthesis-ready format.

## Troubleshooting

### No Primers Found
If a probe returns no forward primers, it is often due to biological constraints (repetitive DNA) in the 750bp window. See `docs/Xist_Primer_Design_Constraints_Note.md` for a detailed case study on this phenomenon.
