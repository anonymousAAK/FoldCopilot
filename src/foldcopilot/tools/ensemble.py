"""Ensemble comparison and cross-model disagreement detection.

The second moat: when two models agree at a residue/interface, surface
high confidence. When they disagree, flag it. No published MCP does this.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from foldcopilot.models.ensemble import (
    AgreementLevel,
    DisagreementSpan,
    EnsembleReport,
    ModelSummary,
    ResidueDisagreement,
)


def parse_ca_coords(pdb_text: str) -> list[tuple[int, str, np.ndarray]]:
    """Extract CA atom coordinates from PDB text.

    Returns list of (residue_index, residue_name, xyz_array).
    """
    cas = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            res_seq = int(line[22:26].strip())
            res_name = line[17:20].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            cas.append((res_seq, res_name, np.array([x, y, z])))
    return cas


def parse_plddt_from_pdb(pdb_text: str) -> dict[int, float]:
    """Extract per-residue pLDDT from PDB B-factor column (CA atoms)."""
    plddt = {}
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            res_seq = int(line[22:26].strip())
            b_factor = float(line[60:66].strip())
            plddt[res_seq] = b_factor
    return plddt


def compute_ca_distances(
    coords_a: list[tuple[int, str, np.ndarray]],
    coords_b: list[tuple[int, str, np.ndarray]],
) -> dict[int, float]:
    """Compute per-residue CA-CA distance between two structures.

    Matches by residue index. Returns {residue_index: distance_angstroms}.
    """
    map_b = {res_idx: xyz for res_idx, _, xyz in coords_b}
    distances = {}
    for res_idx, _, xyz_a in coords_a:
        if res_idx in map_b:
            distances[res_idx] = float(np.linalg.norm(xyz_a - map_b[res_idx]))
    return distances


def compute_rmsd(distances: dict[int, float]) -> float:
    """Compute RMSD from per-residue CA distances."""
    if not distances:
        return float("inf")
    vals = np.array(list(distances.values()))
    return float(np.sqrt(np.mean(vals ** 2)))


def classify_agreement(
    plddt_a: float,
    plddt_b: float,
    ca_distance: float | None,
    plddt_threshold: float = 70.0,
    distance_threshold: float = 3.0,
) -> AgreementLevel:
    """Classify agreement between two models at a single residue."""
    a_confident = plddt_a >= plddt_threshold
    b_confident = plddt_b >= plddt_threshold

    if not a_confident and not b_confident:
        return AgreementLevel.BOTH_UNCERTAIN

    if ca_distance is None:
        # Can't compare structure, use pLDDT agreement only
        if a_confident and b_confident:
            return AgreementLevel.MODERATE_AGREE
        return AgreementLevel.BOTH_UNCERTAIN

    structures_match = ca_distance < distance_threshold

    if a_confident and b_confident:
        return AgreementLevel.STRONG_AGREE if structures_match else AgreementLevel.DISAGREE
    else:
        return AgreementLevel.MODERATE_AGREE if structures_match else AgreementLevel.DISAGREE


def find_disagreement_spans(
    residue_agreements: list[ResidueDisagreement],
    target_level: AgreementLevel = AgreementLevel.DISAGREE,
    min_length: int = 3,
) -> list[DisagreementSpan]:
    """Find contiguous spans of a given agreement level."""
    spans = []
    in_span = False
    span_residues: list[ResidueDisagreement] = []

    for r in residue_agreements:
        if r.agreement == target_level:
            if not in_span:
                in_span = True
                span_residues = []
            span_residues.append(r)
        else:
            if in_span and len(span_residues) >= min_length:
                spans.append(_make_span(span_residues, target_level))
            in_span = False
            span_residues = []

    # Handle span at end
    if in_span and len(span_residues) >= min_length:
        spans.append(_make_span(span_residues, target_level))

    return spans


def _make_span(
    residues: list[ResidueDisagreement], level: AgreementLevel
) -> DisagreementSpan:
    distances = [r.ca_distance for r in residues if r.ca_distance is not None]
    mean_dist = float(np.mean(distances)) if distances else None
    mean_a = float(np.mean([r.plddt_a for r in residues]))
    mean_b = float(np.mean([r.plddt_b for r in residues]))

    if level == AgreementLevel.DISAGREE:
        interpretation = (
            f"Models disagree across {len(residues)} residues "
            f"(mean CA distance: {mean_dist:.1f} A). "
            f"Both models are confident (pLDDT {mean_a:.0f} vs {mean_b:.0f}) "
            f"but predict different structures. Experimental validation recommended."
            if mean_dist
            else f"Models disagree across {len(residues)} residues. "
            f"Experimental validation recommended."
        )
    elif level == AgreementLevel.STRONG_AGREE:
        interpretation = (
            f"Models strongly agree across {len(residues)} residues "
            f"(mean CA distance: {mean_dist:.1f} A, pLDDT {mean_a:.0f}/{mean_b:.0f}). "
            f"High structural confidence."
            if mean_dist
            else f"Models agree across {len(residues)} residues."
        )
    else:
        interpretation = f"{level.value} across {len(residues)} residues."

    return DisagreementSpan(
        start=residues[0].residue_index,
        end=residues[-1].residue_index,
        length=len(residues),
        mean_ca_distance=round(mean_dist, 2) if mean_dist else None,
        mean_plddt_a=round(mean_a, 1),
        mean_plddt_b=round(mean_b, 1),
        agreement=level,
        interpretation=interpretation,
    )


def _generate_interpretation(
    agreement_frac: float,
    disagree_frac: float,
    rmsd: float | None,
    model_a: str,
    model_b: str,
) -> str:
    """Generate human-readable interpretation of ensemble comparison."""
    parts = []

    if agreement_frac > 0.9:
        parts.append(
            f"{model_a} and {model_b} show strong overall agreement "
            f"({agreement_frac:.0%} of residues)."
        )
    elif agreement_frac > 0.7:
        parts.append(
            f"{model_a} and {model_b} show moderate agreement "
            f"({agreement_frac:.0%} of residues)."
        )
    else:
        parts.append(
            f"{model_a} and {model_b} show significant disagreement "
            f"({disagree_frac:.0%} of residues disagree)."
        )

    if rmsd is not None:
        if rmsd < 2.0:
            parts.append(f"Global CA-RMSD is {rmsd:.1f} A (very similar structures).")
        elif rmsd < 5.0:
            parts.append(f"Global CA-RMSD is {rmsd:.1f} A (moderately similar).")
        else:
            parts.append(f"Global CA-RMSD is {rmsd:.1f} A (substantially different structures).")

    if disagree_frac > 0.1:
        parts.append(
            "Regions where both models are confident but predict different "
            "structures should be prioritized for experimental validation."
        )

    return " ".join(parts)


async def compare_structures(
    pdb_a: str,
    pdb_b: str,
    model_a_name: str = "model_a",
    model_b_name: str = "model_b",
    plddt_threshold: float = 70.0,
    distance_threshold: float = 3.0,
) -> dict:
    """Compare two predicted structures and detect disagreement.

    This is the ensemble cross-model disagreement feature.
    Feed it PDB outputs from two different backends (e.g., Boltz-2 and AF3)
    and get a detailed disagreement report.
    """
    # Parse structures
    coords_a = parse_ca_coords(pdb_a)
    coords_b = parse_ca_coords(pdb_b)
    plddt_a = parse_plddt_from_pdb(pdb_a)
    plddt_b = parse_plddt_from_pdb(pdb_b)

    if not coords_a or not coords_b:
        return {"error": "Could not parse CA atoms from one or both PDB inputs."}

    # Compute distances
    ca_distances = compute_ca_distances(coords_a, coords_b)
    rmsd = compute_rmsd(ca_distances) if ca_distances else None

    # Per-residue agreement
    all_residues = sorted(set(plddt_a.keys()) & set(plddt_b.keys()))
    if not all_residues:
        return {"error": "No overlapping residues between structures."}

    residue_agreements = []
    for res_idx in all_residues:
        pa = plddt_a[res_idx]
        pb = plddt_b[res_idx]
        dist = ca_distances.get(res_idx)
        agreement = classify_agreement(
            pa, pb, dist, plddt_threshold, distance_threshold
        )
        residue_agreements.append(ResidueDisagreement(
            residue_index=res_idx,
            model_a=model_a_name,
            model_b=model_b_name,
            plddt_a=round(pa, 1),
            plddt_b=round(pb, 1),
            ca_distance=round(dist, 2) if dist is not None else None,
            agreement=agreement,
        ))

    # Count agreement levels
    total = len(residue_agreements)
    counts = {}
    for level in AgreementLevel:
        counts[level] = sum(1 for r in residue_agreements if r.agreement == level) / total

    # Find spans
    disagree_spans = find_disagreement_spans(
        residue_agreements, AgreementLevel.DISAGREE, min_length=3
    )
    agree_spans = find_disagreement_spans(
        residue_agreements, AgreementLevel.STRONG_AGREE, min_length=5
    )

    # pLDDT correlation
    plddt_vals_a = np.array([plddt_a[r] for r in all_residues])
    plddt_vals_b = np.array([plddt_b[r] for r in all_residues])
    plddt_corr = float(np.corrcoef(plddt_vals_a, plddt_vals_b)[0, 1])

    # Model summaries
    summary_a = ModelSummary(
        model_name=model_a_name,
        mean_plddt=round(float(np.mean(plddt_vals_a)), 1),
        median_plddt=round(float(np.median(plddt_vals_a)), 1),
        residue_count=len(plddt_a),
    )
    summary_b = ModelSummary(
        model_name=model_b_name,
        mean_plddt=round(float(np.mean(plddt_vals_b)), 1),
        median_plddt=round(float(np.median(plddt_vals_b)), 1),
        residue_count=len(plddt_b),
    )

    agreement_frac = counts[AgreementLevel.STRONG_AGREE] + counts[AgreementLevel.MODERATE_AGREE]

    report = EnsembleReport(
        models=[summary_a, summary_b],
        residue_count=total,
        mean_ca_rmsd=round(rmsd, 2) if rmsd is not None else None,
        plddt_correlation=round(plddt_corr, 3),
        agreement_fraction=round(agreement_frac, 3),
        strong_agree_fraction=round(counts[AgreementLevel.STRONG_AGREE], 3),
        moderate_agree_fraction=round(counts[AgreementLevel.MODERATE_AGREE], 3),
        disagree_fraction=round(counts[AgreementLevel.DISAGREE], 3),
        both_uncertain_fraction=round(counts[AgreementLevel.BOTH_UNCERTAIN], 3),
        disagreement_spans=disagree_spans,
        high_confidence_consensus_spans=agree_spans,
        interpretation=_generate_interpretation(
            agreement_frac,
            counts[AgreementLevel.DISAGREE],
            rmsd,
            model_a_name,
            model_b_name,
        ),
    )
    report.add_standard_caveats()

    return report.model_dump()


async def compare_by_paths(
    pdb_path_a: str,
    pdb_path_b: str,
    model_a_name: str = "model_a",
    model_b_name: str = "model_b",
    plddt_threshold: float = 70.0,
    distance_threshold: float = 3.0,
) -> dict:
    """Compare two PDB files by path. Convenience wrapper for compare_structures."""
    path_a = Path(pdb_path_a)
    path_b = Path(pdb_path_b)

    if not path_a.exists():
        return {"error": f"File not found: {pdb_path_a}"}
    if not path_b.exists():
        return {"error": f"File not found: {pdb_path_b}"}

    pdb_a = path_a.read_text()
    pdb_b = path_b.read_text()

    return await compare_structures(
        pdb_a, pdb_b, model_a_name, model_b_name, plddt_threshold, distance_threshold
    )
