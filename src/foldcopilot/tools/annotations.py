"""Annotation tools — AlphaMissense pathogenicity + AlphaFill cofactor transplantation.

v0.6: sequence -> structure -> missense pathogenicity -> cofactor context in one prompt.
"""

from __future__ import annotations

import httpx

from foldcopilot.clients import alphafill_client, alphamissense_client


async def get_missense_landscape(uniprot_id: str) -> dict:
    """Get AlphaMissense pathogenicity landscape for a protein.

    Returns per-residue mean pathogenicity scores and variant-level
    classification (likely_pathogenic, likely_benign, ambiguous).
    """
    return await alphamissense_client.get_missense_predictions(uniprot_id)


async def get_cofactors(uniprot_id: str) -> dict:
    """Get AlphaFill transplanted cofactors, ligands, and metal ions.

    Returns compounds transplanted from experimental PDB structures into
    the AlphaFold model based on structural homology.
    """
    return await alphafill_client.get_alphafill_data(uniprot_id)


async def get_full_annotation(uniprot_id: str) -> dict:
    """Get comprehensive annotation: AlphaMissense + AlphaFill + confidence.

    Combines missense pathogenicity landscape with cofactor context and
    confidence assessment. This is the "one prompt" integration from the plan:
    sequence -> structure -> cofactor transplantation -> missense pathogenicity.
    """
    async with httpx.AsyncClient(timeout=60) as client:
        missense = await alphamissense_client.get_missense_predictions(
            uniprot_id, client=client
        )
        cofactors = await alphafill_client.get_alphafill_data(
            uniprot_id, client=client
        )

    # Cross-reference: flag pathogenic variants near cofactor binding sites
    hotspots = _find_cofactor_pathogenicity_hotspots(missense, cofactors)

    return {
        "uniprot_id": uniprot_id,
        "missense": missense,
        "cofactors": cofactors,
        "cofactor_pathogenicity_hotspots": hotspots,
        "caveats": [
            "AlphaMissense scores are predictions, not experimental measurements. "
            "Validate pathogenic variants with clinical data.",
            "AlphaFill transplantation is based on structural homology. "
            "Transplanted ligands may not reflect true binding.",
            "This annotation is for research use only.",
        ],
    }


def _find_cofactor_pathogenicity_hotspots(
    missense: dict, cofactors: dict
) -> list[dict]:
    """Find residues near cofactor sites with high pathogenicity scores.

    These are likely functionally important residues where mutations
    could disrupt cofactor binding.
    """
    if not missense.get("available") or not cofactors.get("available"):
        return []

    landscape = missense.get("residue_landscape", {})
    transplants = cofactors.get("transplants", [])

    if not transplants:
        return []

    # For now, flag all residues with high pathogenicity that are
    # in proteins with cofactor binding (positional cross-ref needs
    # 3D coordinates which we'd get from the structure)
    hotspots = []
    for pos_str, data in landscape.items():
        if data.get("classification") == "likely_pathogenic":
            hotspots.append({
                "residue_position": int(pos_str),
                "mean_pathogenicity": data["mean_pathogenicity"],
                "n_cofactors_in_protein": len(transplants),
                "note": "Likely pathogenic residue in a protein with "
                        f"{len(transplants)} transplanted cofactor(s). "
                        "Check 3D proximity to binding sites.",
            })

    # Limit to top 20 most pathogenic
    hotspots.sort(key=lambda h: h["mean_pathogenicity"], reverse=True)
    return hotspots[:20]
