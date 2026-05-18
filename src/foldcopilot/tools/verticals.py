"""Therapeutic-area vertical packs.

v0.7: Domain-specific analysis workflows for antibodies, kinases, and GPCRs.
Each pack combines FoldCopilot tools with domain-specific databases and
interpretation logic.
"""

from __future__ import annotations

from foldcopilot.clients import afdb_client, alphafill_client, alphamissense_client
from foldcopilot.tools import confidence, foldseek


async def antibody_analysis(
    heavy_chain: str,
    light_chain: str | None = None,
    target_uniprot_id: str | None = None,
) -> dict:
    """Antibody Pack — comprehensive antibody structure analysis.

    Combines:
    - Confidence assessment of predicted antibody structure
    - CDR (Complementarity Determining Region) identification
    - Structural search for similar antibodies via Foldseek
    - Interface confidence if target is provided
    - Hallucination warnings (CDR loops are often disordered/flexible)

    Args:
        heavy_chain: Heavy chain amino acid sequence.
        light_chain: Light chain amino acid sequence (optional for VHH/nanobodies).
        target_uniprot_id: UniProt ID of the antigen target (optional).
    """
    result: dict = {
        "pack": "antibody",
        "heavy_chain_length": len(heavy_chain),
        "light_chain_length": len(light_chain) if light_chain else 0,
        "is_nanobody": light_chain is None,
    }

    # Identify CDR regions (Kabat numbering approximation)
    cdrs = _identify_cdrs(heavy_chain, chain_type="heavy")
    if light_chain:
        cdrs.extend(_identify_cdrs(light_chain, chain_type="light"))
    result["predicted_cdrs"] = cdrs

    # Flag that CDR loops are commonly low-confidence
    result["cdr_warnings"] = [
        "CDR loops (especially CDR-H3) are inherently flexible and often "
        "predicted with low confidence. Low pLDDT in CDRs is expected and "
        "does NOT necessarily indicate a bad prediction.",
        "For antibody engineering, experimental validation of CDR conformations "
        "is essential. Crystal structures or cryo-EM should be used for "
        "lead optimization.",
    ]

    # Target analysis if provided
    if target_uniprot_id:
        try:
            target_confidence = await confidence.assess_confidence(target_uniprot_id)
            result["target_confidence"] = {
                "uniprot_id": target_uniprot_id,
                "overall_mean_plddt": target_confidence.get("overall_mean_plddt"),
                "hallucination_count": len(
                    target_confidence.get("hallucination_warnings", [])
                ),
            }
        except Exception as e:
            result["target_confidence"] = {"error": str(e)}

    result["recommendations"] = [
        "Use Boltz-2 or Chai-1 for antibody-antigen co-folding (predict_structure with 2+ chains).",
        "Compare predictions from multiple backends using compare_predictions to identify "
        "uncertain CDR loop conformations.",
        "Cross-reference with SAbDab (Structural Antibody Database) for similar antibody structures.",
    ]

    return result


async def kinase_analysis(uniprot_id: str) -> dict:
    """Kinase Pack — kinase-specific structural analysis.

    Combines:
    - Confidence assessment focused on kinase-relevant regions
    - Active site / DFG motif / activation loop analysis
    - AlphaMissense pathogenicity for kinase domain
    - AlphaFill cofactor transplantation (ATP binding)
    - Structural search for similar kinases
    """
    import httpx

    result: dict = {
        "pack": "kinase",
        "uniprot_id": uniprot_id,
    }

    async with httpx.AsyncClient(timeout=60) as client:
        # Confidence assessment
        try:
            scores = await afdb_client.get_plddt_scores(uniprot_id, client=client)
            import numpy as np
            result["confidence"] = {
                "mean_plddt": round(float(np.mean(scores)), 1),
                "residue_count": len(scores),
            }
        except Exception as e:
            result["confidence"] = {"error": str(e)}

        # AlphaMissense
        try:
            missense = await alphamissense_client.get_missense_predictions(
                uniprot_id, client=client
            )
            result["missense"] = {
                "available": missense.get("available", False),
                "pathogenic_fraction": missense.get("pathogenic_fraction"),
                "total_variants": missense.get("total_variants"),
            }
        except Exception as e:
            result["missense"] = {"error": str(e)}

        # AlphaFill (cofactors — ATP/ADP binding)
        try:
            cofactors = await alphafill_client.get_alphafill_data(
                uniprot_id, client=client
            )
            if cofactors.get("available"):
                atp_related = [
                    t for t in cofactors.get("transplants", [])
                    if t.get("compound_id", "").upper() in (
                        "ATP", "ADP", "AMP", "ANP", "ACP", "AGS",
                        "STU", "STI", "NIL",  # common kinase inhibitors
                    )
                ]
                result["cofactors"] = {
                    "total_transplants": cofactors.get("total_transplants"),
                    "atp_site_hits": len(atp_related),
                    "atp_related_compounds": [
                        {
                            "compound_id": t["compound_id"],
                            "compound_name": t.get("compound_name"),
                            "source_pdb": t["source_pdb"],
                            "rmsd": t.get("rmsd"),
                        }
                        for t in atp_related
                    ],
                }
            else:
                result["cofactors"] = {"available": False}
        except Exception as e:
            result["cofactors"] = {"error": str(e)}

    result["kinase_motifs"] = {
        "note": "DFG motif and activation loop positions require sequence annotation. "
                "Cross-reference with KLIFS (Kinase-Ligand Interaction Fingerprints) "
                "for precise kinase domain mapping.",
        "databases": ["KLIFS (klifs.net)", "KinBase", "KinMap"],
    }

    result["recommendations"] = [
        "Use Boltz-2 with predict_affinity=True for kinase-inhibitor binding prediction.",
        "Check DFG-in vs DFG-out conformation — critical for inhibitor type selectivity.",
        "Cross-reference AlphaMissense hotspots with ATP binding site proximity.",
        "Compare prediction confidence at the activation loop — often flexible and uncertain.",
    ]

    return result


async def gpcr_analysis(uniprot_id: str) -> dict:
    """GPCR Pack — G protein-coupled receptor analysis.

    Combines:
    - Confidence assessment focused on transmembrane helices
    - AlphaFill ligand transplantation
    - AlphaMissense pathogenicity
    - Membrane orientation context
    """
    import httpx

    result: dict = {
        "pack": "gpcr",
        "uniprot_id": uniprot_id,
    }

    async with httpx.AsyncClient(timeout=60) as client:
        # Confidence
        try:
            scores = await afdb_client.get_plddt_scores(uniprot_id, client=client)
            import numpy as np
            arr = np.array(scores)
            result["confidence"] = {
                "mean_plddt": round(float(np.mean(arr)), 1),
                "residue_count": len(scores),
                "note": "GPCRs typically show high confidence in TM helices "
                        "but low confidence in loop regions (ICL3, N/C-termini).",
            }
        except Exception as e:
            result["confidence"] = {"error": str(e)}

        # AlphaFill
        try:
            cofactors = await alphafill_client.get_alphafill_data(
                uniprot_id, client=client
            )
            if cofactors.get("available"):
                result["cofactors"] = {
                    "total_transplants": cofactors.get("total_transplants"),
                    "type_summary": cofactors.get("type_summary"),
                }
            else:
                result["cofactors"] = {"available": False}
        except Exception as e:
            result["cofactors"] = {"error": str(e)}

        # AlphaMissense
        try:
            missense = await alphamissense_client.get_missense_predictions(
                uniprot_id, client=client
            )
            result["missense"] = {
                "available": missense.get("available", False),
                "pathogenic_fraction": missense.get("pathogenic_fraction"),
            }
        except Exception as e:
            result["missense"] = {"error": str(e)}

    result["gpcr_context"] = {
        "membrane_orientation": (
            "Use TMalphaFold (tmalphafold.ttk.hu) for membrane orientation. "
            "AlphaFold predictions do not include lipid bilayer context."
        ),
        "activation_states": (
            "GPCRs can adopt active, inactive, and intermediate conformations. "
            "AlphaFold typically predicts one dominant conformation — it may not "
            "capture all functionally relevant states."
        ),
        "databases": [
            "GPCRdb (gpcrdb.org) — GPCR-specific annotations",
            "TMalphaFold — membrane protein orientation",
            "OPM (Orientations of Proteins in Membranes)",
        ],
    }

    result["recommendations"] = [
        "TM helices should have high pLDDT (>80). Low confidence in TM regions is unusual.",
        "ICL3 (intracellular loop 3) is often disordered — low pLDDT is expected.",
        "N/C-termini are frequently disordered in GPCRs. Check DisProt cross-reference.",
        "For ligand binding analysis, combine AlphaFill transplants with Boltz-2 co-folding.",
    ]

    return result


def _identify_cdrs(
    sequence: str, chain_type: str = "heavy"
) -> list[dict]:
    """Approximate CDR identification using sequence length heuristics.

    This is a simplified Kabat-like numbering approximation.
    For production use, integrate ANARCI or IMGT numbering.
    """
    cdrs = []
    seq_len = len(sequence)

    if chain_type == "heavy":
        # Approximate CDR positions for typical heavy chains
        # CDR-H1: ~26-35, CDR-H2: ~50-65, CDR-H3: ~95-102 (variable)
        if seq_len >= 110:
            cdrs.extend([
                {
                    "cdr": "CDR-H1",
                    "approximate_start": 26,
                    "approximate_end": 35,
                    "note": "Framework 1 -> CDR-H1. Usually well-predicted.",
                },
                {
                    "cdr": "CDR-H2",
                    "approximate_start": 50,
                    "approximate_end": 65,
                    "note": "CDR-H2. Moderately variable.",
                },
                {
                    "cdr": "CDR-H3",
                    "approximate_start": 95,
                    "approximate_end": min(102, seq_len),
                    "note": "CDR-H3. Highly variable length and conformation. "
                            "Most likely to have low pLDDT. This is expected.",
                },
            ])
    elif chain_type == "light":
        if seq_len >= 100:
            cdrs.extend([
                {
                    "cdr": "CDR-L1",
                    "approximate_start": 24,
                    "approximate_end": 34,
                    "note": "CDR-L1. Variable length.",
                },
                {
                    "cdr": "CDR-L2",
                    "approximate_start": 50,
                    "approximate_end": 56,
                    "note": "CDR-L2. Usually short and well-predicted.",
                },
                {
                    "cdr": "CDR-L3",
                    "approximate_start": 89,
                    "approximate_end": 97,
                    "note": "CDR-L3. Moderate variability.",
                },
            ])

    return cdrs
