"""Design validation tool for externally designed proteins.

Accepts PDB content from design tools (BindCraft, RFdiffusion3, ProteinMPNN,
LigandMPNN) and runs confidence analysis to validate the design quality.
"""

from __future__ import annotations


def validate_design(
    pdb_content: str,
    design_tool: str = "unknown",
    design_intent: str = "",
) -> dict:
    """Validate a protein design by analyzing its structural confidence.

    Parses the PDB, extracts per-residue B-factors (pLDDT if from AF2/AF3),
    identifies low-confidence regions, and provides design-specific guidance.

    Args:
        pdb_content: PDB format string of the designed structure.
        design_tool: Name of the design tool used (e.g., "bindcraft", "rfdiffusion3").
        design_intent: What the design is for (e.g., "binder for EGFR", "enzyme scaffold").

    Returns:
        Dict with design quality assessment.
    """
    lines = pdb_content.strip().splitlines()

    # Extract CA atoms and B-factors
    residues = {}
    for line in lines:
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            res_seq = int(line[22:26].strip())
            chain = line[21]
            b_factor = float(line[60:66].strip())
            residues[(chain, res_seq)] = b_factor

    if not residues:
        return {
            "status": "error",
            "message": "No CA atoms found in PDB content. Ensure valid PDB format.",
        }

    b_factors = list(residues.values())
    n_residues = len(b_factors)
    mean_confidence = sum(b_factors) / n_residues

    # Confidence bucketing (pLDDT-style)
    very_high = sum(1 for b in b_factors if b > 90) / n_residues * 100
    confident = sum(1 for b in b_factors if 70 < b <= 90) / n_residues * 100
    low = sum(1 for b in b_factors if 50 < b <= 70) / n_residues * 100
    very_low = sum(1 for b in b_factors if b <= 50) / n_residues * 100

    # Find low-confidence spans
    low_conf_regions = []
    current_span = []
    sorted_residues = sorted(residues.items())
    for (chain, res_seq), b in sorted_residues:
        if b < 70:
            current_span.append((chain, res_seq, b))
        else:
            if len(current_span) >= 3:
                low_conf_regions.append({
                    "chain": current_span[0][0],
                    "start": current_span[0][1],
                    "end": current_span[-1][1],
                    "length": len(current_span),
                    "mean_confidence": round(sum(x[2] for x in current_span) / len(current_span), 1),
                })
            current_span = []
    if len(current_span) >= 3:
        low_conf_regions.append({
            "chain": current_span[0][0],
            "start": current_span[0][1],
            "end": current_span[-1][1],
            "length": len(current_span),
            "mean_confidence": round(sum(x[2] for x in current_span) / len(current_span), 1),
        })

    # Design-tool-specific guidance
    guidance = _get_design_guidance(design_tool, mean_confidence, very_low, low_conf_regions)

    # Quality verdict
    if mean_confidence > 85 and very_low < 5:
        verdict = "high_quality"
        verdict_text = "Design has high structural confidence. Suitable for experimental validation."
    elif mean_confidence > 70 and very_low < 15:
        verdict = "moderate_quality"
        verdict_text = "Design has moderate confidence. Consider redesigning low-confidence regions before experimental work."
    else:
        verdict = "low_quality"
        verdict_text = "Design has significant low-confidence regions. Structural reliability is uncertain. Recommend redesign."

    return {
        "status": "success",
        "design_tool": design_tool,
        "design_intent": design_intent,
        "n_residues": n_residues,
        "mean_confidence": round(mean_confidence, 1),
        "confidence_distribution": {
            "very_high_pct": round(very_high, 1),
            "confident_pct": round(confident, 1),
            "low_pct": round(low, 1),
            "very_low_pct": round(very_low, 1),
        },
        "low_confidence_regions": low_conf_regions,
        "verdict": verdict,
        "verdict_text": verdict_text,
        "guidance": guidance,
        "recommendation": "Re-predict with Boltz-2 or Protenix-v2 and compare via compare_predictions() for ensemble validation." if verdict != "high_quality" else "Consider direct experimental validation.",
    }


def _get_design_guidance(tool: str, mean_conf: float, very_low_pct: float, regions: list) -> list[str]:
    """Tool-specific design guidance."""
    guidance = []

    tool_lower = tool.lower()

    if "bindcraft" in tool_lower:
        guidance.append("BindCraft designs: interface residues with pLDDT < 70 may indicate weak binding. Check interface PAE with compare_predictions().")
        if very_low_pct > 10:
            guidance.append("Consider increasing BindCraft optimization cycles or adjusting hotspot residues.")

    elif "rfdiffusion" in tool_lower:
        guidance.append("RFdiffusion3 designs: low-confidence backbone regions may need sequence optimization via ProteinMPNN/LigandMPNN.")
        if regions:
            guidance.append(f"Found {len(regions)} low-confidence spans — consider constraining these regions in the diffusion process.")

    elif "proteinmpnn" in tool_lower or "ligandmpnn" in tool_lower:
        guidance.append("Sequence design output: re-predict structure with Boltz-2 to verify the designed sequence folds as intended.")
        guidance.append("Compare designed vs target structure using compare_predictions() for self-consistency check.")

    else:
        guidance.append("Re-predict this design with at least 2 backends (e.g., Boltz-2 + Protenix-v2) and check ensemble agreement.")

    if mean_conf < 70:
        guidance.append("Overall confidence is low. The design may not fold as predicted. Consider iterating on the design before wet-lab work.")

    return guidance
