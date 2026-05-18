"""Education mode — plain-language explanations for wet-lab biologists.

v0.8: --explain flag equivalent. Emits human-readable interpretation of
pLDDT, PAE, ipTM, and confidence reports with citations.
Plays to LLM strengths and to wet-lab biologist needs.
"""

from __future__ import annotations

from foldcopilot.models.confidence import ConfidenceBucket


def explain_plddt(score: float) -> dict:
    """Generate plain-language explanation of a pLDDT score.

    Targeted at wet-lab biologists who need actionable interpretation,
    not just numbers.
    """
    if score > 90:
        bucket = ConfidenceBucket.VERY_HIGH
        summary = "Very high confidence"
        meaning = (
            "AlphaFold is highly confident in both the backbone AND side-chain "
            "positions at this residue. The local structure is likely correct."
        )
        action = (
            "This region is suitable for structure-based analysis: active site "
            "identification, docking, mutation design, and molecular dynamics."
        )
        analogy = (
            "Think of this like a high-resolution crystal structure region "
            "(~1.5A resolution equivalent). You can trust atom positions."
        )
    elif score > 70:
        bucket = ConfidenceBucket.HIGH
        summary = "High confidence"
        meaning = (
            "AlphaFold is confident in the backbone trace (the overall fold) "
            "but less certain about exact side-chain rotamers."
        )
        action = (
            "The overall fold is trustworthy. Use for domain identification, "
            "fold classification, and general structural analysis. Be cautious "
            "about precise side-chain contacts."
        )
        analogy = (
            "Like a medium-resolution crystal structure (~2-3A). "
            "You know the shape but not every atomic detail."
        )
    elif score > 50:
        bucket = ConfidenceBucket.LOW
        summary = "Low confidence"
        meaning = (
            "AlphaFold is uncertain about this region. The predicted structure "
            "may be partially correct but could also be significantly wrong."
        )
        action = (
            "DO NOT use this region for precise structural analysis. "
            "It might be a flexible loop, a crystal packing artifact, "
            "or an intrinsically disordered region. Validate experimentally."
        )
        analogy = (
            "Like trying to photograph something that's moving. "
            "You get a blur — the general area is right but details are lost."
        )
    else:
        bucket = ConfidenceBucket.VERY_LOW
        summary = "Very low confidence — do not interpret"
        meaning = (
            "AlphaFold essentially does not know what this region looks like. "
            "The predicted coordinates are likely meaningless."
        )
        action = (
            "IGNORE this region structurally. It is almost certainly "
            "intrinsically disordered (no fixed 3D structure), a flexible linker, "
            "or a region that only folds upon binding a partner."
        )
        analogy = (
            "Like a cooked spaghetti noodle — it doesn't have a fixed shape. "
            "AlphaFold had to put atoms somewhere, but those positions are fiction."
        )

    return {
        "score": score,
        "bucket": bucket.value,
        "summary": summary,
        "what_it_means": meaning,
        "what_to_do": action,
        "analogy": analogy,
        "citation": (
            "Jumper et al., Nature 2021. pLDDT is a per-residue confidence "
            "metric (0-100) representing the predicted lDDT-Ca score. "
            "See also: EBI's 'Interpreting AlphaFold structures' guide."
        ),
    }


def explain_pae(mean_pae: float, context: str = "general") -> dict:
    """Generate plain-language explanation of PAE (Predicted Aligned Error).

    Args:
        mean_pae: Mean PAE value in Angstroms.
        context: "general", "interface", or "domain" for context-specific advice.
    """
    if mean_pae < 5:
        quality = "excellent"
        meaning = (
            "The relative positions of residues are predicted with high accuracy. "
            "Domain arrangements and interfaces are likely correct."
        )
        action = "Trust the relative domain positioning and interface geometry."
    elif mean_pae < 10:
        quality = "moderate"
        meaning = (
            "Moderate confidence in relative positioning. Individual domains "
            "are probably correct, but their exact arrangement may be off."
        )
        action = (
            "Trust individual domain structures but be cautious about "
            "inter-domain angles and distances."
        )
    elif mean_pae < 20:
        quality = "poor"
        meaning = (
            "Low confidence in relative positioning. Domains may be "
            "correctly folded individually but their arrangement is uncertain."
        )
        action = (
            "Analyze domains separately. Do NOT trust the relative "
            "orientation of domains or chains. The 'hinge' between domains "
            "is likely flexible in reality."
        )
    else:
        quality = "very poor"
        meaning = (
            "AlphaFold has essentially no idea how these regions relate "
            "to each other spatially. This often indicates intrinsic disorder "
            "or that the regions only adopt fixed positions when bound to partners."
        )
        action = (
            "These regions should be treated as independent structural units. "
            "Their relative positioning in the prediction is meaningless."
        )

    interface_note = ""
    if context == "interface":
        interface_note = (
            "For protein-protein interfaces, PAE > 10A strongly suggests "
            "the predicted interface is unreliable. The proteins may interact "
            "differently or not at all. Validate with co-IP, crosslinking-MS, "
            "or cryo-EM."
        )
    elif context == "domain":
        interface_note = (
            "For multi-domain proteins, high inter-domain PAE is common and "
            "often reflects real flexibility. The domains may 'wobble' relative "
            "to each other. Consider using SAXS or cryo-EM to determine the "
            "ensemble of domain arrangements."
        )

    return {
        "mean_pae_angstroms": mean_pae,
        "quality": quality,
        "what_it_means": meaning,
        "what_to_do": action,
        "context_note": interface_note or None,
        "citation": (
            "PAE (Predicted Aligned Error) represents the expected position "
            "error (in Angstroms) of residue X when the structure is aligned "
            "on residue Y. Evans et al., bioRxiv 2021."
        ),
    }


def explain_hallucination_warning(
    af_plddt: float,
    idr_source: str,
    severity: str,
) -> dict:
    """Explain a hallucination warning in plain language."""
    if severity == "high":
        headline = "AlphaFold is CONFIDENTLY WRONG here"
        explanation = (
            f"AlphaFold predicts an ordered structure (pLDDT {af_plddt:.0f}) "
            f"for this region, but {idr_source} — a curated database of "
            f"experimentally verified disorder — says this region is "
            f"intrinsically disordered (no fixed structure)."
        )
        impact = (
            "This is a hallucination. AlphaFold generated a plausible-looking "
            "structure that does not exist in reality. If you used this region "
            "for docking, active site analysis, or mutation design, your "
            "conclusions may be wrong."
        )
        what_to_do = (
            "1. Do NOT trust the atomic coordinates in this region.\n"
            "2. This region is likely a flexible loop, linker, or IDR.\n"
            "3. It may fold upon binding a partner (coupled folding).\n"
            "4. Use circular dichroism, NMR, or SAXS to confirm disorder.\n"
            "5. Check the literature for known binding-induced folding."
        )
    else:
        headline = "AlphaFold may be partially wrong here"
        explanation = (
            f"AlphaFold predicts moderate structure (pLDDT {af_plddt:.0f}) "
            f"for this region, but {idr_source} indicates disorder. "
            f"The prediction may capture a transient or partially ordered state."
        )
        impact = (
            "This region is ambiguous. It might be partially structured or "
            "conditionally ordered (e.g., folds upon binding). Treat with caution."
        )
        what_to_do = (
            "1. Consider this region uncertain, not definitively structured.\n"
            "2. Look for binding partners that might stabilize it.\n"
            "3. Experimental validation recommended before any design work."
        )

    return {
        "headline": headline,
        "severity": severity,
        "explanation": explanation,
        "real_world_impact": impact,
        "what_to_do": what_to_do,
        "background": (
            "AlphaFold 3 hallucinations in disordered regions: ~22% of IDR "
            "residues are falsely predicted as ordered (arXiv 2510.15939). "
            "This is a known limitation, not a bug."
        ),
    }


def explain_confidence_report(report: dict) -> dict:
    """Generate a plain-language summary of a full ConfidenceReport.

    Translates the structured report into actionable prose for a
    biologist who doesn't want to parse JSON.
    """
    plddt = report.get("overall_mean_plddt", 0)
    hallucinations = report.get("hallucination_warnings", [])
    idr_flags = report.get("idr_flags", [])
    chain_summaries = report.get("chain_summaries", [])
    pae = report.get("pae_summary")

    # Overall verdict
    if plddt > 85 and not hallucinations:
        verdict = "HIGH CONFIDENCE — this prediction is likely reliable"
        verdict_detail = (
            "The overall structure has high confidence with no hallucination "
            "warnings. Individual low-confidence regions (if any) are noted below."
        )
    elif plddt > 70 and len(hallucinations) <= 2:
        verdict = "MODERATE CONFIDENCE — mostly reliable with some caveats"
        verdict_detail = (
            "The core structure is probably correct but some regions need "
            "attention. Check the warnings below."
        )
    elif plddt > 50:
        verdict = "LOW CONFIDENCE — interpret with caution"
        verdict_detail = (
            "Significant portions of this structure are uncertain. "
            "Only use high-confidence regions for analysis."
        )
    else:
        verdict = "VERY LOW CONFIDENCE — this prediction is mostly unreliable"
        verdict_detail = (
            "Most of this protein lacks confident structure prediction. "
            "It may be largely disordered or require binding partners to fold."
        )

    # Low-confidence regions summary
    low_regions = []
    for chain in chain_summaries:
        for span in chain.get("low_confidence_spans", []):
            low_regions.append(
                f"Residues {span['start']}-{span['end']} "
                f"(pLDDT {span['mean_plddt']:.0f}, {span['length']} residues)"
            )

    # Hallucination summary
    hallucination_summary = []
    for h in hallucinations:
        hallucination_summary.append(
            f"Residues {h['start']}-{h['end']}: AF predicts order "
            f"(pLDDT {h['af_mean_plddt']:.0f}) but {h['idr_source']} "
            f"says DISORDERED [{h['severity']} severity]"
        )

    return {
        "verdict": verdict,
        "verdict_detail": verdict_detail,
        "overall_plddt": plddt,
        "plddt_explanation": explain_plddt(plddt),
        "low_confidence_regions": low_regions or ["None — all regions are confident"],
        "hallucination_warnings": hallucination_summary or ["None detected"],
        "idr_database_hits": len(idr_flags),
        "pae_explanation": explain_pae(
            pae["mean_pae"], context="general"
        ) if pae and pae.get("mean_pae") else None,
        "bottom_line": (
            f"{'Trust this structure for analysis.' if plddt > 80 and not hallucinations else ''}"
            f"{'Check flagged regions before using.' if hallucinations else ''}"
            f"{'Consider experimental validation for key regions.' if plddt < 70 else ''}"
        ).strip() or "Review the details above for region-specific guidance.",
    }
