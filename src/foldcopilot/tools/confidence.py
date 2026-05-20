"""The killer feature: assess_confidence tool.

Combines AFDB pLDDT + PAE + IDR flagging from DisProt/MobiDB into a
structured ConfidenceReport with hallucination warnings.
"""

from __future__ import annotations

import numpy as np

from foldcopilot.clients import afdb_client, disprot_client, mobidb_client
from foldcopilot.models.confidence import (
    ChainConfidenceSummary,
    ConfidenceBucket,
    ConfidenceReport,
    HallucinationWarning,
    IDRFlag,
    LowConfidenceSpan,
    PAESummary,
)


def _bucket(plddt: float) -> ConfidenceBucket:
    if plddt > 90:
        return ConfidenceBucket.VERY_HIGH
    elif plddt > 70:
        return ConfidenceBucket.HIGH
    elif plddt > 50:
        return ConfidenceBucket.LOW
    return ConfidenceBucket.VERY_LOW


def _find_low_confidence_spans(
    scores: list[float], threshold: float = 70.0, min_length: int = 3
) -> list[LowConfidenceSpan]:
    """Find contiguous spans where pLDDT < threshold."""
    spans = []
    in_span = False
    start = 0

    for i, s in enumerate(scores):
        if s < threshold:
            if not in_span:
                start = i
                in_span = True
        else:
            if in_span:
                length = i - start
                if length >= min_length:
                    span_scores = scores[start:i]
                    mean = float(np.mean(span_scores))
                    spans.append(LowConfidenceSpan(
                        start=start + 1,  # 1-indexed
                        end=i,            # 1-indexed, inclusive
                        mean_plddt=round(mean, 1),
                        bucket=_bucket(mean),
                        length=length,
                    ))
                in_span = False

    # Handle span at end of sequence
    if in_span:
        length = len(scores) - start
        if length >= min_length:
            span_scores = scores[start:]
            mean = float(np.mean(span_scores))
            spans.append(LowConfidenceSpan(
                start=start + 1,
                end=len(scores),
                mean_plddt=round(mean, 1),
                bucket=_bucket(mean),
                length=length,
            ))

    return spans


def _check_hallucinations(
    scores: list[float],
    idr_regions: list[dict],
    source: str,
) -> list[HallucinationWarning]:
    """Check for regions where AF predicts order but IDR databases say disorder."""
    warnings = []
    for region in idr_regions:
        start = region["start"] - 1  # convert to 0-indexed
        end = region["end"]          # exclusive
        start = max(0, start)
        end = min(len(scores), end)

        if start >= end:
            continue

        region_scores = scores[start:end]
        mean_plddt = float(np.mean(region_scores))

        # Flag if AF is confident (pLDDT > 50) in a region that's known disordered
        if mean_plddt > 50:
            severity = "high" if mean_plddt > 70 else "moderate"
            warnings.append(HallucinationWarning(
                start=region["start"],
                end=region["end"],
                af_mean_plddt=round(mean_plddt, 1),
                idr_source=source,
                idr_annotation=region.get("annotation"),
                severity=severity,
            ))

    return warnings


async def assess_confidence(uniprot_id: str) -> dict:
    """Assess confidence of an AlphaFold DB structure.

    Combines pLDDT scores, PAE matrix analysis, and IDR cross-checking
    against DisProt and MobiDB to produce a comprehensive confidence report
    with hallucination warnings.

    This is the core value proposition of FoldCopilot.
    """
    import httpx

    async with httpx.AsyncClient(timeout=60) as client:
        # Fetch all data in parallel via gather-like sequential calls
        # (httpx doesn't natively support gather, but we reuse the client)
        scores = await afdb_client.get_plddt_scores(uniprot_id, client=client)
        pae_matrix = await afdb_client.get_pae_matrix(uniprot_id, client=client)
        meta = await afdb_client.get_prediction_metadata(uniprot_id, client=client)

        # IDR data — these may 404 gracefully
        disprot_regions = await disprot_client.get_disprot_regions(
            uniprot_id, client=client
        )
        mobidb_regions = await mobidb_client.get_mobidb_regions(
            uniprot_id, client=client
        )

    arr = np.array(scores)

    # Chain summary (single chain for AFDB monomers)
    bucket_dist = {
        "very_high": float(np.mean(arr > 90)),
        "high": float(np.mean((arr > 70) & (arr <= 90))),
        "low": float(np.mean((arr > 50) & (arr <= 70))),
        "very_low": float(np.mean(arr <= 50)),
    }

    chain_summary = ChainConfidenceSummary(
        chain_id="A",
        mean_plddt=round(float(np.mean(arr)), 1),
        median_plddt=round(float(np.median(arr)), 1),
        bucket_distribution={k: round(v, 3) for k, v in bucket_dist.items()},
        low_confidence_spans=_find_low_confidence_spans(scores),
        residue_count=len(scores),
    )

    # PAE summary
    pae_summary = None
    if pae_matrix is not None:
        pae_summary = PAESummary(
            mean_pae=round(float(np.mean(pae_matrix)), 1),
            median_pae=round(float(np.median(pae_matrix)), 1),
            interface_pae=None,  # single chain, no interface
            high_error_fraction=round(float(np.mean(pae_matrix > 10)), 3),
        )

    # IDR flags
    idr_flags = []
    for r in disprot_regions:
        idr_flags.append(IDRFlag(
            source="disprot",
            start=r["start"],
            end=r["end"],
            annotation=r.get("annotation"),
        ))
    for r in mobidb_regions:
        idr_flags.append(IDRFlag(
            source="mobidb",
            start=r["start"],
            end=r["end"],
            annotation=r.get("annotation"),
        ))

    # Hallucination detection
    hallucination_warnings = []
    hallucination_warnings.extend(
        _check_hallucinations(scores, disprot_regions, "disprot")
    )
    hallucination_warnings.extend(
        _check_hallucinations(scores, mobidb_regions, "mobidb")
    )

    # IDR strategy note: AF2 vs AF3 for disordered regions
    idr_strategy_note = (
        "For IDR identification, AF2 is preferred over AF3 — AF2 avoids structural "
        "hallucinations in disordered regions while AF3 creates them (PubMed 41454828). "
        "Consider cross-referencing AF2 pLDDT < 50 as a complementary disorder signal."
    )

    # Benchmarking context
    confidence_interpretation_note = (
        "CASP16 assessment (2026) confirms monomer fold prediction is largely solved. "
        "Remaining challenges: multimers (<25% high quality), irregular secondary "
        "structures, and interaction-induced conformational changes."
    )

    report = ConfidenceReport(
        uniprot_id=uniprot_id,
        source="afdb",
        model_version=str(meta.get("latestVersion", "unknown")),
        chain_summaries=[chain_summary],
        overall_mean_plddt=chain_summary.mean_plddt,
        overall_median_plddt=chain_summary.median_plddt,
        pae_summary=pae_summary,
        idr_flags=idr_flags,
        hallucination_warnings=hallucination_warnings,
        idr_strategy_note=idr_strategy_note,
        confidence_interpretation_note=confidence_interpretation_note,
    )
    report.add_standard_caveats()

    return report.model_dump()
