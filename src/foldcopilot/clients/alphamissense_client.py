"""AlphaMissense API client — missense variant pathogenicity predictions.

AlphaMissense (Cheng et al., Science 2023) predicts pathogenicity of all
possible single amino acid substitutions across the human proteome.
Data available via AFDB 2025 integration and Google Cloud.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import httpx

AFDB_BASE = "https://alphafold.ebi.ac.uk/api"
ALPHAMISSENSE_BASE = "https://alphafold.ebi.ac.uk/files"

_CACHE_DIR = Path.home() / ".cache" / "foldcopilot" / "alphamissense"


def _cache_path(uniprot_id: str) -> Path:
    h = hashlib.sha256(uniprot_id.encode()).hexdigest()[:12]
    return _CACHE_DIR / f"{uniprot_id}_{h}.json"


def _read_cache(path: Path) -> Any | None:
    if path.exists():
        return json.loads(path.read_text())
    return None


def _write_cache(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


async def get_missense_predictions(
    uniprot_id: str, *, client: httpx.AsyncClient | None = None
) -> dict:
    """Fetch AlphaMissense pathogenicity predictions for a UniProt accession.

    Returns per-residue pathogenicity landscape and variant-level scores.
    """
    cp = _cache_path(uniprot_id)
    cached = _read_cache(cp)
    if cached is not None:
        return cached

    own_client = client is None
    if own_client:
        client = httpx.AsyncClient(timeout=60)
    try:
        # Fetch via AFDB AlphaMissense endpoint
        entry_id = f"AF-{uniprot_id}-F1"
        url = f"{ALPHAMISSENSE_BASE}/{entry_id}-aa-substitutions.csv"

        resp = await client.get(url)
        if resp.status_code == 404:
            return {
                "uniprot_id": uniprot_id,
                "available": False,
                "message": "AlphaMissense data not available for this protein.",
            }
        resp.raise_for_status()

        # Parse CSV: columns are protein_variant, am_pathogenicity, am_class
        lines = resp.text.strip().splitlines()
        header_idx = 0
        for i, line in enumerate(lines):
            if line.startswith("#"):
                header_idx = i + 1
                continue
            if "protein_variant" in line.lower() or "am_pathogenicity" in line.lower():
                header_idx = i + 1
                break

        variants = []
        residue_scores: dict[int, list[float]] = {}

        for line in lines[header_idx:]:
            parts = line.strip().split(",")
            if len(parts) < 3:
                continue
            variant = parts[0].strip()
            try:
                score = float(parts[1].strip())
            except ValueError:
                continue
            classification = parts[2].strip() if len(parts) > 2 else ""

            # Extract residue position from variant (e.g., "M1A" -> position 1)
            pos = _extract_position(variant)
            if pos is not None:
                residue_scores.setdefault(pos, []).append(score)

            variants.append({
                "variant": variant,
                "pathogenicity_score": score,
                "classification": classification,
            })

        # Compute per-residue mean pathogenicity
        residue_landscape = {}
        for pos, scores in sorted(residue_scores.items()):
            import numpy as np
            mean_score = float(np.mean(scores))
            residue_landscape[str(pos)] = {
                "mean_pathogenicity": round(mean_score, 4),
                "classification": _classify_score(mean_score),
                "n_variants": len(scores),
            }

        result = {
            "uniprot_id": uniprot_id,
            "available": True,
            "total_variants": len(variants),
            "residue_count": len(residue_landscape),
            "residue_landscape": residue_landscape,
            "pathogenic_fraction": _compute_fraction(variants, "pathogenic"),
            "benign_fraction": _compute_fraction(variants, "benign"),
            "ambiguous_fraction": _compute_fraction(variants, "ambiguous"),
        }

        _write_cache(cp, result)
        return result

    finally:
        if own_client:
            await client.aclose()


def _extract_position(variant: str) -> int | None:
    """Extract residue position from variant string like 'M1A', 'G123R'."""
    digits = ""
    started = False
    for c in variant:
        if c.isdigit():
            digits += c
            started = True
        elif started:
            break
    try:
        return int(digits) if digits else None
    except ValueError:
        return None


def _classify_score(score: float) -> str:
    """Classify AlphaMissense pathogenicity score."""
    if score >= 0.564:
        return "likely_pathogenic"
    elif score <= 0.34:
        return "likely_benign"
    return "ambiguous"


def _compute_fraction(variants: list[dict], classification: str) -> float:
    """Compute fraction of variants with a given classification."""
    if not variants:
        return 0.0
    count = sum(1 for v in variants if classification in v.get("classification", "").lower())
    return round(count / len(variants), 4)
