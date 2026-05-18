"""AlphaFill API client — cofactor and ligand transplantation.

AlphaFill (Hekkelman et al., Nature Methods 2023) transplants ligands,
cofactors, and ions from experimental PDB structures into AlphaFold models
based on structural homology.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import httpx

ALPHAFILL_API = "https://alphafill.eu/v1"

_CACHE_DIR = Path.home() / ".cache" / "foldcopilot" / "alphafill"


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


async def get_alphafill_data(
    uniprot_id: str, *, client: httpx.AsyncClient | None = None
) -> dict:
    """Fetch AlphaFill transplanted ligands/cofactors for a UniProt accession.

    Returns transplanted compounds with their source PDB, identity,
    RMSD, and confidence metrics.
    """
    cp = _cache_path(uniprot_id)
    cached = _read_cache(cp)
    if cached is not None:
        return cached

    own_client = client is None
    if own_client:
        client = httpx.AsyncClient(timeout=60)
    try:
        entry_id = f"AF-{uniprot_id}-F1"
        resp = await client.get(f"{ALPHAFILL_API}/aff/{entry_id}")

        if resp.status_code == 404:
            return {
                "uniprot_id": uniprot_id,
                "available": False,
                "message": "AlphaFill data not available for this protein.",
            }
        resp.raise_for_status()
        data = resp.json()

        # Parse transplanted compounds
        transplants = []
        for hit in data.get("hits", []):
            transplant = {
                "compound_id": hit.get("compound_id", ""),
                "compound_name": hit.get("compound_name", ""),
                "compound_type": _classify_compound(hit),
                "source_pdb": hit.get("pdb_id", ""),
                "source_chain": hit.get("chain_id", ""),
                "identity": hit.get("identity"),
                "rmsd": hit.get("rmsd"),
                "clash_count": hit.get("clash_count", 0),
                "analogue_count": hit.get("analogue_count", 0),
            }
            transplants.append(transplant)

        # Summarize by type
        type_counts: dict[str, int] = {}
        for t in transplants:
            ct = t["compound_type"]
            type_counts[ct] = type_counts.get(ct, 0) + 1

        result = {
            "uniprot_id": uniprot_id,
            "available": True,
            "total_transplants": len(transplants),
            "type_summary": type_counts,
            "transplants": transplants,
            "cif_url": data.get("cif_url"),
        }

        _write_cache(cp, result)
        return result

    finally:
        if own_client:
            await client.aclose()


def _classify_compound(hit: dict) -> str:
    """Classify a transplanted compound by type."""
    compound_id = hit.get("compound_id", "").upper()
    compound_name = (hit.get("compound_name") or "").lower()

    # Common cofactors
    cofactors = {
        "NAD", "NAP", "FAD", "FMN", "COA", "SAM", "SAH", "PLP",
        "TPP", "HEM", "HEC", "BCL", "CLA", "FES", "SF4", "F3S",
    }
    if compound_id in cofactors:
        return "cofactor"

    # Metal ions
    metals = {
        "ZN", "MG", "CA", "FE", "MN", "CU", "CO", "NI", "CD",
        "NA", "K", "FE2", "MG2", "CA2", "ZN2",
    }
    if compound_id in metals:
        return "metal_ion"

    # Nucleotides
    nucleotides = {
        "ATP", "ADP", "AMP", "GTP", "GDP", "GMP",
        "CTP", "UTP", "TTP",
    }
    if compound_id in nucleotides:
        return "nucleotide"

    # Water / solvent
    if compound_id in ("HOH", "DOD", "WAT"):
        return "solvent"

    # Check name-based hints
    if any(kw in compound_name for kw in ("inhibitor", "drug", "substrate")):
        return "ligand"

    return "ligand"
