"""Foldseek structural similarity search tools.

First MCP wrapper for Foldseek. Supports:
- Direct PDB search
- UniProt ID search (fetches AFDB structure, then searches Foldseek)
- AFDB <-> Foldseek agent loop with confidence filtering
"""

from __future__ import annotations

import httpx

from foldcopilot.clients import afdb_client, foldseek_client


async def search_structure(
    pdb_content: str,
    databases: list[str] | None = None,
    mode: str = "3diaa",
    max_hits: int = 20,
) -> dict:
    """Search Foldseek with a PDB structure.

    Returns top structural hits across selected databases.
    """
    async with httpx.AsyncClient(timeout=120) as client:
        raw = await foldseek_client.search_pdb(
            pdb_content, databases=databases, mode=mode, client=client
        )

    hits = foldseek_client.parse_alignments(raw, max_hits=max_hits)
    return {
        "hit_count": len(hits),
        "mode": mode,
        "databases": databases or ["afdb50", "afdb-swissprot", "pdb100"],
        "hits": hits,
    }


async def search_by_uniprot(
    uniprot_id: str,
    databases: list[str] | None = None,
    mode: str = "3diaa",
    max_hits: int = 20,
) -> dict:
    """Fetch an AFDB structure by UniProt ID, then search Foldseek for structural homologs.

    Convenience tool that chains AFDB lookup -> Foldseek search.
    """
    async with httpx.AsyncClient(timeout=120) as client:
        # Get PDB from AFDB
        meta = await afdb_client.get_prediction_metadata(uniprot_id, client=client)
        pdb_url = meta.get("pdbUrl")
        if not pdb_url:
            entry_id = meta.get("entryId", f"AF-{uniprot_id}-F1")
            pdb_url = f"https://alphafold.ebi.ac.uk/files/{entry_id}-model_v4.pdb"

        resp = await client.get(pdb_url)
        resp.raise_for_status()
        pdb_content = resp.text

        # Search Foldseek
        raw = await foldseek_client.search_pdb(
            pdb_content, databases=databases, mode=mode, client=client
        )

    hits = foldseek_client.parse_alignments(raw, max_hits=max_hits)
    return {
        "query_uniprot_id": uniprot_id,
        "hit_count": len(hits),
        "mode": mode,
        "databases": databases or ["afdb50", "afdb-swissprot", "pdb100"],
        "hits": hits,
    }


async def search_confident_homologs(
    uniprot_id: str,
    min_plddt: float = 70.0,
    min_tm_score: float = 0.5,
    max_hits: int = 10,
) -> dict:
    """AFDB <-> Foldseek agent loop: find structurally similar proteins
    that have confident interfaces.

    1. Fetches AFDB structure for uniprot_id
    2. Searches Foldseek for structural homologs
    3. Filters hits by TM-score threshold
    4. For top hits from AFDB, checks their pLDDT confidence
    5. Returns only hits where both query and target have confident structures

    This is the AFDB<->Foldseek loop described in the plan — no existing MCP does this.
    """
    async with httpx.AsyncClient(timeout=180) as client:
        # Step 1: Get query pLDDT
        query_scores = await afdb_client.get_plddt_scores(uniprot_id, client=client)
        import numpy as np
        query_mean_plddt = float(np.mean(query_scores))

        # Step 2: Get PDB and search Foldseek
        meta = await afdb_client.get_prediction_metadata(uniprot_id, client=client)
        pdb_url = meta.get("pdbUrl")
        if not pdb_url:
            entry_id = meta.get("entryId", f"AF-{uniprot_id}-F1")
            pdb_url = f"https://alphafold.ebi.ac.uk/files/{entry_id}-model_v4.pdb"

        resp = await client.get(pdb_url)
        resp.raise_for_status()
        pdb_content = resp.text

        raw = await foldseek_client.search_pdb(
            pdb_content,
            databases=["afdb-swissprot"],
            mode="3diaa",
            client=client,
        )
        all_hits = foldseek_client.parse_alignments(raw, max_hits=50)

        # Step 3: Filter by TM-score
        tm_filtered = [
            h for h in all_hits
            if (h.get("tm_score") or 0) >= min_tm_score
        ]

        # Step 4: Check confidence of AFDB hits
        confident_hits = []
        for hit in tm_filtered[:max_hits * 2]:  # check extra to account for filtering
            target_id = _extract_uniprot_from_target(hit.get("target", ""))
            if not target_id:
                # Non-AFDB hit, include but mark confidence as unknown
                hit["target_confidence"] = "unknown"
                confident_hits.append(hit)
                continue

            try:
                target_scores = await afdb_client.get_plddt_scores(
                    target_id, client=client
                )
                target_mean = float(np.mean(target_scores))
                hit["target_uniprot_id"] = target_id
                hit["target_mean_plddt"] = round(target_mean, 1)
                hit["target_confidence"] = (
                    "high" if target_mean >= min_plddt else "low"
                )
                if target_mean >= min_plddt:
                    confident_hits.append(hit)
            except Exception:
                hit["target_confidence"] = "fetch_error"
                confident_hits.append(hit)

            if len(confident_hits) >= max_hits:
                break

    return {
        "query_uniprot_id": uniprot_id,
        "query_mean_plddt": round(query_mean_plddt, 1),
        "min_plddt_filter": min_plddt,
        "min_tm_score_filter": min_tm_score,
        "total_foldseek_hits": len(all_hits),
        "tm_filtered_hits": len(tm_filtered),
        "confident_hit_count": len(confident_hits),
        "hits": confident_hits[:max_hits],
    }


def _extract_uniprot_from_target(target: str) -> str | None:
    """Try to extract a UniProt accession from a Foldseek target identifier.

    AFDB targets look like: AF-P12345-F1 or AF-P12345-F1-model_v4
    """
    if not target:
        return None
    # Pattern: AF-{UNIPROT}-F{N}
    if target.startswith("AF-"):
        parts = target.split("-")
        if len(parts) >= 3:
            return parts[1]
    # Could also be raw UniProt accession
    # UniProt accession pattern: 6-10 alphanumeric chars
    if 6 <= len(target) <= 10 and target.isalnum():
        return target
    return None
