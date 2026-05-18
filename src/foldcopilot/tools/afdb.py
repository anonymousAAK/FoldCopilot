"""AFDB lookup tools — fetch structure metadata, pLDDT, and PAE from AlphaFold DB."""

from __future__ import annotations

from foldcopilot.clients import afdb_client


async def lookup_structure(uniprot_id: str) -> dict:
    """Look up an AlphaFold DB entry by UniProt accession.

    Returns metadata including structure URLs, model version, sequence length,
    and available confidence data.
    """
    meta = await afdb_client.get_prediction_metadata(uniprot_id)
    return {
        "uniprot_id": uniprot_id,
        "entry_id": meta.get("entryId"),
        "gene": meta.get("gene"),
        "organism": meta.get("organismScientificName"),
        "sequence_length": meta.get("uniprotEnd", 0) - meta.get("uniprotStart", 0) + 1,
        "model_version": meta.get("latestVersion"),
        "pdb_url": meta.get("pdbUrl"),
        "cif_url": meta.get("cifUrl"),
        "pae_image_url": meta.get("paeImageUrl"),
        "pae_json_url": meta.get("paeDocUrl"),
    }


async def get_plddt(uniprot_id: str) -> dict:
    """Fetch per-residue pLDDT scores for an AFDB entry.

    Returns scores and basic statistics.
    """
    import numpy as np

    scores = await afdb_client.get_plddt_scores(uniprot_id)
    arr = np.array(scores)
    return {
        "uniprot_id": uniprot_id,
        "residue_count": len(scores),
        "mean_plddt": float(np.mean(arr)),
        "median_plddt": float(np.median(arr)),
        "min_plddt": float(np.min(arr)),
        "max_plddt": float(np.max(arr)),
        "very_high_fraction": float(np.mean(arr > 90)),
        "high_fraction": float(np.mean((arr > 70) & (arr <= 90))),
        "low_fraction": float(np.mean((arr > 50) & (arr <= 70))),
        "very_low_fraction": float(np.mean(arr <= 50)),
        "scores": scores,
    }


async def get_pae(uniprot_id: str) -> dict:
    """Fetch PAE (Predicted Aligned Error) matrix summary for an AFDB entry.

    Returns summary statistics. Full matrix not returned to keep response compact.
    """
    import numpy as np

    matrix = await afdb_client.get_pae_matrix(uniprot_id)
    if matrix is None:
        return {
            "uniprot_id": uniprot_id,
            "available": False,
            "message": "PAE data not available for this entry.",
        }

    return {
        "uniprot_id": uniprot_id,
        "available": True,
        "matrix_size": list(matrix.shape),
        "mean_pae": float(np.mean(matrix)),
        "median_pae": float(np.median(matrix)),
        "max_pae": float(np.max(matrix)),
        "high_error_fraction": float(np.mean(matrix > 10)),
        "very_high_error_fraction": float(np.mean(matrix > 20)),
    }
