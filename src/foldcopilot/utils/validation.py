"""Input validation utilities for FoldCopilot.

Validates at system boundaries: user-provided sequences, UniProt IDs, PDB content.
"""

from __future__ import annotations

import re

# Standard amino acid one-letter codes
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
# Including ambiguous codes
AMINO_ACIDS_EXTENDED = AMINO_ACIDS | set("BJOUXZ")

# UniProt accession pattern: 6-10 alphanumeric
UNIPROT_PATTERN = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")


class ValidationError(Exception):
    """Raised when input validation fails."""
    pass


def validate_sequence(sequence: str, max_length: int = 5000) -> str:
    """Validate and normalize a protein sequence.

    Returns the cleaned sequence (uppercase, whitespace stripped).
    Raises ValidationError on invalid input.
    """
    if not sequence or not sequence.strip():
        raise ValidationError("Sequence cannot be empty.")

    cleaned = re.sub(r"\s+", "", sequence.upper())

    if len(cleaned) < 10:
        raise ValidationError(
            f"Sequence too short ({len(cleaned)} residues). "
            f"Minimum 10 residues required for meaningful prediction."
        )

    if len(cleaned) > max_length:
        raise ValidationError(
            f"Sequence too long ({len(cleaned)} residues). "
            f"Maximum {max_length} residues. Split into domains or use "
            f"a chunking strategy."
        )

    invalid = set(cleaned) - AMINO_ACIDS_EXTENDED
    if invalid:
        raise ValidationError(
            f"Invalid characters in sequence: {sorted(invalid)}. "
            f"Expected standard amino acid one-letter codes."
        )

    return cleaned


def validate_sequences(sequences: list[str], max_chains: int = 10) -> list[str]:
    """Validate a list of protein sequences for multi-chain prediction."""
    if not sequences:
        raise ValidationError("At least one sequence is required.")

    if len(sequences) > max_chains:
        raise ValidationError(
            f"Too many chains ({len(sequences)}). Maximum {max_chains}."
        )

    return [validate_sequence(seq) for seq in sequences]


def validate_uniprot_id(uniprot_id: str) -> str:
    """Validate a UniProt accession ID.

    Accepts standard UniProt accessions (P12345, A0A0A0MRZ7) and
    common aliases like entry names (P53_HUMAN).
    """
    if not uniprot_id or not uniprot_id.strip():
        raise ValidationError("UniProt ID cannot be empty.")

    cleaned = uniprot_id.strip().upper()

    # Allow entry names like P53_HUMAN
    if "_" in cleaned:
        parts = cleaned.split("_")
        if len(parts) == 2 and all(p.isalnum() for p in parts):
            return cleaned

    # Standard accession validation
    if not cleaned.isalnum():
        raise ValidationError(
            f"Invalid UniProt ID: '{uniprot_id}'. "
            f"Expected alphanumeric accession (e.g., P04637, A0A0A0MRZ7)."
        )

    if len(cleaned) < 6 or len(cleaned) > 10:
        raise ValidationError(
            f"Invalid UniProt ID length: '{uniprot_id}'. "
            f"Expected 6-10 characters."
        )

    return cleaned


def validate_pdb_content(pdb_content: str) -> str:
    """Basic validation of PDB content."""
    if not pdb_content or not pdb_content.strip():
        raise ValidationError("PDB content cannot be empty.")

    lines = pdb_content.strip().splitlines()
    atom_lines = [l for l in lines if l.startswith(("ATOM", "HETATM"))]

    if not atom_lines:
        raise ValidationError(
            "No ATOM or HETATM records found in PDB content. "
            "Expected standard PDB format."
        )

    return pdb_content
