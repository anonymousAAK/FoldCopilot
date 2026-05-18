"""Data models for cross-model ensemble comparison and disagreement detection."""

from __future__ import annotations

from enum import Enum

from pydantic import BaseModel, Field


class AgreementLevel(str, Enum):
    STRONG_AGREE = "strong_agree"       # Both confident, structures match
    MODERATE_AGREE = "moderate_agree"    # One confident, structures match
    DISAGREE = "disagree"               # Both confident, structures differ
    BOTH_UNCERTAIN = "both_uncertain"   # Neither model confident
    SINGLE_MODEL = "single_model"       # Only one model available


class ResidueDisagreement(BaseModel):
    """Per-residue disagreement between two models."""

    residue_index: int  # 1-indexed
    residue_name: str | None = None
    model_a: str
    model_b: str
    plddt_a: float
    plddt_b: float
    ca_distance: float | None = None  # Angstroms between CA atoms
    agreement: AgreementLevel


class DisagreementSpan(BaseModel):
    """Contiguous region of disagreement between models."""

    start: int  # 1-indexed
    end: int    # 1-indexed, inclusive
    length: int
    mean_ca_distance: float | None = None
    mean_plddt_a: float
    mean_plddt_b: float
    agreement: AgreementLevel
    interpretation: str


class ModelSummary(BaseModel):
    """Summary of a single model's prediction in ensemble context."""

    model_name: str
    mean_plddt: float
    median_plddt: float
    residue_count: int
    pdb_path: str | None = None


class EnsembleReport(BaseModel):
    """Cross-model comparison report — the ensemble disagreement feature."""

    models: list[ModelSummary]
    residue_count: int

    # Global metrics
    mean_ca_rmsd: float | None = None
    plddt_correlation: float | None = None
    agreement_fraction: float  # fraction of residues where models agree

    # Per-residue breakdown
    strong_agree_fraction: float
    moderate_agree_fraction: float
    disagree_fraction: float
    both_uncertain_fraction: float

    # Flagged regions
    disagreement_spans: list[DisagreementSpan] = Field(default_factory=list)
    high_confidence_consensus_spans: list[DisagreementSpan] = Field(default_factory=list)

    interpretation: str
    caveats: list[str] = Field(default_factory=list)

    def add_standard_caveats(self) -> None:
        self.caveats.extend([
            "Cross-model agreement increases confidence but does not guarantee "
            "correctness — models share training data biases.",
            "Disagreement regions warrant experimental validation. Both models "
            "may be wrong in the same way for some fold classes.",
            "RMSD comparison assumes aligned sequences. Insertions/deletions "
            "inflate RMSD at boundaries.",
            "This assessment is for research use only.",
        ])
