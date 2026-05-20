"""Data models for confidence reports and structure assessment."""

from __future__ import annotations

from enum import Enum
from pydantic import BaseModel, Field


class ConfidenceBucket(str, Enum):
    VERY_HIGH = "very_high"  # pLDDT > 90
    HIGH = "high"            # 70 < pLDDT <= 90
    LOW = "low"              # 50 < pLDDT <= 70
    VERY_LOW = "very_low"    # pLDDT <= 50


class ResidueConfidence(BaseModel):
    residue_index: int
    residue_name: str | None = None
    plddt: float
    bucket: ConfidenceBucket


class LowConfidenceSpan(BaseModel):
    start: int
    end: int
    mean_plddt: float
    bucket: ConfidenceBucket
    length: int


class IDRFlag(BaseModel):
    """A region flagged as intrinsically disordered by DisProt/MobiDB."""

    source: str  # "disprot" or "mobidb"
    start: int
    end: int
    annotation: str | None = None


class HallucinationWarning(BaseModel):
    """Region where AF predicts order but IDR databases say disorder."""

    start: int
    end: int
    af_mean_plddt: float
    idr_source: str
    idr_annotation: str | None = None
    severity: str  # "high" if pLDDT > 70 on a known IDR


class PAESummary(BaseModel):
    mean_pae: float | None = None
    median_pae: float | None = None
    interface_pae: float | None = None
    high_error_fraction: float | None = None  # fraction of pairs with PAE > 10


class ChainConfidenceSummary(BaseModel):
    chain_id: str | None = None
    mean_plddt: float
    median_plddt: float
    bucket_distribution: dict[str, float]  # bucket -> fraction of residues
    low_confidence_spans: list[LowConfidenceSpan]
    residue_count: int


class ConfidenceReport(BaseModel):
    """The killer feature: structured confidence assessment for a protein structure."""

    uniprot_id: str | None = None
    source: str  # "afdb", "prediction", etc.
    model_version: str | None = None

    chain_summaries: list[ChainConfidenceSummary]
    overall_mean_plddt: float
    overall_median_plddt: float

    pae_summary: PAESummary | None = None

    idr_flags: list[IDRFlag] = Field(default_factory=list)
    hallucination_warnings: list[HallucinationWarning] = Field(default_factory=list)

    caveats: list[str] = Field(default_factory=list)

    af2_disorder_regions: list[dict] | None = None
    idr_strategy_note: str | None = None
    confidence_interpretation_note: str | None = None

    def add_standard_caveats(self) -> None:
        self.caveats.extend([
            "pLDDT > 70 does not guarantee correctness; it indicates the model's "
            "own confidence. Always validate critical regions experimentally.",
            "AlphaFold can hallucinate ordered structure in intrinsically disordered "
            "regions (IDRs). ~22% of IDR residues may be falsely predicted as ordered "
            "(arXiv 2510.15939).",
            "PAE values indicate predicted alignment error between residue pairs. "
            "High PAE (>10 A) at interfaces suggests unreliable domain/chain positioning.",
            "This assessment is for research use only. Do not use for clinical decisions.",
        ])
