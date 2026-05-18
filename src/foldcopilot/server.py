"""FoldCopilot MCP Server — confidence-aware protein structure interpretation.

v0.1: Confidence Copilot MVP
- AFDB lookup (structure metadata, pLDDT, PAE)
- Confidence interpretation with IDR cross-checking
- Hallucination detection against DisProt + MobiDB

v0.2: Foldseek MCP (first MCP wrapper for Foldseek)
- Structural similarity search via Foldseek webserver API
- AFDB <-> Foldseek agent loop with confidence filtering

v0.3: Boltz-2 prediction backend
- Local structure prediction via Boltz-2 CLI (MIT, 20s/GPU)
- License routing (commercial-safe by default)
- Reproducibility manifests with every prediction

v0.4: Ensemble + cross-model disagreement
- Compare structures from multiple backends
- Per-residue agreement classification
- Disagreement span detection — the second moat

v0.6: AlphaMissense + AlphaFill
- Missense variant pathogenicity landscape
- Cofactor/ligand transplantation from experimental structures
- Combined annotation pipeline

v0.7: Therapeutic vertical packs
- Antibody Pack (CDR identification, interface confidence)
- Kinase Pack (ATP site, DFG motif, KLIFS cross-reference)
- GPCR Pack (TM helix confidence, membrane orientation context)

v0.8: Education mode
- Plain-language explanations of pLDDT, PAE, hallucinations
- Targeted at wet-lab biologists, not bioinformaticians

v0.9: Benchmarking harness
- Evaluate predictions against experimental structures
- GDT-TS, CA-RMSD, pLDDT calibration per bucket
- Batch evaluation and publication-ready reports
"""

from __future__ import annotations

from typing import Annotated

from fastmcp import FastMCP

from foldcopilot.tools import (
    afdb, annotations, benchmarks, confidence, education, ensemble,
    foldseek, predict, verticals,
)

mcp = FastMCP(
    "FoldCopilot",
    instructions=(
        "FoldCopilot is a confidence-aware interpretation copilot for protein "
        "structure predictions. Use these tools to look up AlphaFold DB structures, "
        "assess prediction confidence, and detect potential hallucinations in "
        "disordered regions. Always present confidence results with appropriate "
        "caveats — pLDDT is a model's self-assessment, not ground truth."
    ),
)


# --- AFDB Lookup Tools ---


@mcp.tool()
async def lookup_structure(uniprot_id: str) -> dict:
    """Look up an AlphaFold DB entry by UniProt accession.

    Returns structure metadata including download URLs, model version,
    organism, and gene name.

    Example: lookup_structure("P00520") for human ABL1.
    """
    return await afdb.lookup_structure(uniprot_id)


@mcp.tool()
async def get_plddt_scores(uniprot_id: str) -> dict:
    """Get per-residue pLDDT confidence scores for an AFDB structure.

    Returns scores (0-100) and summary statistics. pLDDT buckets:
    - Very high (>90): High confidence in backbone AND side-chain
    - High (70-90): High confidence in backbone
    - Low (50-70): Low confidence — treat with caution
    - Very low (<50): Should not be interpreted — likely disordered or wrong

    Example: get_plddt_scores("P00520")
    """
    return await afdb.get_plddt(uniprot_id)


@mcp.tool()
async def get_pae_summary(uniprot_id: str) -> dict:
    """Get PAE (Predicted Aligned Error) matrix summary for an AFDB structure.

    PAE measures the expected position error of residue X when aligned on
    residue Y. High PAE (>10 angstroms) between domains suggests their
    relative positioning is unreliable.

    Returns summary statistics — not the full matrix (to keep responses compact).

    Example: get_pae_summary("P00520")
    """
    return await afdb.get_pae(uniprot_id)


# --- The Killer Feature ---


@mcp.tool()
async def assess_confidence(uniprot_id: str) -> dict:
    """Comprehensive confidence assessment for an AlphaFold DB structure.

    This is FoldCopilot's core tool. It combines:
    1. pLDDT score analysis with bucketing and low-confidence span detection
    2. PAE matrix summary for domain/interface reliability
    3. IDR (intrinsically disordered region) cross-check against DisProt and MobiDB
    4. Hallucination warnings where AlphaFold predicts order in known disordered regions

    Returns a structured ConfidenceReport with:
    - Chain-level summaries (mean/median pLDDT, bucket distribution, low-confidence spans)
    - PAE summary (mean, high-error fraction)
    - IDR flags from DisProt and MobiDB
    - Hallucination warnings with severity ratings
    - Standard caveats and disclaimers

    IMPORTANT: ~22% of IDR residues may be falsely predicted as ordered by AlphaFold
    (arXiv 2510.15939). This tool flags such cases.

    Example: assess_confidence("P53_HUMAN") or assess_confidence("P04637")
    """
    return await confidence.assess_confidence(uniprot_id)


# --- Foldseek Structural Search Tools ---


@mcp.tool()
async def foldseek_search(
    pdb_content: Annotated[str, "PDB file content as string"],
    databases: Annotated[
        list[str] | None,
        "Databases to search. Options: afdb50, afdb-swissprot, afdb-proteome, "
        "afdb-uniprot50, pdb100, gmgcl_id, bfmd, cath50, mgnify_esm30. "
        "Default: [afdb50, afdb-swissprot, pdb100]"
    ] = None,
    mode: Annotated[
        str, "Search mode: '3diaa' (fast, default) or 'tmalign' (slower, more sensitive)"
    ] = "3diaa",
    max_hits: Annotated[int, "Maximum number of hits to return"] = 20,
) -> dict:
    """Search for structurally similar proteins using Foldseek.

    Submit a PDB structure and find structural homologs across multiple databases.
    This is the first MCP wrapper for Foldseek.

    Returns hits sorted by E-value, each with: target ID, description,
    TM-score, E-value, sequence identity, aligned length.

    Example: foldseek_search(pdb_content="ATOM 1 CA ...", databases=["pdb100"])
    """
    return await foldseek.search_structure(
        pdb_content, databases=databases, mode=mode, max_hits=max_hits
    )


@mcp.tool()
async def foldseek_search_uniprot(
    uniprot_id: Annotated[str, "UniProt accession (e.g., P00520)"],
    databases: Annotated[list[str] | None, "Databases to search. Default: [afdb50, afdb-swissprot, pdb100]"] = None,
    mode: Annotated[str, "Search mode: '3diaa' (fast) or 'tmalign' (sensitive)"] = "3diaa",
    max_hits: Annotated[int, "Maximum hits to return"] = 20,
) -> dict:
    """Search Foldseek using an AlphaFold DB structure by UniProt ID.

    Convenience tool that fetches the AFDB structure for the given UniProt
    accession and submits it to Foldseek for structural similarity search.

    Example: foldseek_search_uniprot("P00520") to find proteins structurally
    similar to human ABL1.
    """
    return await foldseek.search_by_uniprot(
        uniprot_id, databases=databases, mode=mode, max_hits=max_hits
    )


@mcp.tool()
async def find_confident_homologs(
    uniprot_id: Annotated[str, "UniProt accession to search from"],
    min_plddt: Annotated[float, "Minimum mean pLDDT for target structures"] = 70.0,
    min_tm_score: Annotated[float, "Minimum TM-score for structural similarity"] = 0.5,
    max_hits: Annotated[int, "Maximum hits to return"] = 10,
) -> dict:
    """Find structurally similar proteins with confident AlphaFold structures.

    AFDB <-> Foldseek agent loop (no existing MCP does this):
    1. Fetches AFDB structure for the query UniProt ID
    2. Searches Foldseek for structural homologs
    3. Filters by TM-score threshold
    4. Checks pLDDT confidence of each AFDB target hit
    5. Returns only hits where target structures are confident

    Use this to find reliable structural homologs — avoiding hits where
    the target structure itself is low-confidence.

    Example: find_confident_homologs("P04637", min_plddt=80, min_tm_score=0.6)
    """
    return await foldseek.search_confident_homologs(
        uniprot_id,
        min_plddt=min_plddt,
        min_tm_score=min_tm_score,
        max_hits=max_hits,
    )


# --- Prediction Tools ---


@mcp.tool()
async def predict_structure(
    sequences: Annotated[list[str], "Amino acid sequences to fold (one per chain)"],
    backend: Annotated[str, "Prediction backend: 'boltz2' (MIT), 'openfold3' (Apache-2.0), 'chai1' (Apache-2.0)"] = "boltz2",
    commercial_use: Annotated[bool, "Set True if results will be used commercially. Enforces license routing."] = False,
    recycling_steps: Annotated[int, "Number of recycling steps (Boltz-2)"] = 3,
    sampling_steps: Annotated[int, "Number of diffusion sampling steps (Boltz-2)"] = 200,
    diffusion_samples: Annotated[int, "Number of diffusion samples to generate"] = 1,
    use_msa: Annotated[bool, "Use MSA (multiple sequence alignment). Slower but more accurate."] = True,
    predict_affinity: Annotated[bool, "Predict binding affinity (Boltz-2, requires 2+ chains)"] = False,
) -> dict:
    """Predict protein structure using a supported backend.

    Runs structure prediction locally using user-supplied compute.
    Requires the chosen backend installed with GPU access.

    Backends:
    - boltz2 (MIT): ~20s/GPU, affinity prediction, fast default
    - openfold3 (Apache-2.0): AF3 reproduction, commercial-safe
    - chai1 (Apache-2.0): multi-chain, ligand support

    Every prediction ships with a reproducibility manifest containing:
    model version, weights hash, input hash, parameters, runtime env.

    Output files (PDB/CIF) are saved locally — paths returned in response.
    Raw structures are NOT returned in MCP response (too large).

    Example: predict_structure(["MKFLILLFNILCLFPVLAADNHGVS..."])
    """
    return await predict.predict_structure(
        sequences=sequences,
        backend=backend,
        commercial_use=commercial_use,
        recycling_steps=recycling_steps,
        sampling_steps=sampling_steps,
        diffusion_samples=diffusion_samples,
        use_msa=use_msa,
        predict_affinity=predict_affinity,
    )


@mcp.tool()
def check_backend_status(
    backend: Annotated[str, "Backend to check: 'boltz2', 'openfold3', 'chai1'"] = "boltz2",
) -> dict:
    """Check if a prediction backend is installed and ready.

    Returns installation status, GPU availability, and setup instructions
    if the backend is not yet configured.

    Example: check_backend_status("boltz2")
    """
    return predict.get_backend_status(backend)


@mcp.tool()
def list_prediction_backends() -> dict:
    """List all supported prediction backends with license and status info.

    Shows which backends are implemented, their license type (commercial_ok
    vs non_commercial), and whether they're installed locally.

    License routing: when commercial_use=True, only MIT/Apache-2.0 backends
    are available. AF3 weights are non-commercial — use Boltz-2 or OpenFold3.
    """
    return predict.list_backends()


# --- Ensemble Comparison Tools ---


@mcp.tool()
async def compare_predictions(
    pdb_content_a: Annotated[str, "PDB content from first model prediction"],
    pdb_content_b: Annotated[str, "PDB content from second model prediction"],
    model_a_name: Annotated[str, "Name of first model (e.g., 'boltz2')"] = "model_a",
    model_b_name: Annotated[str, "Name of second model (e.g., 'alphafold3')"] = "model_b",
    plddt_threshold: Annotated[float, "pLDDT threshold for 'confident' classification"] = 70.0,
    distance_threshold: Annotated[float, "CA distance threshold (Angstroms) for structural agreement"] = 3.0,
) -> dict:
    """Compare two structure predictions and detect cross-model disagreement.

    This is FoldCopilot's second moat — no other MCP does this.

    Feeds PDB outputs from two different backends (e.g., Boltz-2 vs AlphaFold 3)
    and produces a detailed ensemble report:

    - Per-residue agreement classification:
      * strong_agree: both confident, structures match (<3A CA distance)
      * moderate_agree: one confident, structures match
      * disagree: both confident but structures DIFFER (>3A) — flag for validation
      * both_uncertain: neither model confident — likely disordered

    - Disagreement spans: contiguous regions where models disagree
    - Consensus spans: regions of strong agreement (high structural confidence)
    - Global metrics: CA-RMSD, pLDDT correlation

    When models agree, trust is higher. When they disagree, validate experimentally.

    Example workflow:
    1. predict_structure(sequences, backend="boltz2") -> get pdb_path_a
    2. Run AF3/OpenFold3 separately -> get pdb_path_b
    3. compare_predictions(pdb_a, pdb_b, "boltz2", "alphafold3")
    """
    return await ensemble.compare_structures(
        pdb_content_a, pdb_content_b,
        model_a_name, model_b_name,
        plddt_threshold, distance_threshold,
    )


@mcp.tool()
async def compare_prediction_files(
    pdb_path_a: Annotated[str, "Path to first PDB file"],
    pdb_path_b: Annotated[str, "Path to second PDB file"],
    model_a_name: Annotated[str, "Name of first model"] = "model_a",
    model_b_name: Annotated[str, "Name of second model"] = "model_b",
    plddt_threshold: Annotated[float, "pLDDT threshold for confidence"] = 70.0,
    distance_threshold: Annotated[float, "CA distance threshold (Angstroms)"] = 3.0,
) -> dict:
    """Compare two PDB files from different prediction backends.

    Same as compare_predictions but reads PDB content from file paths.
    Useful after running predict_structure which returns output file paths.

    Example: compare_prediction_files("~/.cache/.../boltz2/out.pdb",
                                       "~/.cache/.../af3/out.pdb",
                                       "boltz2", "alphafold3")
    """
    return await ensemble.compare_by_paths(
        pdb_path_a, pdb_path_b,
        model_a_name, model_b_name,
        plddt_threshold, distance_threshold,
    )


# --- Annotation Tools (v0.6) ---


@mcp.tool()
async def get_missense_landscape(
    uniprot_id: Annotated[str, "UniProt accession (e.g., P04637)"],
) -> dict:
    """Get AlphaMissense pathogenicity landscape for a protein.

    AlphaMissense predicts pathogenicity of all possible single amino acid
    substitutions. Returns per-residue mean pathogenicity scores and
    variant-level classification (likely_pathogenic / likely_benign / ambiguous).

    Use this to identify mutation-sensitive hotspots in a protein structure.

    Example: get_missense_landscape("P04637") for human p53.
    """
    return await annotations.get_missense_landscape(uniprot_id)


@mcp.tool()
async def get_cofactors(
    uniprot_id: Annotated[str, "UniProt accession"],
) -> dict:
    """Get AlphaFill transplanted cofactors, ligands, and metal ions.

    AlphaFill transplants ligands from experimental PDB structures into
    AlphaFold models based on structural homology. Returns compound type
    (cofactor, metal_ion, nucleotide, ligand), source PDB, RMSD, and identity.

    Use this to understand the biochemical context of a predicted structure.

    Example: get_cofactors("P00520") for human ABL1 (expect ATP-site ligands).
    """
    return await annotations.get_cofactors(uniprot_id)


@mcp.tool()
async def get_full_annotation(
    uniprot_id: Annotated[str, "UniProt accession"],
) -> dict:
    """Comprehensive annotation: AlphaMissense + AlphaFill combined.

    One-prompt pipeline: structure -> cofactor transplantation -> missense
    pathogenicity landscape. Also identifies cofactor-pathogenicity hotspots
    where pathogenic variants cluster near cofactor binding sites.

    Example: get_full_annotation("P00520")
    """
    return await annotations.get_full_annotation(uniprot_id)


# --- Therapeutic Vertical Packs (v0.7) ---


@mcp.tool()
async def analyze_antibody(
    heavy_chain: Annotated[str, "Heavy chain amino acid sequence"],
    light_chain: Annotated[str | None, "Light chain sequence (omit for nanobodies/VHH)"] = None,
    target_uniprot_id: Annotated[str | None, "UniProt ID of the antigen target"] = None,
) -> dict:
    """Antibody Pack — comprehensive antibody structure analysis.

    Domain-specific analysis for antibody engineering:
    - CDR identification (approximate Kabat numbering)
    - CDR-specific confidence warnings (CDR-H3 is expected to be low-confidence)
    - Target antigen confidence assessment (if UniProt ID provided)
    - Recommendations for co-folding and experimental validation

    Supports: IgG heavy+light, VHH/nanobodies (heavy only), Fab fragments.

    Example: analyze_antibody("EVQLVES...", "DIQMTQS...", target_uniprot_id="P01308")
    """
    return await verticals.antibody_analysis(heavy_chain, light_chain, target_uniprot_id)


@mcp.tool()
async def analyze_kinase(
    uniprot_id: Annotated[str, "UniProt accession of a kinase"],
) -> dict:
    """Kinase Pack — kinase-specific structural analysis.

    Domain-specific analysis for kinase drug discovery:
    - Confidence assessment focused on catalytic domain
    - AlphaFill: ATP-site ligand transplantation (ATP, ADP, kinase inhibitors)
    - AlphaMissense: mutation sensitivity of kinase domain
    - DFG motif and activation loop context
    - KLIFS cross-reference recommendations

    Use with predict_affinity=True in predict_structure for binding prediction.

    Example: analyze_kinase("P00519") for human ABL1.
    """
    return await verticals.kinase_analysis(uniprot_id)


@mcp.tool()
async def analyze_gpcr(
    uniprot_id: Annotated[str, "UniProt accession of a GPCR"],
) -> dict:
    """GPCR Pack — G protein-coupled receptor analysis.

    Domain-specific analysis for GPCR pharmacology:
    - Transmembrane helix confidence (should be high; low is unusual)
    - Loop region warnings (ICL3, N/C-termini often disordered)
    - AlphaFill ligand transplantation
    - AlphaMissense pathogenicity
    - Membrane orientation context (TMalphaFold, OPM references)
    - Activation state limitations of AF predictions

    Example: analyze_gpcr("P08172") for human muscarinic M2 receptor.
    """
    return await verticals.gpcr_analysis(uniprot_id)


# --- Education Mode (v0.8) ---


@mcp.tool()
def explain_score(
    score: Annotated[float, "pLDDT score (0-100) to explain"],
) -> dict:
    """Explain a pLDDT score in plain language for wet-lab biologists.

    Returns an actionable interpretation with:
    - What the score means in practical terms
    - What you should and shouldn't do with this region
    - A real-world analogy (crystal structure resolution equivalent)
    - Citation for further reading

    Use this when a biologist asks "what does pLDDT 73 mean?" or
    "can I trust this region?"

    Example: explain_score(73.5)
    """
    return education.explain_plddt(score)


@mcp.tool()
def explain_pae_score(
    mean_pae: Annotated[float, "Mean PAE value in Angstroms"],
    context: Annotated[str, "Context: 'general', 'interface', or 'domain'"] = "general",
) -> dict:
    """Explain a PAE (Predicted Aligned Error) value in plain language.

    Translates the abstract Angstrom value into actionable guidance:
    - Can you trust domain arrangements?
    - Is this protein-protein interface reliable?
    - Should you analyze domains separately?

    Example: explain_pae_score(15.3, context="interface")
    """
    return education.explain_pae(mean_pae, context)


@mcp.tool()
def explain_hallucination(
    af_plddt: Annotated[float, "AlphaFold pLDDT in the flagged region"],
    idr_source: Annotated[str, "Source database: 'disprot' or 'mobidb'"],
    severity: Annotated[str, "Severity: 'high' or 'moderate'"],
) -> dict:
    """Explain a hallucination warning in plain language.

    When assess_confidence flags a region where AlphaFold predicts order
    but IDR databases say disorder, use this tool to get a clear,
    actionable explanation of what went wrong and what to do about it.

    Example: explain_hallucination(85.0, "disprot", "high")
    """
    return education.explain_hallucination_warning(af_plddt, idr_source, severity)


@mcp.tool()
def explain_report(
    confidence_report: Annotated[dict, "Output from assess_confidence tool"],
) -> dict:
    """Translate a full ConfidenceReport into plain-language summary.

    Takes the structured JSON output from assess_confidence and produces
    a human-readable interpretation with:
    - Overall verdict (HIGH / MODERATE / LOW / VERY LOW confidence)
    - Plain-language pLDDT and PAE explanations
    - Hallucination warning summaries
    - Bottom-line recommendation

    Use this after assess_confidence when the user is a wet-lab biologist
    who needs actionable guidance, not raw numbers.

    Example: explain_report(assess_confidence("P04637"))
    """
    return education.explain_confidence_report(confidence_report)


# --- Benchmarking Harness (v0.9) ---


@mcp.tool()
def list_benchmarks() -> dict:
    """List available benchmark datasets for prediction evaluation.

    Returns built-in datasets (DisProt hallucination set, CASP16 monomers)
    and the option to use custom user-provided PDB pairs.
    """
    return benchmarks.list_benchmark_datasets()


@mcp.tool()
def benchmark_prediction(
    predicted_pdb: Annotated[str, "Predicted structure PDB content"],
    reference_pdb: Annotated[str, "Experimental reference PDB content"],
    target_name: Annotated[str, "Name/ID of the target"] = "target",
) -> dict:
    """Evaluate a single prediction against an experimental structure.

    Computes:
    - CA-RMSD (global backbone accuracy)
    - GDT-TS (Global Distance Test — fraction within 1/2/4/8 Angstroms)
    - pLDDT calibration (is high-confidence actually accurate?)
    - Per-residue distance distribution
    - Accuracy breakdown by pLDDT bucket

    Use this to validate predictions before publishing or making decisions.

    Example: benchmark_prediction(boltz2_output, crystal_structure, "my_kinase")
    """
    return benchmarks.evaluate_structure_pair(predicted_pdb, reference_pdb, target_name)


@mcp.tool()
def benchmark_batch(
    pairs: Annotated[
        list[dict],
        "List of {predicted_pdb: str, reference_pdb: str, target_name: str} dicts"
    ],
) -> dict:
    """Evaluate a batch of predictions against experimental structures.

    Runs benchmark_prediction on each pair and computes aggregate statistics:
    - Mean/median RMSD and GDT-TS across all targets
    - Best and worst targets
    - pLDDT calibration across the dataset

    Use this for systematic backend evaluation or paper benchmarks.
    """
    return benchmarks.evaluate_batch(pairs)


@mcp.tool()
def generate_report(
    batch_results: Annotated[dict, "Output from benchmark_batch"],
    dataset_name: Annotated[str, "Dataset name (e.g., 'casp16_monomers', 'custom')"] = "custom",
    backend_name: Annotated[str, "Backend that produced the predictions"] = "unknown",
) -> dict:
    """Generate a publication-ready benchmark report.

    Formats batch evaluation results into a structured report suitable for:
    - JOSS paper supplementary data
    - bioRxiv preprint tables
    - Zenodo dataset metadata
    - Public leaderboard entries

    Includes pLDDT calibration analysis (does pLDDT actually predict accuracy?).
    """
    return benchmarks.generate_benchmark_report(batch_results, dataset_name, backend_name)


def main():
    """Run the FoldCopilot MCP server."""
    mcp.run()


if __name__ == "__main__":
    main()
