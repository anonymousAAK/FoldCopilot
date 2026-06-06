"""Tests for Foldseek client and tools."""


from foldcopilot.clients.foldseek_client import parse_alignments
from foldcopilot.tools.foldseek import _extract_uniprot_from_target


class TestParseAlignments:
    def test_empty_results(self):
        hits = parse_alignments({})
        assert hits == []

    def test_empty_results_list(self):
        hits = parse_alignments({"results": []})
        assert hits == []

    def test_single_hit(self):
        raw = {
            "results": [
                {
                    "db": "pdb100",
                    "alignments": [
                        {
                            "target": "1abc_A",
                            "tDescription": "Some protein",
                            "prob": 99.0,
                            "eval": 1e-50,
                            "score": 500,
                            "tmScore": 0.95,
                            "alnLength": 200,
                            "seqId": 0.85,
                            "qStartPos": 1,
                            "qEndPos": 200,
                            "tStartPos": 1,
                            "tEndPos": 200,
                        }
                    ],
                }
            ]
        }
        hits = parse_alignments(raw)
        assert len(hits) == 1
        assert hits[0]["database"] == "pdb100"
        assert hits[0]["target"] == "1abc_A"
        assert hits[0]["tm_score"] == 0.95
        assert hits[0]["evalue"] == 1e-50

    def test_multiple_databases(self):
        raw = {
            "results": [
                {
                    "db": "pdb100",
                    "alignments": [
                        {"target": "1abc_A", "eval": 1e-50, "tmScore": 0.9},
                    ],
                },
                {
                    "db": "afdb50",
                    "alignments": [
                        {"target": "AF-P12345-F1", "eval": 1e-30, "tmScore": 0.8},
                    ],
                },
            ]
        }
        hits = parse_alignments(raw)
        assert len(hits) == 2
        # Should be sorted by evalue
        assert hits[0]["evalue"] == 1e-50

    def test_max_hits_limit(self):
        raw = {
            "results": [
                {
                    "db": "pdb100",
                    "alignments": [
                        {"target": f"hit_{i}", "eval": float(i), "tmScore": 0.5}
                        for i in range(50)
                    ],
                }
            ]
        }
        hits = parse_alignments(raw, max_hits=5)
        assert len(hits) == 5

    def test_sort_by_evalue_then_tm(self):
        raw = {
            "results": [
                {
                    "db": "pdb100",
                    "alignments": [
                        {"target": "a", "eval": 1e-10, "tmScore": 0.5},
                        {"target": "b", "eval": 1e-10, "tmScore": 0.9},
                        {"target": "c", "eval": 1e-50, "tmScore": 0.7},
                    ],
                }
            ]
        }
        hits = parse_alignments(raw)
        assert hits[0]["target"] == "c"  # best evalue
        assert hits[1]["target"] == "b"  # same evalue as a, better TM


class TestExtractUniprot:
    def test_afdb_format(self):
        assert _extract_uniprot_from_target("AF-P12345-F1") == "P12345"

    def test_afdb_with_version(self):
        assert _extract_uniprot_from_target("AF-Q9NZJ5-F1-model_v4") == "Q9NZJ5"

    def test_raw_accession(self):
        assert _extract_uniprot_from_target("P12345") == "P12345"

    def test_pdb_id(self):
        # Too short for UniProt, no AF- prefix
        assert _extract_uniprot_from_target("1abc") is None

    def test_empty(self):
        assert _extract_uniprot_from_target("") is None

    def test_none(self):
        assert _extract_uniprot_from_target("") is None

    def test_long_accession(self):
        assert _extract_uniprot_from_target("A0A0A0MRZ7") == "A0A0A0MRZ7"
