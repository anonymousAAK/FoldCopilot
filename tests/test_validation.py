"""Tests for input validation utilities."""

import pytest

from foldcopilot.utils.validation import (
    ValidationError,
    validate_pdb_content,
    validate_sequence,
    validate_sequences,
    validate_uniprot_id,
)


class TestValidateSequence:
    def test_valid_sequence(self):
        result = validate_sequence("MKFLILLFNILCLFPVLAAD")
        assert result == "MKFLILLFNILCLFPVLAAD"

    def test_lowercase_normalized(self):
        result = validate_sequence("mkflillfnilclfpvlaad")
        assert result == "MKFLILLFNILCLFPVLAAD"

    def test_whitespace_stripped(self):
        result = validate_sequence("MKFL ILLF\nNILC LFPV LAAD")
        assert result == "MKFLILLFNILCLFPVLAAD"

    def test_empty_raises(self):
        with pytest.raises(ValidationError, match="empty"):
            validate_sequence("")

    def test_too_short(self):
        with pytest.raises(ValidationError, match="too short"):
            validate_sequence("MKFL")

    def test_too_long(self):
        with pytest.raises(ValidationError, match="too long"):
            validate_sequence("A" * 6000)

    def test_custom_max_length(self):
        with pytest.raises(ValidationError, match="too long"):
            validate_sequence("A" * 100, max_length=50)

    def test_invalid_chars(self):
        with pytest.raises(ValidationError, match="Invalid characters"):
            validate_sequence("MKFL123ABCDE")

    def test_ambiguous_codes_allowed(self):
        result = validate_sequence("MKFLXBZOUJAAA")
        assert "X" in result


class TestValidateSequences:
    def test_valid(self):
        result = validate_sequences(["A" * 20, "G" * 30])
        assert len(result) == 2

    def test_empty_list(self):
        with pytest.raises(ValidationError, match="At least one"):
            validate_sequences([])

    def test_too_many_chains(self):
        with pytest.raises(ValidationError, match="Too many"):
            validate_sequences(["A" * 20] * 15, max_chains=10)


class TestValidateUniprotId:
    def test_standard_accession(self):
        assert validate_uniprot_id("P04637") == "P04637"

    def test_long_accession(self):
        assert validate_uniprot_id("A0A0A0MRZ7") == "A0A0A0MRZ7"

    def test_entry_name(self):
        assert validate_uniprot_id("P53_HUMAN") == "P53_HUMAN"

    def test_lowercase_normalized(self):
        assert validate_uniprot_id("p04637") == "P04637"

    def test_empty(self):
        with pytest.raises(ValidationError, match="empty"):
            validate_uniprot_id("")

    def test_too_short(self):
        with pytest.raises(ValidationError, match="length"):
            validate_uniprot_id("P04")

    def test_special_chars(self):
        with pytest.raises(ValidationError, match="Invalid"):
            validate_uniprot_id("P046@7")


class TestValidatePdbContent:
    def test_valid(self):
        pdb = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 85.50           C\nEND"
        result = validate_pdb_content(pdb)
        assert "ATOM" in result

    def test_empty(self):
        with pytest.raises(ValidationError, match="empty"):
            validate_pdb_content("")

    def test_no_atoms(self):
        with pytest.raises(ValidationError, match="No ATOM"):
            validate_pdb_content("HEADER test\nEND\n")
