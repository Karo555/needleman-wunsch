import pytest
from src.aligner.io import read_fasta
from src.aligner.models import Sequence


def test_read_fasta_single(tmp_path):
    content = ">seq1\nACGTACGT\n"
    path = tmp_path / "single.fasta"
    path.write_text(content)
    records = read_fasta(str(path))
    assert len(records) == 1
    seq = records[0]
    assert isinstance(seq, Sequence)
    assert seq.id == "seq1"
    assert seq.sequence == "ACGTACGT"


def test_read_fasta_multiple(tmp_path):
    content = ">s1\nAAA\n>s2\nTTT\n"
    path = tmp_path / "multi.fasta"
    path.write_text(content)
    records = read_fasta(str(path))
    assert [rec.id for rec in records] == ["s1", "s2"]
    assert [rec.sequence for rec in records] == ["AAA", "TTT"]


def test_read_fasta_missing_file():
    with pytest.raises(FileNotFoundError):
        read_fasta("nonexistent.fasta")


def test_read_fasta_malformed(tmp_path):
    # Sequence data before any header
    path = tmp_path / "bad.fasta"
    path.write_text("ACGTACGT\n")
    with pytest.raises(ValueError):
        read_fasta(str(path))


def test_read_fasta_empty(tmp_path):
    path = tmp_path / "empty.fasta"
    path.write_text("")
    with pytest.raises(ValueError):
        read_fasta(str(path))


def test_read_fasta_protein(tmp_path):
    content = ">prot1\nARNDCQ\n"
    path = tmp_path / "protein.fasta"
    path.write_text(content)
    records = read_fasta(str(path), alphabet="protein")
    seq = records[0]
    assert seq.id == "prot1"
    assert seq.sequence == "ARNDCQ"
