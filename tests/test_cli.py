import pytest
from src.aligner.cli import parse_args


def test_parse_args_input():
    args = parse_args([
        "--input", "file1.fasta", "file2.fasta",
        "--match", "2", "--mismatch", "-1", "--gap", "-3"
    ])
    assert args.input == ["file1.fasta", "file2.fasta"]
    assert not args.manual
    assert args.match == 2
    assert args.mismatch == -1
    assert args.gap == -3
    assert args.output is None
    assert args.plot is None


def test_parse_args_manual():
    args = parse_args(["--manual"])
    assert args.manual
    assert args.input is None
    assert args.output is None
    assert args.plot is None


def test_help_shows_usage(capsys):
    with pytest.raises(SystemExit) as excinfo:
        parse_args(["--help"])
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "usage:" in captured.out.lower()

