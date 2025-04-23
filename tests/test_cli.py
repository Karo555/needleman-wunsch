import pytest
import sys
import json
from pathlib import Path
from src.aligner.cli import parse_args


def test_parse_args_input():
    args = parse_args(
        [
            "--input",
            "file1.fasta",
            "file2.fasta",
            "--match",
            "2",
            "--mismatch",
            "-1",
            "--gap",
            "-3",
        ]
    )
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


def test_cli_matrix_export(tmp_path, monkeypatch):
    # Prep src on path
    project = Path(__file__).parent.parent
    monkeypatch.syspath_prepend(str(project / "src"))

    # Create minimal FASTAs
    data = tmp_path / "data"
    data.mkdir()
    (data / "s1.fasta").write_text(">s1\nA\n")
    (data / "s2.fasta").write_text(">s2\nA\n")

    # Move into tmp dir so relative paths work
    monkeypatch.chdir(tmp_path)

    # Set argv
    out_csv = tmp_path / "matrix_out.csv"
    sys.argv = [
        "aligner.cli",
        "--input",
        "data/s1.fasta",
        "data/s2.fasta",
        "--matrix-out",
        str(out_csv),
    ]

    # Run
    import aligner.cli as cli

    cli.main()

    # Validate the CSV for a simple match/mismatch/gap setup
    text = out_csv.read_text().splitlines()
    # For two single 'A's with default match=1,mismatch=-1,gap=-2, matrix is 2×2:
    # [0, -2]
    # [-2, 1]
    assert text == ["0,-2", "-2,1"]


def test_cli_json_output(tmp_path, monkeypatch):
    # prep src on path
    project = Path(__file__).parent.parent
    monkeypatch.syspath_prepend(str(project / "src"))

    # make FASTAs
    data = tmp_path / "data"
    data.mkdir()
    (data / "s1.fasta").write_text(">s1\nA\n")
    (data / "s2.fasta").write_text(">s2\nA\n")

    # make output folder and cwd
    monkeypatch.chdir(tmp_path)
    reports = tmp_path / "reports"
    reports.mkdir()
    out_json = reports / "out.json"

    # simulate CLI
    sys.argv = [
        "aligner.cli",
        "--input",
        "data/s1.fasta",
        "data/s2.fasta",
        "--json",
        str(out_json),
    ]
    import aligner.cli as cli

    cli.main()

    # assert JSON exists and has expected keys
    assert out_json.exists()
    payload = json.loads(out_json.read_text())
    assert payload["parameters"]["match"] == 1
    assert payload["sequences"]["s1"] == "A"
    assert isinstance(payload["matrix"], list)
    assert isinstance(payload["alignments"], list)


def test_cli_html_output(tmp_path, monkeypatch):
    import sys
    from pathlib import Path

    # Prep src/ on path
    project = Path(__file__).parent.parent
    monkeypatch.syspath_prepend(str(project / "src"))

    # Write minimal FASTAs
    data = tmp_path / "data"
    data.mkdir()
    (data / "s1.fasta").write_text(">s1\nA\n")
    (data / "s2.fasta").write_text(">s2\nA\n")

    # CWD = tmp to use relative paths
    monkeypatch.chdir(tmp_path)

    # Set argv
    out_html = tmp_path / "report.html"
    sys.argv = [
        "aligner.cli",
        "--input",
        "data/s1.fasta",
        "data/s2.fasta",
        "--html",
        str(out_html),
    ]

    import aligner.cli as cli

    cli.main()

    # Check file exists and contains our header
    text = out_html.read_text()
    assert "<h1>Needleman" in text
    assert "Sequence 1: s1: A" in text
