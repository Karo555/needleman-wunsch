from typing import List, Tuple
from src.aligner.models import Sequence


def read_fasta(path: str, alphabet: str = "dna") -> List[Sequence]:
    """
    Read a FASTA file and return a list of Sequence objects.
    :param path: Path to the FASTA file
    :param alphabet: "dna" or "protein"
    :raises FileNotFoundError: if the file does not exist
    :raises ValueError: on malformed FASTA or no records
    """
    records: List[Sequence] = []
    seq_id = None
    seq_lines: List[str] = []

    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save previous record
                if seq_id is not None:
                    seq_str = "".join(seq_lines)
                    records.append(Sequence(seq_id, seq_str, alphabet))
                seq_id = line[1:].strip()
                seq_lines = []
            else:
                if seq_id is None:
                    raise ValueError("FASTA format error: sequence data before header")
                seq_lines.append(line)
        # End of file: flush last record
        if seq_id is not None:
            seq_str = "".join(seq_lines)
            records.append(Sequence(seq_id, seq_str, alphabet))

    if not records:
        raise ValueError("No FASTA records found in file")

    return records


def read_manual(alphabet: str = "dna") -> Tuple[Sequence, Sequence]:
    """
    Prompt the user to manually enter two sequences and return as Sequence objects.
    :param alphabet: "dna" or "protein"
    """
    print("Enter two sequences for alignment (alphabet: {}):".format(alphabet))
    id1 = input("Identifier for sequence 1: ").strip()
    seq1 = input("Sequence 1: ").strip()
    id2 = input("Identifier for sequence 2: ").strip()
    seq2 = input("Sequence 2: ").strip()

    # Validate using Sequence class
    s1 = Sequence(id1, seq1, alphabet)
    s2 = Sequence(id2, seq2, alphabet)
    return s1, s2


def write_report(path: str, report: str) -> None:
    """
    Write the alignment report to a text file.
    """
    with open(path, "w") as f:
        f.write(report)


def format_report(
    seq1: Sequence,
    seq2: Sequence,
    aligned1: str,
    aligned2: str,
    match: int,
    mismatch: int,
    gap: int,
) -> str:
    """
    Format alignment parameters, sequences, alignment, and metrics into a report string.
    """
    length = len(aligned1)
    identical = sum(1 for a, b in zip(aligned1, aligned2) if a == b)
    percent = (identical / length * 100) if length > 0 else 0.0
    total_gaps = aligned1.count("-") + aligned2.count("-")

    lines = [
        "Parameters:",
        f"  Match score: {match}",
        f"  Mismatch score: {mismatch}",
        f"  Gap penalty: {gap}",
        "",
        "Sequences:",
        f"  {seq1.id}: {seq1.sequence}",
        f"  {seq2.id}: {seq2.sequence}",
        "",
        "Alignment:",
        f"  {aligned1}",
        f"  {aligned2}",
        "",
        "Statistics:",
        f"  Alignment length: {length}",
        f"  Identical positions: {identical}",
        f"  Percentage identity: {percent:.2f}%",
        f"  Total gaps: {total_gaps}",
    ]
    return "\n".join(lines)
