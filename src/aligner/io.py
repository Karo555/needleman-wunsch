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
