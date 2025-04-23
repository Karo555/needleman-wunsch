from typing import List
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
