from typing import ClassVar, Set


class Sequence:
    """
    Represents a biological sequence with validation against DNA or protein alphabets.
    """

    dna_alphabet: ClassVar[Set[str]] = set("ACGT")
    protein_alphabet: ClassVar[Set[str]] = set("ARNDCQEGHILKMFPSTWYV")

    def __init__(self, identifier: str, sequence: str, alphabet: str = "dna"):
        """
        :param identifier: a unique name or header for the sequence
        :param sequence: the raw sequence string (will be uppercased)
        :param alphabet: "dna" or "protein"
        """
        self.id = identifier
        self.sequence = sequence.strip().upper()

        if alphabet == "dna":
            valid_chars = self.dna_alphabet
        elif alphabet == "protein":
            valid_chars = self.protein_alphabet
        else:
            raise ValueError(f"Unknown alphabet: '{alphabet}'")

        invalid = set(self.sequence) - valid_chars
        if invalid:
            raise ValueError(f"Invalid characters for {alphabet} alphabet: {invalid}")

    def __len__(self) -> int:
        return len(self.sequence)

    def __repr__(self) -> str:
        return f"Sequence(id={self.id!r}, length={len(self)})"
