import re
from typing import Literal


class Sequence:
    """
    A simple container for an ID and a biological sequence string.
    """

    def __init__(
        self,
        identifier: str,
        sequence: str,
        alphabet: Literal["dna", "protein"] = "dna",
    ):
        self.id = identifier
        self.sequence = sequence.upper()
        self.alphabet = alphabet
        self._validate()

    def _validate(self):
        """
        Ensure the sequence contains only valid characters for the chosen alphabet.
        Empty sequences are now allowed.
        """
        if self.alphabet.lower() == "dna":
            pat = r"^[ACGTacgt]*$"
        elif self.alphabet.lower() == "protein":
            pat = r"^[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]*$"
        else:
            raise ValueError(f"Unknown alphabet: {self.alphabet}")

        if not re.fullmatch(pat, self.sequence):
            raise ValueError(
                f"Invalid characters in sequence {self.id!r} for {self.alphabet}"
            )

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        return f"Sequence(id={self.id!r}, sequence={self.sequence!r})"
