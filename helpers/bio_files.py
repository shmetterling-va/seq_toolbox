"""Helper functions for working with fasta files and BLAST text output."""

import os
from typing import Iterator


def make_oneline_fasta_name(input_fasta: str) -> str:
    """
    Builds a default name for the oneline fasta file.

    Arguments:
    input_fasta: path to the input fasta file

    Returns a path like "<input name>_oneline.fasta" next to the input file.
    """
    base, _ = os.path.splitext(input_fasta)
    return base + "_oneline.fasta"


def read_fasta(input_fasta: str) -> Iterator[tuple[str, str]]:
    """
    Reads a fasta file where sequences may span several lines.

    Arguments:
    input_fasta: path to the input fasta file

    Yields tuples (header, sequence) one record at a time,
    the sequence is joined into a single line.
    """
    header = None
    sequence_parts: list[str] = []
    with open(input_fasta) as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)
                header = line
                sequence_parts = []
            else:
                sequence_parts.append(line)
    if header is not None:
        yield header, "".join(sequence_parts)
