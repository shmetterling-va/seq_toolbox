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


def read_best_blast_hits(input_file: str) -> Iterator[str]:
    """
    Reads a BLAST text output and finds the best hit for every query.

    Arguments:
    input_file: path to the BLAST results txt file

    Yields the Description of the first hit in each
    "Sequences producing significant alignments" table.
    """
    description_width = None
    with open(input_file) as blast_file:
        for line in blast_file:
            if line.startswith("Description"):
                # the Description column ends where the next column ("Name") starts
                description_width = line.index("Name")
            elif description_width is not None:
                yield line[:description_width].strip()
                description_width = None


def write_lines(output_file: str, lines: list[str]) -> None:
    """
    Writes strings to a file, one per line.

    Arguments:
    output_file: path to the output file
    lines: strings to write
    """
    with open(output_file, "w") as file:
        for line in lines:
            file.write(line + "\n")
