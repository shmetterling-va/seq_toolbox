"""Helper functions for filter_fastq: reading, checking and writing reads."""

import os
from typing import Iterator, TextIO

Bounds = tuple[float, float]
FastqRecord = tuple[str, str, str, str]


def normalize_bounds(bounds: float | Bounds) -> Bounds:
    """
    Converts filter bounds to a (lower, upper) tuple.

    Arguments:
    bounds: a single number (treated as the upper bound) or a (lower, upper) tuple

    Returns a tuple (lower, upper).
    """
    if isinstance(bounds, (int, float)):
        return (0, bounds)
    lower, upper = bounds
    return (lower, upper)


def read_fastq(input_fastq: str) -> Iterator[FastqRecord]:
    """
    Reads a fastq file one record at a time.

    Arguments:
    input_fastq: path to the input fastq file

    Yields tuples (name, sequence, comment, quality) without newline characters.
    """
    with open(input_fastq) as fastq_file:
        while True:
            name = fastq_file.readline().rstrip("\n")
            if not name:
                break
            sequence = fastq_file.readline().rstrip("\n")
            comment = fastq_file.readline().rstrip("\n")
            quality = fastq_file.readline().rstrip("\n")
            yield name, sequence, comment, quality


def write_record(output_file: TextIO, record: FastqRecord) -> None:
    """
    Writes one fastq record (4 lines) to an already opened file.

    Arguments:
    output_file: file object opened for writing
    record: tuple (name, sequence, comment, quality)
    """
    output_file.write("\n".join(record) + "\n")


def make_output_path(output_fastq: str, folder: str = "filtered") -> str:
    """
    Creates the output folder if needed and returns the path to the output file.

    Arguments:
    output_fastq: name of the output fastq file
    folder: folder to save the file to, "filtered" by default

    Returns the path folder/output_fastq.
    """
    os.makedirs(folder, exist_ok=True)
    return os.path.join(folder, output_fastq)


def gc_content(sequence: str) -> float:
    """
    Computes GC content of a sequence in percent.

    Arguments:
    sequence: nucleotide sequence

    Returns GC content from 0 to 100 (0 for an empty sequence).
    """
    if not sequence:
        return 0
    gc_count = sequence.upper().count("G") + sequence.upper().count("C")
    return gc_count / len(sequence) * 100


def mean_quality(quality: str) -> float:
    """
    Computes the mean phred33 quality of a read.

    Arguments:
    quality: quality string in phred33 encoding

    Returns the mean quality score (0 for an empty string).
    """
    if not quality:
        return 0
    return sum(ord(symbol) - 33 for symbol in quality) / len(quality)


def is_in_bounds(value: float, bounds: Bounds) -> bool:
    """
    Checks whether a value lies within bounds, both bounds included.

    Arguments:
    value: number to check
    bounds: tuple (lower, upper)

    Returns True if lower <= value <= upper.
    """
    lower, upper = bounds
    return lower <= value <= upper


def passes_gc(sequence: str, gc_bounds: Bounds) -> bool:
    """Returns True if the GC content of the sequence is within gc_bounds."""
    return is_in_bounds(gc_content(sequence), gc_bounds)


def passes_length(sequence: str, length_bounds: Bounds) -> bool:
    """Returns True if the length of the sequence is within length_bounds."""
    return is_in_bounds(len(sequence), length_bounds)


def passes_quality(quality: str, quality_threshold: float) -> bool:
    """Returns True if the mean read quality is not below quality_threshold."""
    return mean_quality(quality) >= quality_threshold
