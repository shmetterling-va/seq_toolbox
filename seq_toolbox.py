"""Seq toolbox: utilities for fastq, fasta and BLAST output files."""

from helpers.fastq import (
    make_output_path,
    normalize_bounds,
    passes_gc,
    passes_length,
    passes_quality,
    read_fastq,
    write_record,
)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: float | tuple[float, float] = (0, 100),
    length_bounds: float | tuple[float, float] = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Filters reads from a fastq file by GC content, length and mean quality.
    Reads are processed one at a time, so the whole file is never kept in memory.

    Arguments:
    input_fastq: path to the input fastq file
    output_fastq: name of the output file, it is saved to the "filtered" folder
    gc_bounds: GC content interval in percent, (0, 100) by default;
        a single number is treated as the upper bound
    length_bounds: read length interval, (0, 2**32) by default;
        a single number is treated as the upper bound
    quality_threshold: minimal mean read quality (phred33), 0 by default

    All bounds are inclusive. Returns None, the result is written to a file.
    """
    gc_bounds = normalize_bounds(gc_bounds)
    length_bounds = normalize_bounds(length_bounds)
    output_path = make_output_path(output_fastq)

    with open(output_path, "w") as output_file:
        for record in read_fastq(input_fastq):
            _, sequence, _, quality = record
            if (
                passes_gc(sequence, gc_bounds)
                and passes_length(sequence, length_bounds)
                and passes_quality(quality, quality_threshold)
            ):
                write_record(output_file, record)
