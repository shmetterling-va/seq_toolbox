# Seq toolbox 🧬
Seq toolbox is a small python package with utilities for common bioinformatics files: fastq, fasta and BLAST text output.
Only standart python library is used, so no extra installing is needed.
## Installation 💾
Clone the repository:
```bash
git clone git@github.com:shmetterling-va/seq_toolbox.git
cd seq_toolbox
```
## Usage 🚀

Import functions from main script:

```python
from seq_toolbox import (
    filter_fastq,
    convert_multiline_fasta_to_oneline,
    parse_blast_output,
)
```

### `filter_fastq` 🧹

Works like a strict bouncer for your reads: only good ones get inside 🕶️
It checks GC content, length and mean quality (phred33) of every read.
File is read one read at a time, so even a huge file will not eat all your RAM.
Good reads are saved to `filtered/` folder (it is created automaticly if it not exist).

```python
filter_fastq(
    "reads.fastq",
    "good_reads.fastq",
    gc_bounds=(20, 80),      # GC from 20 to 80 %
    length_bounds=(50, 150), # length from 50 to 150
    quality_threshold=30,    # mean quality 30 or more
)
```

If only one number is given to `gc_bounds` or `length_bounds` it means upper bound. For example `gc_bounds=44.4` keeps reads with GC up to 44.4 %.
All bounds are included ✅

### `convert_multiline_fasta_to_oneline` 📏

Sometimes fasta files split one sequence into many lines and nobody knows why 🤷‍♀️
This function puts every sequence back into one line.
If `output_fasta` is not given, file `<input name>_oneline.fasta` is created next to input file.

```python
convert_multiline_fasta_to_oneline("genes.fasta")
convert_multiline_fasta_to_oneline("genes.fasta", "genes_oneline.fasta")
```

### `parse_blast_output` 🔍

Takes BLAST results downloaded as txt file (the one that looks like the matrix 😵‍💫).
For every query it takes the best hit (first line in `Description` column), removes duplicates and saves all proteins sorted by alphabet, one per line.

```python
parse_blast_output("blast_results.txt", "best_hits.txt")
```

## Structure 📂

```
seq_toolbox/
|- README.md
|- seq_toolbox.py      # main script with all utilities
|- helpers/
    |- fastq.py        # helper functions for filter_fastq
    |- bio_files.py    # helper functions for fasta and BLAST
```
