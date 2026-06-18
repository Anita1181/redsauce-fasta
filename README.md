# RedSauce FASTA - Sequence Analyser for .fasta/.fa Files

It reports basic sequence statistics, finds open reading frames (ORFs), and detects repeated subsequences — built and tested against the *Enterobacteria phage lambda* genome.

## Features

- **Record counting** (counts the number of sequences in the file)
- **Length analysis** (identifies the longest/shortest sequences in the file)
- **ORF detection** (scans a chosen reading frame for open-reading-frames (start codon `ATG` to the nearest in-frame stop codon `TAA`/`TAG`/`TGA`) and reports the longest ORF per sequence and overall)
- **Detects repeat seqs** (finds all repeated subsequences of k-mers and reports the most frequent one)

## Requirements

- Python 3.7+
- [Biopython](https://biopython.org/)

Install the dependency with:

```bash
pip install biopython
```

## Usage

Place a FASTA file (e.g. used in the .py file `lambda_virus.fa`) in the same directory as the script, then run:

```bash
python3 fasta_analyzer.py
```

By default the script runs against `lambda_virus.fa` and:
1. Prints the number of sequence records in the file
2. Prints the longest and shortest sequence length(s), with IDs
3. Finds and prints the longest ORF in reading frame 1
4. Finds and prints all repeated 15-base subsequences

To analyze a different file or reading frame, edit the variables at the bottom of the script:

```python
if __name__ == "__main__":
    file = "lambda_virus.fa" # add file name here

    print(record_counter(file))
    analyze_fasta_lengths(file)
    analyze_orfs(file, 1)      # insert 1/2/3 in second param
    find_repeats(file, 15)     # insert chosen count in second param
```

## Example Output

Running the script against the included lambda phage genome (a single 48,502 bp sequence) produces output like:

*Longest seq length: 48502*
*Number of longest seqs: 1*
  *ID: gi|9626243|ref|NC_001416.1|*

*Shortest seq length: 48502*
*Number of shortest seqs: 1*
  *ID: gi|9626243|ref|NC_001416.1|*

*Reading Frame: 1
  Longest ORF in file: 3399 bp*
    *Sequence ID: gi|9626243|ref|NC_001416.1|
  Start position: 15505
  ORF: ATGGGTAAAGGAAGCAGT...*

*Repeats of length 15:
  CATGACGGAGGATGA: 2
Most frequent: CATGACGGAGGATGA (2 times)*

**Note:** the lambda genome FASTA file contains only one sequence record, so the longest and shortest lengths being identical is expected - not a bug.
