from Bio import SeqIO
from collections import Counter

#find the amount of records in a file - works
def record_counter(file_path):
    count = 0
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip()[0:1]=='>':
                count+=1
    return count

#find max, min lengths of sequences - works
def analyze_fasta_lengths(file_path):
    sequences = list(SeqIO.parse(file_path, "fasta"))

    if not sequences:
        return "No sequences found."

    #lengths
    lengths = [len(seq.seq) for seq in sequences]
    max_len = max(lengths)
    min_len = min(lengths)

    #max and min
    longest_seqs = [seq for seq in sequences if len(seq.seq) == max_len]
    shortest_seqs = [seq for seq in sequences if len(seq.seq) == min_len]

    print(f"Longest seq length: {max_len}")
    print(f"Number of longest seqs: {len(longest_seqs)}")

    for seq in longest_seqs:
        print(f" - ID: {seq.id}")

    print(f"\nShortest seq length: {min_len}")
    print(f"Number of shortest seqs: {len(shortest_seqs)}")
    for seq in shortest_seqs:
        print(f" - ID: {seq.id}")

#find orfs in reading frame - works
def find_orfs(seq, frame):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []

    start_pos = frame-1
    seq = seq[start_pos:]
    i=0
    while i < len(seq) - 2:
        codon = seq[i:i+3]
        if codon == start_codon:
            # potential ORF start?
            for j in range(i + 3, len(seq) - 2, 3):
                stop_codon = seq[j:j+3]
                if stop_codon in stop_codons:
                    orf_seq = seq[i:j+3]
                    orfs.append((start_pos + i + 1, orf_seq))
                    break
        i += 3
    return orfs


#analyse orfs in reading frame - works
def analyze_orfs(file_path, frame):
    if frame not in {1, 2, 3}:
        raise ValueError("Frame must be 1, 2, or 3.")
    all_orfs = []
    longest_orf_overall = ("", 0, 0, "")  # (seq_id, length, position, orf_seq)
    longest_orf_per_seq = {}

    for record in SeqIO.parse(file_path, "fasta"):
        seq_id = record.id
        seq_str = str(record.seq).upper()
        orfs = find_orfs(seq_str, frame)
        
        if orfs:
                # sort ORFs by length - tested, works
                orfs.sort(key=lambda x: len(x[1]), reverse=True)
                start_pos, orf_seq = orfs[0]
                longest_orf_per_seq[seq_id] = {
                    "length": len(orf_seq),
                    "start_position": start_pos,
                    "orf_sequence": orf_seq
                }

                if len(orf_seq) > longest_orf_overall[1]:
                    longest_orf_overall = (seq_id, len(orf_seq), start_pos, orf_seq)
                    
    print(f"Reading Frame: {frame}")
    if longest_orf_overall[1] > 0:
        print(f"\nLongest ORF in file: {longest_orf_overall[1]} bp")
        print(f" - Sequence ID: {longest_orf_overall[0]}")
        print(f" - Start position: {longest_orf_overall[2]}")
        print(f" - ORF: {longest_orf_overall[3]}")
    else:
        print("No valid ORFs found.")

    return longest_orf_per_seq

#find dna bp repeats with given length - works
def find_repeats(file_path, n):
    counts = Counter()
    for record in SeqIO.parse(file_path, "fasta"):
        seq = str(record.seq).upper()
        for i in range(len(seq) - n + 1):
            counts[seq[i:i+n]] += 1

    repeats = {k: v for k, v in counts.items() if v > 1}
    if not repeats:
        print(f"No repeats of length {n}.")
        return

    most = max(repeats.values())
    most_freq = [k for k, v in repeats.items() if v == most]

    print(f"Repeats of length {n}:")
    for k, v in repeats.items():
        print(f"{k}: {v}")
    print(f"\nMost frequent: {', '.join(most_freq)} ({most} times)")


if __name__ == "__main__":
    file = "lambda_virus.fa"

    print(record_counter(file))
    analyze_fasta_lengths(file)
    analyze_orfs(file, 1)
    find_repeats(file, 15)