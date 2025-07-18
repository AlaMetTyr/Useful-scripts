from Bio import SeqIO
from Bio.Seq import Seq

IUPAC_CODES = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"}, "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"}, "N": {"A", "C", "G", "T"}
}

import sys

fasta_file = sys.argv[1]
forward_primer = sys.argv[2].upper()
reverse_primer = sys.argv[3].upper()
output_file = sys.argv[4]
three_prime_length = 5

rev_primer_rc = str(Seq(reverse_primer).reverse_complement())

def count_mismatches(seq_segment, primer_segment):
    return sum(
        seq_base not in IUPAC_CODES[primer_base]
        for seq_base, primer_base in zip(seq_segment, primer_segment)
    )

def get_best_match(seq, primer):
    best_idx = None
    least_mismatches = len(primer) + 1
    for i in range(len(seq) - len(primer) + 1):
        segment = seq[i:i+len(primer)]
        mismatches = count_mismatches(segment, primer)
        if mismatches < least_mismatches:
            least_mismatches = mismatches
            best_idx = i
    return best_idx, least_mismatches

with open(fasta_file, "r") as infile, open(output_file, "w") as outfile:
    kept, discarded = 0, 0

    for record in SeqIO.parse(infile, "fasta"):
        seq = str(record.seq).upper()

        fwd_idx, _ = get_best_match(seq, forward_primer)
        rev_idx, _ = get_best_match(seq, rev_primer_rc)

        if fwd_idx is None or rev_idx is None or fwd_idx >= rev_idx:
            discarded += 1
            continue

        # Check 3' ends
        fwd_seq_end = seq[fwd_idx + len(forward_primer) - three_prime_length : fwd_idx + len(forward_primer)]
        fwd_primer_end = forward_primer[-three_prime_length:]
        fwd_mismatches = count_mismatches(fwd_seq_end, fwd_primer_end)

        rev_seq_start = seq[rev_idx : rev_idx + three_prime_length]
        rev_primer_start = rev_primer_rc[:three_prime_length]
        rev_mismatches = count_mismatches(rev_seq_start, rev_primer_start)

        if fwd_mismatches == 0 and rev_mismatches == 0:
            outfile.write(f">{record.description}\n{seq}\n")
            kept += 1
        else:
            discarded += 1

print(f"Kept: {kept} sequences")
print(f"Discarded: {discarded} sequences (3' mismatches)")
