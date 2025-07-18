import sys
from Bio import SeqIO
from Bio.Seq import Seq

# IUPAC nucleotide codes for degeneracy
IUPAC_CODES = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"}, "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"}, "N": {"A", "C", "G", "T"}
}

def matches_with_mismatches(seq, primer, max_mismatches=4):
    """Finds the best match for a primer within a sequence, allowing up to max_mismatches mismatches.
       Returns the match position and number of mismatches.
    """
    primer_length = len(primer)
    best_match = None
    best_mismatch_count = max_mismatches + 1  # Start higher than allowed mismatches

    for i in range(len(seq) - primer_length + 1):
        mismatches = sum(
            seq[i + j] not in IUPAC_CODES[primer[j]]
            for j in range(primer_length)
        )

        if mismatches <= max_mismatches and mismatches < best_mismatch_count:
            best_match = i
            best_mismatch_count = mismatches

    if best_match is not None:
        return best_match, best_mismatch_count
    else:
        return None, None

def find_amplicon(seq, fwd_primer, rev_primer, max_mismatches=4):
    """Finds amplicons allowing up to max_mismatches mismatches per primer.
       Returns the amplicon and mismatch counts for forward and reverse primers.
    """
    rev_primer_rc = str(Seq(rev_primer).reverse_complement())

    fwd_match, fwd_mismatches = matches_with_mismatches(seq, fwd_primer, max_mismatches)
    rev_match, rev_mismatches = matches_with_mismatches(seq, rev_primer_rc, max_mismatches)

    if fwd_match is not None and rev_match is not None and fwd_match < rev_match:
        amplicon = seq[fwd_match:rev_match + len(rev_primer_rc)]
        return amplicon, fwd_mismatches, rev_mismatches

    return None, None, None

# --- MAIN SCRIPT ---

# Read input arguments
fasta_file = sys.argv[1]
forward_primer = sys.argv[2]
reverse_primer = sys.argv[3]
output_file = sys.argv[4]

# Initialize mismatch summary counter
summary_counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}

# Process FASTA file
with open(fasta_file, "r") as fasta, open(output_file, "w") as output:
    for record in SeqIO.parse(fasta, "fasta"):
        sequence = str(record.seq)
        amplicon, fwd_mm, rev_mm = find_amplicon(sequence, forward_primer, reverse_primer)

        if amplicon:
            output.write(f">{record.id} | FWD mismatches: {fwd_mm}, REV mismatches: {rev_mm}\n{amplicon}\n")
            total_mm = fwd_mm + rev_mm
            if total_mm <= 4:
                summary_counts[total_mm] += 1
        else:
            output.write(f">{record.id}\nNo amplicon found\n")

    # Write summary at the end
    output.write("\n# Summary of mismatch counts:\n")
    for mm in range(5):
        output.write(f"# Total amplicons with {mm} mismatches: {summary_counts[mm]}\n")

print(f"Processing complete! Results saved to {output_file}")

