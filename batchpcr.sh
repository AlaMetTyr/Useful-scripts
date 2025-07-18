#!/bin/bash -e
#SBATCH --job-name=PCR
#SBATCH --time=0:10:00
#SBATCH --mem=4GB
#SBATCH --output=pcr.out
#SBATCH --error=pcr.err

module purge
module load Python/3.11.6-foss-2023a

shopt -s nullglob

# The forward and reverse primer we wish to do our in silico PCR with
FORWARD_PRIMER={"insert primer seq here"}
REVERSE_PRIMER={"insert primer seq here"}

# Loop through all FASTA files in the directory
for file in ./high_barcode_frew/*.{fas,fasta}; do
    filename=$(basename "$file")  # just the filename
    species_name=$(echo "$filename" | awk -F'_' '{print $1 "_" $2}')

    # Output file name and location
    output_file="./insilico_pcr/${species_name}_pcrresult.txt"

    # Run the in silico PCR script
    python in_silico_pcr.py "$file" "$FORWARD_PRIMER" "$REVERSE_PRIMER" "$output_file"

    echo "Processed $file -> $output_file"
done

echo "PCR complete! Prase be the PCR gods"
