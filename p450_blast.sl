#!/bin/bash -e

#SBATCH --job-name      BLAST
#SBATCH --time          08:30:00
#SBATCH --mem           30G  
#SBATCH --ntasks        1
#SBATCH --account=ga03488

module load BLAST/2.13.0-GCC-11.3.0 
module load BLASTDB/2023-01

# Set the name of the directory containing the query files
query_dir="../."

# Loop over all files in the query directory that end in .faa
for query_file in "$query_dir"/*.faa; do
    # Get the filename without the extension
    filename=$(basename "$query_file" .faa)
    # Run the blastp command with the appropriate arguments
    blastp -query "$query_file" -db clan4_db -num_threads 8 -evalue 1.0e-6 -out "clan4$filename".out -outfmt 6
done
