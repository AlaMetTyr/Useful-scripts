#!/bin/bash -e

#SBATCH --job-name      BLAST
#SBATCH --time          02:30:00
#SBATCH --mem           150G  # 30 GB plus the database
#SBATCH --ntasks        1
#SBATCH --cpus-per-task 36    # half a node
#SBATCH --account=ga03488

module load BLAST/2.13.0-GCC-11.3.0 
module load BLASTDB/2023-01

# This script takes one argument, the FASTA file of query sequences.
QUERIES=$1
FORMAT="6 qseqid qstart qend qseq sseqid sgi sacc sstart send staxids sscinames stitle length evalue bitscore"
BLASTOPTS="-task blastp"
BLASTAPP=blastp
DB=refseq_protein

# Keep the database in RAM
cp $BLASTDB/{$DB,taxdb}* $TMPDIR/ 
export BLASTDB=$TMPDIR

$BLASTAPP $BLASTOPTS -db $DB -query $QUERIES -outfmt "$FORMAT" \
    -out $QUERIES.$DB.$BLASTAPP -num_threads $SLURM_CPUS_PER_TASK
 