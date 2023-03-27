#!/bin/bash -e

#SBATCH --job-name=Interproscan # job name (shows up in the queue)
#SBATCH --time=24:01:00      # Walltime (HH:MM:SS)
#SBATCH --mem=70GB          # Memory in MB
#SBATCH --account=ga03488
#SBATCH --cpus-per-task=4  

module purge
module load InterProScan/5.51-85.0-gimkl-2020a-Perl-5.30.1-Python-3.8.2

for f in *.fa
do
  interproscan.sh -i ./${f%*} -dp -iprlookup --goterms --appl Pfam,PRINTS,PANTHER
done

