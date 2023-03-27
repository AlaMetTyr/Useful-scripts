#!/bin/bash -e

#SBATCH --job-name=Interproscan # job name (shows up in the queue)
#SBATCH --time=24:01:00      # Walltime (HH:MM:SS)
#SBATCH --mem=70GB          # Memory in MB
#SBATCH --account=ga03488
#SBATCH --cpus-per-task=4  cafe5

for f in *.tsv.txt; do
base_name=${f%_Orthogroups.GeneCount.tsv.txt}
ID=${base_name#*_}
/nesi/project/ga03488/software/CAFE5/bin/cafe5 -i "${f%}" -t "${ID}_SpeciesTree_rooted.tree.ultrametric.tre" > "output_${ID}.out"
done
