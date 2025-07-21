# **eDNA surveillance barcodes project**  
`in_silico_pcr.py` and `check_mismatch.py` were both used to check COI barcodes downloaded from available repositories against universal barcodes to identify unique barcodes for an alert watchlist.  
## 1. Run PCR
`in_silico_pcr.py` runs an in silico PCR and records the numbers of mismatches ot the f&r primer of each 'reaction. This was submitted as a batch job using `batchpcr.sh`, which also defines the primers.  
## 2. Check mismatches
`check_mismtach.py` looked at any mismatches on the primers and removed sequences that had mismatches within the first 5bp of the 3' end of each primer.  
## 4. Remove sequence qwith no amplicon  
`remove_seq` is an awk one liner to remove any of the sequences that are in the resulting .txt file that say 'No amplicon found'
## 4. Random scripts  
`rename_bulk_barcodes` is an awk oneliner to change the naming conformity of any barcodes compiled together in one fasta alignment file.  
## 5. Random scripts  	
`fuzznuc` was also used for forward and reverse primers to find the sequences and mismatches individually (used before python script was written).  



# **Other scripts for multiple projects to make life easier**  
## 1. Count SNPs
`count_snps.sh` uses bcftools to count the number of single nucleotide polymorphisms in a given vcf file.
## 2. Zip file corrupt?
`corrupt_check.sh` will check if a zipped file is corrupt (used this when trying to acquire data from overseas collaborators).   
## 3. BLAST  
Speaks for itself. `blast.sl` is a slurm script to blastp a reference fasta you provide   
## 3. P450s  
I wanted to find annotated P450s in sequences, so used `p450_grep.sh` and  `P450_blast.sh` to find them (though now I dont remember why).  
## 4. Cafe5  
Script that runs cafe5 is `cafe5.sl`. Needs to be ran afetr orthofinder `orthofinder_faa.sl`. More instructinos in the invasion genomics repository readme.
