`in_silico_pcr.py` and `check_mismatch.py` were both used to check COI barcodes downloaded from available repositories against universal barcodes to identify unique barcodes for an alert watchlist.  
`in_silico_pcr.py` runs an in silico PCR and records the numbers of mismatches ot the f&r primer of each 'reaction. This was submitted as a batch job using `batchpcr.sh`, which also defines the primers.  
`check_mismtach.py` looked at any mismatches on the primers and removed sequences that had mismatches within the first 5bp of the 3' end of each primer.  

