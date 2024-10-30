#!/bin/bash -e

#check whether the zipped file is valid or not
for file in /path/to/your/.directory/*.g.vcf.gz; do
    gunzip -t "$file" || echo "Corrupted file: $file"
done
