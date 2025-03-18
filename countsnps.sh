#!/bin/bash -e
bcftools view -v snps -H {filename}.vcf | wc -l
