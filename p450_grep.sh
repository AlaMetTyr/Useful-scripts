#!/bin/bash -e


input_dir=/nesi/nobackup/ga03488/Amy/genome_datasets/2023_datasetfaa
output_dir=/nesi/nobackup/ga03488/Amy/genome_datasets/2023_datasetfaa/full_subset_p450_faa

# Check if input directory exists
if [ ! -d "${input_dir}" ]; then
  echo "Input directory ${input_dir} not found"
  exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "${output_dir}"

for file in "${input_dir}"/*.faa
do
  # Get the filename without the extension
  filename=$(basename "${file}" .faa)
  # Find the lines containing the P450 annotation and print the subsequent sequence lines until the next header
  perl -ne 'if (/^>/) { if($seq){print "$seq\n";} $seq=""; print $_; } else { s/\s+//g; $seq .= $_; } END { print "$seq\n"; }' "${file}" | awk '/P450/{flag=1;print;next}/^>/{flag=0}flag' > "${output_dir}/${filename}_p450.faa"
done