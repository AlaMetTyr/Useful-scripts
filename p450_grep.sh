#!/bin/bash -e


input_dir=/nesi/nobackup/ga03488/Amy/genome_datasets/2023_datasetfaa
output_dir=/nesi/nobackup/ga03488/Amy/genome_datasets/2023_datasetfaa/subset_p450_faa

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
  # Find the lines containing the P450 annotation
  grep -A1 "P450" "${file}" | grep -v "^--$" > "${output_dir}/${filename}_p450.faa"
done
