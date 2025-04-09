#!/bin/bash

# Check for required arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <file_extension>"
    echo "Example: $0 out_sheep_gut/BIN fa"
    exit 1
fi

input_dir="$1"
file_ext="$2"
output_csv="binning_result.csv"

# Start CSV with header
echo "Bin index,Contig index" > "$output_csv"

# Loop over files with given extension in the input directory
for file in "$input_dir"/*."$file_ext"; do
    bin_name=$(basename "$file" ."$file_ext")
    
    # Extract contig names
    grep "^>" "$file" | sed 's/^>//' | cut -d ' ' -f1 | while read -r contig; do
        echo "${bin_name},${contig}" >> "$output_csv"
    done
done

echo "CSV file generated: $output_csv"