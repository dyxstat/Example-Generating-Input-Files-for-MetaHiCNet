#!/bin/bash

# Usage: ./extract_mge.sh mge_results.csv assembly.fa cutoff

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <mge_results.csv> <assembly.fa> <cutoff>"
    exit 1
fi

csv_file="$1"
assembly_file="$2"
cutoff="$3"

# Temporary files
phage_ids="phage_ids.txt"
plasmid_ids="plasmid_ids.txt"

# Extract contigs with phage_score or plasmid_score >= cutoff
awk -F',' -v c="$cutoff" 'NR>1 && $3+0 >= c {print $1}' "$csv_file" > "$phage_ids"
awk -F',' -v c="$cutoff" 'NR>1 && $5+0 >= c {print $1}' "$csv_file" > "$plasmid_ids"

# Function to extract FASTA entries for given IDs
extract_contigs () {
    id_file="$1"
    input_fasta="$2"
    output_fasta="$3"

    awk 'BEGIN {
             while ((getline < "'$id_file'") > 0) ids[$1]=1
         }
         /^>/ {
             h=substr($0,2); keep=ids[h]
         }
         keep' "$input_fasta" > "$output_fasta"
}

# Extract matching contigs
extract_contigs "$phage_ids" "$assembly_file" "viral_contig.fa"
extract_contigs "$plasmid_ids" "$assembly_file" "plasmid_contig.fa"

# Clean up
rm "$phage_ids" "$plasmid_ids"

echo "Done. Output: viral_contig.fa and plasmid_contig.fa"