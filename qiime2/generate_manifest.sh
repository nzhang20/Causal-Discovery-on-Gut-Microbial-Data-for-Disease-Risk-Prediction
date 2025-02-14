#!/bin/bash

manifest_file="manifest.tsv"

#header
echo -e "sample-id\tabsolute-filepath" > "$manifest_file"

#Keeping forward reads only; end in _1.fastq
find "$PWD" -type f -name "*_1.fastq" | while read -r filepath; do
    #sample ID (remove path and _1.fastq suffix)
    sample_id=$(basename "$filepath" | sed 's/_1.fastq//')

    echo -e "${sample_id}\t${filepath}" >> "$manifest_file"
done

echo "Manifest file created: $manifest_file"

