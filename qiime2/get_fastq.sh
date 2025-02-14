#!/usr/bin/env bash

#Requires prefecth and fasterq-dump from sra-toolkit
if ! command -v prefetch &> /dev/null; then
    echo "prefetch could not be found. Please install the SRA Toolkit."
    exit 1
fi

if ! command -v fasterq-dump &> /dev/null; then
    echo "fasterq-dump could not be found. Please install the SRA Toolkit."
    exit 1
fi

#Insert CSV file with accession codes
csv_file="sra-explorer.csv"

if [[ ! -f "$csv_file" ]]; then
    echo "The file $csv_file does not exist."
    exit 1
fi

while IFS=',' read -r order url accession; do
    if [[ "$accession" == "accession" ]]; then
        continue
    fi

    #Get accession code
    accession=$(echo "$accession" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//;s/"//g;s/'"'"'//g')
    echo "Processing accession code: ${accession}"

    #Use prefetch to download the SRA file
    if prefetch "$accession"; then
        echo "Downloaded $accession successfully."

        #Convert the downloaded SRA file to FASTQ format
        if fasterq-dump --split-3 "$accession"; then
            echo "Converted $accession to FASTQ format successfully."
        else
            echo "Failed to convert $accession to FASTQ format."
        fi

        #Remove the SRA file after conversion
        rm "${accession}.sra" 2>/dev/null

    else
        echo "Failed to download $accession."
    fi
done < "$csv_file"

echo "All accessions processed."
