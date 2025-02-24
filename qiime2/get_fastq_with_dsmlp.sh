#!/usr/bin/env bash
SRA_PATH="$HOME/teams/a09/team1/kondo_et_al_2021/sra-toolkit/sratoolkit.3.2.0-ubuntu64/bin"

#Requires prefecth and fasterq-dump from sra-toolkit
if [[ !  -x "$SRA_PATH/prefetch" ]]; then
    echo "prefetch could not be found. Please install the SRA Toolkit."
    exit 1
fi

if [[ ! -x "$SRA_PATH/fasterq-dump" ]]; then
    echo "fasterq-dump could not be found. Please install the SRA Toolkit."
    exit 1
fi

#Insert CSV file with accession codes
csv_file="sra-explorer-all.csv"

if [[ ! -f "$csv_file" ]]; then
    echo "The file $csv_file does not exist."
    exit 1
fi

tail -n +2 "csv_file" | while IFS=',' read -r _ accession; do
    #if [[ "$accession" == "accession" ]]; then
    #    continue
    #fi

    #Get accession code
    accession=$(echo "$accession" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//;s/"//g;s/'"'"'//g')
    echo "Processing accession code: ${accession}"

    #Use prefetch to download the SRA file
    if "$SRA_PATH/prefetch" "$accession"; then
        echo "Downloaded $accession successfully."

        #Convert the downloaded SRA file to FASTQ format
        if "$SRA_PATH/fasterq-dump" --split-files "$accession"; then
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
