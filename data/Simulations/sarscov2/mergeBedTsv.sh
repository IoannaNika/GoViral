#!/bin/bash

# Example input files
fileBed="nCoV-2019.bed"
fileTsv="nCoV-2019.tsv"

output_file="merged.bed"

> "$output_file"

while IFS=$'\t' read -r seqId start end info ph strand ; do

    # Extract the number (the part between the last underscores)
    numberBed=$(echo "$info" | sed 's/.*_\([0-9]*\)_[A-Za-z]*$/\1/')
    echo $info
    echo "Number bed ${numberBed}"

    # Extract the LEFT/RIGHT suffix (the part after the last underscore)
    suffixBed=$(echo "$info" | sed 's/.*_\([A-Za-z]*\)$/\1/')

    echo $suffixBed

    while IFS=$'\t' read -r name pool seq length gc tm ; do

        # Extract the number (the part between the last underscores)
        numberTsv=$(echo "$name" | sed 's/.*_\([0-9]*\)_[A-Za-z]*$/\1/')

        echo "Number TSV ${numberTsv}"

        # Extract the LEFT/RIGHT suffix (the part after the last underscore)
        suffixTsv=$(echo "$name" | sed 's/.*_\([A-Za-z]*\)$/\1/')

        if [[ "$numberTsv" == "$numberBed" && "$suffixTsv" == "$suffixBed" ]]; then

            bed_line="$seqId\t$start\t$end\t$info\t$ph\t$strand\t$seq"

            # Append the line to the output BED file
            echo -e "$bed_line" >> "$output_file"
        fi

    done < "$fileTsv" 

done < "$fileBed"  

echo "BED file created: $output_file"
