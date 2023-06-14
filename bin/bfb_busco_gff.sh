#!/bin/bash
name=$1
busco_ref=$2
aligned=$3
# This version of the busco gene extractor from aligned bam files
#   requires the output of the busco_to_gff.sh script

while IFS= read -r line
    do
        fields=($(echo $line | awk '! /#/ {print $1,$3,$4,$5,$6,$7}'))
        if [ -z "${fields[4]}" ]; then
            continue
        fi
        region=${fields[1]}:${fields[2]}-${fields[3]}
        length=$((fields[3]-fields[2]))
        header=">${name}_${fields[0]} ${fields[5]} ${fields[4]}"
        echo "$header GFF_length:$length" > \
        "${name}"_"${fields[0]}".fasta
        samtools view "$aligned" "$region" | \
        awk '{print $10}'>> "${name}"_"${fields[0]}".fasta
    done < "$busco_ref"
