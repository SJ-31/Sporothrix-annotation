#!/bin/bash
name=$1
busco_ref=$2
aligned=$3
while read line
    do
    fields=($(echo $line | awk '! /#/ && $2 ~ /Complete/ {print $1,$3,$4,$5,$6,$10}'))
    case ${fields[4]} in
    '-')
        region=${fields[1]}:${fields[3]}-${fields[2]};;
    '+')
        region=${fields[1]}:${fields[2]}-${fields[3]};;
    *)
        continue;;
    esac
    echo ">${name}_${fields[0]}:${fields[5]}" > \
    "${name}"_"${fields[0]}".fasta
    samtools view "$aligned" "$region" | \
    awk '{print $10}'>> "${name}"_"${fields[0]}".fasta
    done < "$busco_ref"
