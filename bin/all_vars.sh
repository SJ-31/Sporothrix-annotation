#!/bin/bash

sample=$1
to_read=$2
vcf=$3

counts="${sample}_num_vars.tsv"
while IFS=$'\t' read -r name loc
do
    filename="${sample}_${name}.vcf"
    if [[ -z "$loc" ]]
    then
        touch "${filename}EMPTY"
        continue
    fi
    bcftools view -r "$loc" "$vcf" > "$filename"
    echo -n "$name" >> "temp"
    bcftools plugin counts "$filename" | awk -F ":" \
        '{gsub(/ /, "", $2)}
        $1 !~ /samples/ {printf "\t%s", $2 >> "temp"}
        END{print "" >> "temp"}'
done<"$to_read"
echo -e "Gene ID\tSNPs\tINDELs\tMNPs\tOthers\tTotal" > "$counts"
sort -nrk6 temp >> "$counts"
rm temp

