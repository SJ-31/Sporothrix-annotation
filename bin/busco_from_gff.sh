#!/bin/bash
# Extract gff regions as fasta files
while getopts 's:b:f:g:' option
do 
    case "${option}"
        in 
        s) sample=${OPTARG};;
        b) busco_gff_tsv=${OPTARG};;
        f) sample_fasta=${OPTARG};;
        g) sample_gff=${OPTARG};;
        *) echo "Invalid"
            exit 1;;
    esac
done

while IFS= read -r line
    do
        fields=($(echo $line | awk '! /#/ {print $3,$4,$5,$1,$2,$7,$6}'))
        if [ -z "${fields[4]}" ]; then
            continue
        fi
        region="${fields[0]}:${fields[1]}-${fields[2]}"
        header=">${sample}|${fields[3]}|${fields[4]}|${fields[5]}"
        output="${sample}_${fields[3]}.fasta"
        gffread "$sample_gff" -r "$region" -w "temp" -g "$sample_fasta"
        sed "s/>.*/$header/" "temp" > "${output}"
    done < "$busco_gff_tsv"
rm temp

# fields[0] Sequence ID
# fields[1] Start
# fields[2] Stop
# fields[3] Busco ID
# fields[4] NCBI descriptor
# fields[5] Busco descriptor
