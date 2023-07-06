#!/bin/bash

snpeff=$1

warnings=("INFO_RELALIGN_3_PRIME"
        "WARNING_SEQUENCE_NOT_AVAILABLE"
        "WARNING_REF_DOES_NOT_MATCH_GENOME"
        "WARNING_TRANSCRIPT_INCOMPLETE"
        "WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS"
        "WARNING_TRANSCRIPT_NO_START_CODON"
        "ERROR_CHROMOSOME_NOT_FOUND"
        "ERROR_OUT_OF_CHROMOSOME_RANGE"
        "ERROR_OUT_OF_EXON"
        "ERROR_MISSING_CDS_SEQUENCE"
        )

warning_count=0
output="$(echo "$snpeff" | sed 's/.vcf//')"
sample="$(echo "$snpeff" | sed 's/_.*//' )"
echo -e "sample:\t${sample}\nwarning\tcount" > "${output}_warnings.txt"
for warning in "${warnings[@]}"
do
    echo -en "${warning}\t" >> "${output}_warnings.txt"
    count=$(bcftools view "$1" | grep "${warning}" | wc -l)
    warning_count=$((warning_count+count))
    echo $count >> "${output}_warnings.txt"
done
total="$(bcftools view -H "$1" | wc -l)"
percent=$(awk "BEGIN {print $warning_count/$total}")
echo -en "Percent warnings:\t${percent}" >> "${output}_warnings.txt"
