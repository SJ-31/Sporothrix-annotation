#!/bin/bash

# bcftools view "$1" | sed 's/NW_[[:digit:]]*\([[:digit:]]\{2\}\)\.1/Cont\1/' | \
#     java -jar /home/sc31/Bio_SDD/tools/snpEff/snpEff.jar -d \
#     -dataDir "$2" \
#     -nodownload \
#     "$3" > "${4}_annotated_vars.vcf" # Specify the genome

bcftools view "$1" | \
    java -jar /home/sc31/Bio_SDD/tools/snpEff/snpEff.jar -d \
    -dataDir "$2" \
    -nodownload \
    "$3" > "${4}_annotated_vars.vcf" # Specify the genome


