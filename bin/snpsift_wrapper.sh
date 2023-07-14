#!/bin/bash

# name="$1"
path="$2"
snpsift="$3"

# Sort by high impact and high quality
cat "$path" | java -jar "$snpsift" filter "( ANN[*].IMPACT has 'HIGH') \
    & !( exists ANN[*].ERRORS ) \
    & ( QUAL >= 30 )" \
    > "${1}"_filtered.vcf
