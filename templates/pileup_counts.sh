#!/bin/bash

set -e -o pipefail


for group in `sed 1d !{amplicon_groups} | cut -f8 | sort -un`
do
    awk -v group=${group} 'BEGIN { FS = "\t"; OFS = "\t" } $8 == group { print $2, $5 - 1, $6, $1 }' !{amplicon_groups} > targets.bed

    JAVA_OPTS="-Xmx4g" pileup-counts \
        --id !{id} \
        --input !{id}.${group}.bam \
        --intervals targets.bed \
        --reference-sequence !{reference_sequence} \
        --output !{id}.${group}.pileup.txt
done

awk 'NR == 1 || FNR > 1' !{id}.*.pileup.txt > !{id}.pileup.txt

