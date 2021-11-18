#!/bin/bash

set -e -o pipefail

for group in `sed 1d !{amplicon_groups} | cut -f7 | sort -un`
do
    awk -v group=${group} 'BEGIN { FS = "\t"; OFS = "\t" } $7 == group { print $2, $3 - 1, $4, $1 }' !{amplicon_groups} > amplicons.bed

    JAVA_OPTS="-Xmx!{java_mem}m" extract-amplicon-regions \
        --id "!{id}" \
        --input !{bam} \
        --amplicon-intervals amplicons.bed \
        --output "!{id}.${group}.bam" \
        --coverage "amplicon_coverage.${group}.txt" \
        --maximum-distance !{params.maxDistanceFromAmpliconEnd} \
        --require-both-ends-anchored=!{params.requireBothEndsAnchored} \
        --unmark-duplicate-reads
done

collate_and_sort_amplicon_coverage.R \
    --sample "!{sample}" \
    --amplicons !{amplicon_groups} \
    --output "!{amplicon_coverage}" \
    amplicon_coverage.*.txt

