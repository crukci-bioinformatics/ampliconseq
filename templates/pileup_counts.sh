#!/bin/bash

set -e -o pipefail

for group in `sed 1d !{amplicon_groups} | cut -f8 | sort -un`
do
    awk -v group=${group} 'BEGIN { FS = "\t"; OFS = "\t" } $8 == group { print $2, $5 - 1, $6, $1 }' !{amplicon_groups} > targets.bed

    JAVA_OPTS="-Xmx!{java_mem}m" pileup-counts \
        --id "!{id}" \
        --input "!{id}.${group}.bam" \
        --amplicon-intervals targets.bed \
        --reference-sequence !{reference_sequence} \
        --output "pileup_counts.${group}.txt" \
        --minimum-mapping-quality !{params.minimumMappingQualityForPileup} \
        --minimum-base-quality !{params.minimumBaseQualityForPileup}
done

collate_and_sort_amplicon_coverage.R \
    --sample "!{sample}" \
    --amplicons !{amplicon_groups} \
    --output "!{pileup}" \
    pileup_counts.*.txt

