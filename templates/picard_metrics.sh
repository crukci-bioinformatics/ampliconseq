#!/bin/bash

set -e -o pipefail


# Picard CollectAlignmentSummaryMetrics bundled with GATK

gatk --java-options "-Xmx!{java_mem}m" CollectAlignmentSummaryMetrics \
    --INPUT !{bam} \
    --REFERENCE_SEQUENCE !{reference_sequence} \
    --OUTPUT alignment_metrics.txt

extract_picard_metrics.R \
    --id "!{id}" \
    --metrics alignment_metrics.txt \
    --output "!{alignment_metrics}"


# extract amplicon and target intervals in BED format and convert to Picard
# interval list format

awk 'BEGIN { FS = "\t"; OFS = "\t" } FNR > 1 { print $2, $3, $4, $1 }' !{amplicon_groups} > amplicons.bed
gatk --java-options "-Xmx!{java_mem}m" BedToIntervalList \
    --INPUT amplicons.bed \
    --SEQUENCE_DICTIONARY !{reference_sequence_dictionary} \
    --OUTPUT amplicons.interval_list.txt

awk 'BEGIN { FS = "\t"; OFS = "\t" } FNR > 1 { print $2, $5, $6, $1 }' !{amplicon_groups} > targets.bed
gatk --java-options "-Xmx!{java_mem}m" BedToIntervalList \
    --INPUT targets.bed \
    --SEQUENCE_DICTIONARY !{reference_sequence_dictionary} \
    --OUTPUT targets.interval_list.txt


# Picard CollectTargetedPcrMetrics bundled with GATK

gatk --java-options "-Xmx!{java_mem}m" CollectTargetedPcrMetrics \
    --INPUT !{bam} \
    --REFERENCE_SEQUENCE !{reference_sequence} \
    --AMPLICON_INTERVALS amplicons.interval_list.txt \
    --TARGET_INTERVALS targets.interval_list.txt \
    --OUTPUT targeted_pcr_metrics.txt

extract_picard_metrics.R \
    --id "!{id}" \
    --metrics targeted_pcr_metrics.txt \
    --output "!{targeted_pcr_metrics}"

