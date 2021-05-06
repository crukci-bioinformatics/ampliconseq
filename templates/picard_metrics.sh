#!/bin/bash

set -e -o pipefail


# Picard CollectAlignmentSummaryMetrics bundled with GATK

gatk CollectAlignmentSummaryMetrics \
    --INPUT !{bam} \
    --REFERENCE_SEQUENCE !{reference_sequence} \
    --OUTPUT alignment_metrics.txt

extract_picard_metrics.R \
    --id=!{id} \
    --metrics=alignment_metrics.txt \
    --output=!{id}.alignment_metrics.txt


# extract amplicon and target intervals in BED format and convert to Picard
# interval list format

awk 'BEGIN { FS = "\t"; OFS = "\t" } FNR > 1 { print $2, $3, $4, $1 }' amplicon_groups.txt > amplicons.bed
gatk BedToIntervalList \
    --INPUT amplicons.bed \
    --SEQUENCE_DICTIONARY !{reference_sequence_dictionary} \
    --OUTPUT amplicons.interval_list.txt

awk 'BEGIN { FS = "\t"; OFS = "\t" } FNR > 1 { print $2, $5, $6, $1 }' amplicon_groups.txt > targets.bed
gatk BedToIntervalList \
    --INPUT targets.bed \
    --SEQUENCE_DICTIONARY !{reference_sequence_dictionary} \
    --OUTPUT targets.interval_list.txt


# Picard CollectTargetPcrMetrics bundled with GATK

gatk CollectTargetedPcrMetrics \
    --INPUT !{bam} \
    --REFERENCE_SEQUENCE !{reference_sequence} \
    --AMPLICON_INTERVALS amplicons.interval_list.txt \
    --TARGET_INTERVALS targets.interval_list.txt \
    --OUTPUT targeted_pcr_metrics.txt

extract_picard_metrics.R \
    --id=!{id} \
    --metrics=targeted_pcr_metrics.txt \
    --output=!{id}.targeted_pcr_metrics.txt

