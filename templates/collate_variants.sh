#!/bin/bash

set -e -o pipefail

for amplicon_vcf in !{amplicon_vcfs}
do
    echo ${amplicon_vcf} >> vcf_list.txt
done

gatk --java-options "-Xmx!{java_mem}m" MergeVcfs \
    --INPUT vcf_list.txt \
    --SEQUENCE_DICTIONARY !{reference_sequence_dictionary} \
    --OUTPUT "!{vcf}"

add-assorted-annotations-to-vcf \
    --input "!{vcf}" \
    --output merged_and_annotated_variants.vcf \
    --reference-sequence !{reference_sequence} \
    --sequence-context-length 1

gatk --java-options "-Xmx!{java_mem}m" VariantsToTable \
    --variant merged_and_annotated_variants.vcf \
    --output merged_and_annotated_variants.txt \
    --show-filtered \
    --split-multi-allelic \
    --fields AMPLICON \
    --fields CHROM \
    --fields POS \
    --fields REF \
    --fields ALT \
    --fields MULTI-ALLELIC \
    --fields TYPE \
    --fields FILTER \
    --fields QUAL \
    --fields FivePrimeContext \
    --genotype-fields DP \
    --asGenotypeFieldsToTake AD \
    --asGenotypeFieldsToTake AF

tidy_variant_table.R \
    --input merged_and_annotated_variants.txt \
    --id "!{id}" \
    --output "!{variants}"

