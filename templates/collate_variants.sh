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

gatk --java-options "-Xmx!{java_mem}m" VariantsToTable \
    --variant "!{vcf}" \
    --output variant_table.txt \
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
    --genotype-fields DP \
    --asGenotypeFieldsToTake AD \
    --asGenotypeFieldsToTake AF

tidy_variant_table.R \
    --input variant_table.txt \
    --id "!{id}" \
    --output "!{variants}"

