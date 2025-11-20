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
    --output merged_annotated.vcf \
    --reference-sequence !{reference_sequence} \
    --sequence-context-length 1

# left-align and normalize indels and split multi-allelic variants
# note that this can still result in deletions that have an '*' for the
# alt allele which will need fixing
bcftools norm \
    --multiallelics -both \
    --fasta-ref !{reference_sequence} \
    merged_annotated.vcf \
    --output merged_annotated_normalized.vcf

gatk --java-options "-Xmx!{java_mem}m" VariantsToTable \
    --variant merged_annotated_normalized.vcf \
    --output merged_annotated_normalized.txt \
    --show-filtered \
    --fields AMPLICON \
    --fields CHROM \
    --fields POS \
    --fields REF \
    --fields ALT \
    --fields MULTIALLELIC \
    --fields FILTER \
    --fields QUAL \
    --fields FivePrimeContext \
    --genotype-fields DP \
    --asGenotypeFieldsToTake AD \
    --asGenotypeFieldsToTake AF

tidy_variant_table.R \
    --input merged_annotated_normalized.txt \
    --id "!{id}" \
    --output "!{variants}"

