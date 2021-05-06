#!/bin/bash

set -e -o pipefail

for group in `sed 1d !{amplicon_groups} | cut -f8 | sort -un`
do
    awk -v group=${group} 'BEGIN { FS = "\t"; OFS = "\t" } $8 == group { print $2, $5 - 1, $6, $1 }' !{amplicon_groups} > targets.bed

    JAVA_OPTS="-Xmx!{java_mem}m" vardict-java \
        -b "!{id}.${group}.bam" \
        -G !{reference_sequence} \
        -N "!{id}" \
        -f !{params.minimumAlleleFraction} \
        -z -c 1 -S 2 -E 3 -g 4 targets.bed \
         > vardict.${group}.txt

    cat vardict.${group}.txt \
        | teststrandbias.R \
        > vardict.teststrandbias.${group}.txt

    cat vardict.teststrandbias.${group}.txt \
        | var2vcf_valid.pl -N "!{id}" -E -f !{params.minimumAlleleFraction} \
        > vardict.${group}.vcf

    JAVA_OPTS="-Xmx!{java_mem}m" annotate-vcf-with-amplicon-ids \
        --input vardict.${group}.vcf \
        --amplicon-intervals targets.bed \
        --output vardict.annotated.${group}.vcf

    echo vardict.annotated.${group}.vcf >> vcf_list.txt
done

gatk --java-options "-Xmx!{java_mem}m" MergeVcfs \
    --INPUT vcf_list.txt \
    --SEQUENCE_DICTIONARY !{reference_sequence_dictionary} \
    --OUTPUT "!{id}.vcf"

gatk --java-options "-Xmx!{java_mem}m" VariantsToTable \
    --variant "!{id}.vcf" \
    --output "!{id}.variant_table.txt" \
    --show-filtered \
    --split-multi-allelic \
    --fields AMPLICON \
    --fields CHROM \
    --fields POS \
    --fields REF \
    --fields ALT \
    --fields QUAL \
    --fields FILTER \
    --fields TYPE \
    --fields MULTI-ALLELIC \
    --genotype-fields DP \
    --asGenotypeFieldsToTake AD

echo -e "ID\tAmplicon\tChromosome\tPosition\tRef\tAlt\tQuality\tFilters\tType\tMultiallelic\tDepth\tRef count\tAlt count" > "!{id}.variants.txt"
awk 'BEGIN { FS = "\t"; OFS = "\t" } FNR > 1 && $1 != "NA" { split($11, ad, ","); print "!{id}", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, ad[1], ad[2] }' "!{id}.variant_table.txt" >> "!{id}.variants.txt"

