#!/bin/bash

set -e -o pipefail

for group in `sed 1d !{amplicon_groups} | cut -f7 | sort -un`
do
    awk -v group=${group} 'BEGIN { FS = "\t"; OFS = "\t" } $7 == group { print $2, $5 - 1, $6, $1 }' !{amplicon_groups} > targets.bed

    if [ "!{variant_caller}" == "vardict" ]; then

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
            | var2vcf_valid.pl -N "!{id}" -E -P 0 -f !{params.minimumAlleleFraction} \
            > variants.${group}.vcf

    elif [ "!{variant_caller}" == "haplotypecaller" ]; then

        gatk --java-options "-Xmx!{java_mem}m" HaplotypeCaller \
            --input "!{id}.${group}.bam" \
            --intervals targets.bed \
            --reference !{reference_sequence} \
            --output haplotypecaller.${group}.vcf \
            --max-reads-per-alignment-start !{params.maximumReadsPerAlignmentStart} \
            --native-pair-hmm-threads 1 \
            --force-active

        gatk --java-options "-Xmx!{java_mem}m" VariantFiltration \
            --variant haplotypecaller.${group}.vcf \
            --reference !{reference_sequence} \
            --filter-name "QualByDepth" --filter-expression "QD < 2.0" \
            --filter-name "StrandOddsRatio" --filter-expression "SOR > 3.0" \
            --filter-name "RMSMappingQuality" --filter-expression "MQ < 40.0" \
            --filter-name "MappingQualityRankSumTest" --filter-expression "MQRankSum < -12.5" \
            --output variants.${group}.vcf

    else
        echo "Unrecognized variant caller: !{variant_caller}" >&2
        exit 1
    fi

    JAVA_OPTS="-Xmx!{java_mem}m" annotate-vcf-with-amplicon-ids \
        --input variants.${group}.vcf \
        --amplicon-intervals targets.bed \
        --output variants.annotated.${group}.vcf

    echo variants.annotated.${group}.vcf >> vcf_list.txt
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

