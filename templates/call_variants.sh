#!/bin/bash

set -e -o pipefail

if [ "!{variant_caller}" == "vardict" ]; then

    JAVA_OPTS="-Xmx!{java_mem}m" vardict-java \
        -b !{amplicon_bam} \
        -G !{reference_sequence} \
        -N "!{id}" \
        -f !{params.minimumAlleleFraction} \
        -z -c 1 -S 2 -E 3 -g 4 !{target_bed} \
        > vardict.txt

    cat vardict.txt \
        | teststrandbias.R \
        > vardict.teststrandbias.txt

    cat vardict.teststrandbias.txt \
        | var2vcf_valid.pl -N "!{id}" -E -P 0 -f !{params.minimumAlleleFraction} \
        > variants.vcf

elif [ "!{variant_caller}" == "haplotypecaller" ]; then

    gatk --java-options "-Xmx!{java_mem}m" HaplotypeCaller \
        --input !{amplicon_bam} \
        --intervals !{target_bed} \
        --reference !{reference_sequence} \
        --output haplotypecaller.vcf \
        --max-reads-per-alignment-start !{params.maximumReadsPerAlignmentStart} \
        --native-pair-hmm-threads 1 \
        --force-active

    gatk --java-options "-Xmx!{java_mem}m" VariantFiltration \
        --variant haplotypecaller.vcf \
        --reference !{reference_sequence} \
        --filter-name "QualByDepth" --filter-expression "QD < 2.0" \
        --filter-name "StrandOddsRatio" --filter-expression "SOR > 3.0" \
        --filter-name "RMSMappingQuality" --filter-expression "MQ < 40.0" \
        --filter-name "MappingQualityRankSumTest" --filter-expression "MQRankSum < -12.5" \
        --output variants.vcf

else
    echo "Unrecognized variant caller: !{variant_caller}" >&2
    exit 1
fi

JAVA_OPTS="-Xmx!{java_mem}m" annotate-vcf-with-amplicon-ids \
    --input variants.vcf \
    --target-intervals !{target_bed} \
    --output "!{vcf}"

