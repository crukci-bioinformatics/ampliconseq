#!/bin/bash

set -e -o pipefail

create_distinct_vcf.R \
    --variants !{variants} \
    --output distinct_variants.vcf

add-assorted-annotations-to-vcf \
    --input distinct_variants.vcf \
    --output distinct_variants.annotated.vcf \
    --reference-sequence !{reference_sequence} \
    --sequence-context-length !{params.sequenceContextLength}

gatk --java-options "-Xmx!{java_mem}m" VariantsToTable \
    --variant distinct_variants.annotated.vcf \
    --output "!{variant_annotations}" \
    --fields CHROM \
    --fields POS \
    --fields REF \
    --fields ALT \
    --fields FivePrimeContext \
    --fields ThreePrimeContext \
    --fields IndelLength

