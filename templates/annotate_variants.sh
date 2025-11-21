#!/bin/bash

set -e -o pipefail

create_distinct_vcf.R \
    --input !{variants} \
    --reference-sequence-index !{reference_sequence_index} \
    --output distinct_variants.vcf

add-assorted-annotations-to-vcf \
    --input distinct_variants.vcf \
    --output distinct_variants.annotated.vcf \
    --reference-sequence !{reference_sequence} \
    --sequence-context-length !{params.sequenceContextLength}

gatk --java-options "-Xmx!{java_mem}m" VariantsToTable \
    --variant distinct_variants.annotated.vcf \
    --output distinct_variants.annotated.txt \
    --fields CHROM \
    --fields POS \
    --fields REF \
    --fields ALT \
    --fields FivePrimeContext \
    --fields ThreePrimeContext \
    --fields IndelLength

echo -e "Chromosome\tPosition\tRef\tAlt\t5' context\tAlleles\t3' context\tIndel length" > !{other_annotations}
awk 'BEGIN { FS = "\t"; OFS = "\t" } FNR > 1 { print $1, $2, $3, $4, $5, $3"/"$4, $6, $7 }' distinct_variants.annotated.txt >> !{other_annotations}

add_offset_from_primer_end.R \
    --variants !{variants} \
    --amplicons !{amplicons} \
    --output !{offset_from_primer_end_annotations}
