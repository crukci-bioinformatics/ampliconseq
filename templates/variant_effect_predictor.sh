#!/bin/bash

set -e -o pipefail

create_distinct_vcf.R \
    --input !{variants} \
    --reference-sequence-index !{reference_sequence_index} \
    --output distinct_variants.vcf

vep \
    --input_file distinct_variants.vcf \
    --format vcf \
    --output_file distinct_variants.vep.vcf \
    --vcf \
    --offline \
    --dir_cache !{vep_cache_dir} \
    --species !{params.vepSpecies} \
    --assembly !{params.vepAssembly} \
    --buffer_size 100 \
    --no_stats \
    --dont_skip \
    !{vep_pick_option} \
    --everything \
    --minimal \
    --allele_number 1 \
    --no_escape

vep_vcf_to_tabular.R \
    --vcf distinct_variants.vep.vcf \
    --output !{variant_annotations}

