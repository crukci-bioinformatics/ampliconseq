#!/bin/bash

JAVA_OPTS="-Xmx!{java_mem}m" pileup-counts \
    --id "!{id}" \
    --input !{amplicon_bam} \
    --amplicon-intervals !{target_bed} \
    --reference-sequence !{reference_sequence} \
    --output "!{pileup_counts}" \
    --minimum-mapping-quality !{params.minimumMappingQualityForPileup} \
    --minimum-base-quality !{params.minimumBaseQualityForPileup}

