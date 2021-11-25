#!/bin/bash

set -e -o pipefail

JAVA_OPTS="-Xmx!{java_mem}m" extract-amplicon-regions \
    --id "!{id}" \
    --input !{bam} \
    --amplicon-intervals !{amplicon_bed} \
    --output "!{amplicon_bam}" \
    --coverage "!{amplicon_coverage}" \
    --maximum-distance !{params.maxDistanceFromAmpliconEnd} \
    --require-both-ends-anchored=!{params.requireBothEndsAnchored} \
    --unmark-duplicate-reads

