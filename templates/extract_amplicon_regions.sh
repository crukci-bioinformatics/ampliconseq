#!/bin/bash

for group in `sed 1d !{amplicon_groups} | cut -f8 | sort -un`
do
	awk -v group=${group} 'BEGIN { FS = "\t"; OFS = "\t" } $8 == group { print $2, $3 - 1, $4, $1 }' !{amplicon_groups} > amplicons.bed

	extract-amplicon-regions \
		--id !{id} \
		--input !{bam} \
		--intervals amplicons.bed \
		--output !{id}.${group}.bam \
		--coverage !{id}.${group}.amplicon_coverage.txt \
		--maximum-distance 0 \
		--require-both-ends-anchored \
		--unmark-duplicate-reads
done

awk 'NR == 1 || FNR > 1' ${id}.*.amplicon_coverage.txt > ${id}.amplicon_coverage.txt
