#!/bin/bash

for group in `sed 1d !{amplicon_groups} | cut -f8 | sort -un`
do
	awk -v group=${group} 'BEGIN { FS = "\t"; OFS = "\t" } $8 == group { print $2, $5 - 1, $6, $1 }' !{amplicon_groups} > targets.bed

    vardict-java \
        -b !{id}.${group}.bam \
        -G !{reference_sequence} \
        -N !{id} \
        -f !{params.minimumAlleleFraction} \
        -z -c 1 -S 2 -E 3 -g 4 targets.bed \
         > !{id}.${group}.vardict.txt

    cat !{id}.${group}.vardict.txt \
        | teststrandbias.R \
        > !{id}.${group}.vardict.teststrandbias.txt

    cat !{id}.${group}.vardict.teststrandbias.txt \
        | var2vcf_valid.pl -N !{id} -E -f !{params.minimumAlleleFraction} \
        > !{id}.${group}.vardict.vcf

    annotate-vcf-with-amplicon-ids \
        --input !{id}.${group}.vardict.vcf \
        --intervals targets.bed \
        --output !{id}.${group}.vardict.annotated.vcf

    echo "!{id}.${group}.vardict.annotated.vcf" >> vcf_list.txt
done

gatk MergeVcfs \
    --INPUT vcf_list.txt \
    --SEQUENCE_DICTIONARY !{reference_sequence_dictionary} \
    --OUTPUT !{id}.vardict.vcf
