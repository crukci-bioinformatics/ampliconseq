#!/usr/bin/env nextflow

// enable DSL 2 syntax
nextflow.enable.dsl = 2

// -----------------------------------------------------------------------------
// default parameter settings
// -----------------------------------------------------------------------------

params.help                  = false
params.sampleSheet           = "${launchDir}/samples.csv"
params.bamDir                = "${launchDir}/bam"
params.ampliconDetails       = "${launchDir}/reference_data/amplicons.csv"
params.referenceGenomeFasta  = "${launchDir}/reference_data/GRCh37.fa"
params.minimumAlleleFraction = 0.01
params.outputDir             = "${launchDir}"
params.outputPrefix          = ""


printParameterSummary()

if (params.help) {
    helpMessage()
    exit 0
}


// -----------------------------------------------------------------------------
// check/derive parameters
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
// processes
// -----------------------------------------------------------------------------

// check sample sheet is valid a create a validated version to be used in
// subsequent processes
process check_samples {
    executor = 'local'

    input:
        path sample_sheet

    output:
        path "samples.checked.csv"

    script:
        """
        check_samples_file.R ${sample_sheet} samples.checked.csv
        """
}


// create non-overlapping amplicon groups where none of the amplicons overlap
// with another amplicon from the same group
process create_non_overlapping_amplicon_groups {
    executor = 'local'

    input:
        path amplicon_details
        path reference_sequence_index

    output:
        path "amplicon_groups.txt", emit: amplicons
        path "amplicon_groups.*.txt", emit: amplicon_groups

    script:
        """
        create_non_overlapping_amplicon_groups.R ${amplicon_details} ${reference_sequence_index}
        """
}


// extract reads that correspond to a set of amplicon intervals to create a
// subset BAM file
process extract_amplicon_regions {
    tag "${prefix}"

    memory = { 2.GB * task.attempt }
    time = { 1.hour * task.attempt }
    maxRetries = 2

    input:
        tuple val(id), path(bam), val(amplicon_group), path(amplicons)

    output:
        tuple val(id), path(output_bam), path(output_bai), val(amplicon_group), path(amplicons), emit: bam
        path(coverage), emit: coverage

    script:
        prefix = "${id}.${amplicon_group}"
        output_bam = "${prefix}.bam"
        output_bai = "${prefix}.bai"
        coverage = "${prefix}.amplicon_coverage.txt"
        """
        awk 'BEGIN { FS = "\t"; OFS = "\t" } FNR > 1 { print \$2, \$3 - 1, \$4, \$1 }' ${amplicons} > amplicons.bed

        extract-amplicon-regions \
            --id ${id} \
            --input ${bam} \
            --intervals amplicons.bed \
            --output ${output_bam} \
            --coverage ${coverage} \
            --maximum-distance 0 \
            --require-both-ends-anchored \
            --unmark-duplicate-reads
        """
}


// generates pileup counts table for each position in the given set of amplicon
// intervals
process pileup_counts {
    tag "${prefix}"

    cpus = 2
    memory = { 4.GB * task.attempt }
    time = { 1.hour * task.attempt }
    maxRetries = 2

    input:
        tuple val(id), path(bam), path(bai), val(amplicon_group), path(amplicons), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        path(pileup)

    script:
        prefix = "${id}.${amplicon_group}"
        pileup = "${prefix}.pileup.txt"
        """
        awk 'BEGIN { FS = "\t"; OFS = "\t" } FNR > 1 { print \$2, \$5 - 1, \$6, \$1 }' ${amplicons} > targets.bed

        pileup-counts \
            --id ${id} \
            --input ${bam} \
            --intervals targets.bed \
            --reference-sequence ${reference_sequence} \
            --output ${pileup}
        """
}


// runs the specified variant caller
process call_variants {
    tag "${prefix}"

    cpus = 2
    memory = { 4.GB * task.attempt }
    time = { 1.hour * task.attempt }
    maxRetries = 2

    input:
        tuple val(id), path(bam), path(bai), val(amplicon_group), path(amplicons), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        tuple val(id), path(vcf)

    script:
        prefix = "${id}.${amplicon_group}"
        vcf = "${prefix}.vcf"
        """
        awk 'BEGIN { FS = "\t"; OFS = "\t" } FNR > 1 { print \$2, \$5 - 1, \$6, \$1 }' ${amplicons} > targets.bed

        vardict-java \
            -b ${bam} \
            -G ${reference_sequence} \
            -N ${id} \
            -f ${params.minimumAlleleFraction} \
            -z -c 1 -S 2 -E 3 -g 4 targets.bed \
            > ${prefix}.vardict.txt

        cat ${prefix}.vardict.txt \
            | teststrandbias.R \
            > ${prefix}.vardict.teststrandbias.txt

        cat ${prefix}.vardict.teststrandbias.txt \
            | var2vcf_valid.pl -N ${id} -E -f ${params.minimumAlleleFraction} \
            > ${prefix}.vardict.vcf

        annotate-vcf-with-amplicon-ids \
            --input ${prefix}.vardict.vcf \
            --intervals targets.bed \
            --output ${vcf}
        """
}


// merge VCF files for a sample from groups of non-overlapping amplicons
process merge_sample_vcfs {
    tag "${id}"
    publishDir "${params.outputDir}/vcf", mode: 'copy'

    input:
        tuple val(id), path(vcfs), path(reference_sequence_dictionary)

    output:
        tuple val(id), path(merged_vcf)

    script:
        merged_vcf = "${id}.vcf"
        """
        echo ${vcfs} | tr ' ' '\n' > input_vcf_list.txt
        gatk MergeVcfs \
            --INPUT input_vcf_list.txt \
            --SEQUENCE_DICTIONARY ${reference_sequence_dictionary} \
            --OUTPUT ${merged_vcf}
        """
}


// Picard alignment summary metrics
process alignment_summary_metrics {
    tag "${id}"

    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(id), path(bam), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        path(metrics)

    script:
        metrics = "${id}.alignment_summary_metrics.txt"
        """
        gatk CollectAlignmentSummaryMetrics \
            --INPUT ${bam} \
            --REFERENCE_SEQUENCE ${reference_sequence} \
            --OUTPUT alignment_summary_metrics.txt

        extract_picard_metrics.R ${id} alignment_summary_metrics.txt ${metrics}
        """
}


// Picard targeted PCR metrics
process targeted_pcr_metrics {
    tag "${id}"

    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(id), path(bam), path(amplicons), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        path(metrics)

    script:
        metrics = "${id}.targeted_pcr_metrics.txt"
        """
        cp ${reference_sequence_dictionary} amplicon_intervals.txt
        awk 'BEGIN { FS = "\t"; OFS = "\t" } FNR > 1 { print \$2, \$3, \$4, "+", \$1 }' ${amplicons} >> amplicon_intervals.txt

        cp ${reference_sequence_dictionary} target_intervals.txt
        awk 'BEGIN { FS = "\t"; OFS = "\t" } FNR > 1 { print \$2, \$5, \$6, "+", \$1 }' ${amplicons} >> target_intervals.txt

        gatk CollectTargetedPcrMetrics \
            --INPUT ${bam} \
            --REFERENCE_SEQUENCE ${reference_sequence} \
            --AMPLICON_INTERVALS amplicon_intervals.txt \
            --TARGET_INTERVALS target_intervals.txt \
            --OUTPUT targeted_pcr_metrics.txt

        extract_picard_metrics.R ${id} targeted_pcr_metrics.txt ${metrics}
        """
}


// -----------------------------------------------------------------------------
// workflow
// -----------------------------------------------------------------------------

workflow {

    sample_sheet = channel.fromPath(params.sampleSheet, checkIfExists: true)
    amplicon_details = channel.fromPath(params.ampliconDetails, checkIfExists: true)

    reference_sequence_fasta_file = file("${params.referenceGenomeFasta}")
    reference_sequence_fasta = channel.fromPath(params.referenceGenomeFasta, checkIfExists: true)
    reference_sequence_index = channel.fromPath("${params.referenceGenomeFasta}.fai", checkIfExists: true)
    reference_sequence_dictionary = channel.fromPath("${params.referenceGenomeFasta}".replaceFirst("${reference_sequence_fasta_file.extension}\$", "dict"), checkIfExists: true)
    reference_sequence = reference_sequence_fasta
        .combine(reference_sequence_index)
        .combine(reference_sequence_dictionary) 

    // check sample sheet
    check_samples(sample_sheet)

    // create groups of non-overlapping amplicons
    create_non_overlapping_amplicon_groups(
        amplicon_details,
        reference_sequence_index
    )

    // amplicon group channel in which each value is a tuple that includes
    // the amplicon file for that group and the group number extracted from the
    // amplicon file name
    amplicon_groups = create_non_overlapping_amplicon_groups.out.amplicon_groups
        .flatten()
        .map { file -> tuple( (file.name =~ /(\d+)/)[0][1] as Integer, file ) }

    // BAM file channel created by reading the validated sample sheet
    bam = check_samples.out
        .splitCsv(header: true, quote: '"')
        .map { row -> tuple("${row.ID}", file("${params.bamDir}/${row.ID}.bam", checkIfExists: true)) }

    // Picard alignment summary metrics
    alignment_summary_metrics(bam.combine(reference_sequence))
    collected_alignment_summary_metrics = alignment_summary_metrics.out
        .collectFile(name: "alignment_summary_metrics.txt", keepHeader: true)

    // Picard targeted PCR metrics
    targeted_pcr_metrics(
        bam
            .combine(create_non_overlapping_amplicon_groups.out.amplicons)
            .combine(reference_sequence)
    )
    collected_targeted_pcr_metrics = targeted_pcr_metrics.out
        .collectFile(name: "targeted_pcr_metrics.txt", keepHeader: true)

    // extract reads matching amplicons into a subset BAM file
    // for all pairs of BAM files and amplicon groups
    extract_amplicon_regions(bam.combine(amplicon_groups))
    collected_amplicon_coverage = extract_amplicon_regions.out.coverage
        .collectFile(name: "amplicon_coverage.txt", keepHeader: true)

    // generate pileup counts
    pileup_counts(
        extract_amplicon_regions.out.bam
            .combine(reference_sequence)
    )
    collected_pileup_counts = pileup_counts.out
        .collectFile(name: "pileup.txt", keepHeader: true)

    // VarDict variant calling
    call_variants(
        extract_amplicon_regions.out.bam
            .combine(reference_sequence)
    )

    // merge VCF files for each sample
    merge_sample_vcfs(
        call_variants.out
            .groupTuple()
            .combine(reference_sequence_dictionary)
    )

    // concatenate coverage files for each BAM file
    // amplicon_coverage = extract_amplicon_regions.out.coverage
    //     .collectFile(keepHeader: true) { id, coverage -> [ "${id}.amplicon_coverage.txt", coverage ] }
}


// -----------------------------------------------------------------------------
// summary of configuration parameters
// -----------------------------------------------------------------------------

def printParameterSummary() {
    log.info ""
    log.info """
        Variant calling pipeline for amplicon sequencing data
        =====================================================

        Sample sheet              : ${params.sampleSheet}
        BAM directory             : ${params.bamDir}
        Amplicon details          : ${params.ampliconDetails}
        Reference genome sequence : ${params.referenceGenomeFasta}
        Minimum allele fraction   : ${params.minimumAlleleFraction}
        Output directory          : ${params.outputDir}
        Output prefix             : ${params.outputPrefix}
    """.stripIndent()
    log.info ""
}


// ----------------------------------------------------------------------------
// help/usage
// ----------

def helpMessage() {
    log.info """
        Usage:
            nextflow run crukci-bioinformatics/ampliconseq

        Options:
            --help                        Show this message and exit
            --sample-sheet                CSV/TSV file containing details of sample datasets (ID and Sample columns required)
            --bam-dir                     Directory in which BAM files are located
            --amplicon-details            CSV/TSV file containing details of the amplicons (ID, Chromosome, AmpliconStart, AmpliconEnd, TargetStart, TargetEnd, Gene columns required)
            --reference-genome-fasta      FASTA file containing the reference genome sequence (must be indexed, i.e. have an accompanying .fai file)
            --minimum-allele-fraction     Lower allele fraction limit for detection of variants (for variant callers that provide this option only)
            --output-dir                  Directory to which output files are written
            --output-prefix               Prefix for output file names

        Alternatively, override settings using a configuration file such as the
        following, in which parameter names used are the camelCase equivalent of the
        above options:

        params {
            sampleSheet           = "samples.csv"
            bamDir                = "bam"
            ampliconDetails       = "amplicons.csv"
            referenceGenomeFasta  = "/reference_data/GRCh37.fa"
            minimumAlleleFraction = 0.01
            outputDir             = "results"
            outputPrefix          = ""
        }

        and run as follows:
            nextflow run crukci-bioinformatics/ampliconseq -c ampliconseq.config
    """.stripIndent()
    log.info ""
}
