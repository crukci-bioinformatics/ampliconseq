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
    executor "local"

    input:
        path sample_sheet

    output:
        path checked_sample_sheet

    script:
        checked_sample_sheet = "samples.checked.csv"
        """
        check_samples_file.R ${sample_sheet} ${checked_sample_sheet}
        """
}


// create non-overlapping amplicon groups where none of the amplicons overlap
// with another amplicon from the same group
process create_non_overlapping_amplicon_groups {
    executor "local"

    input:
        path amplicon_details
        path reference_sequence_index

    output:
        path amplicon_groups

    script:
        amplicon_groups = "amplicon_groups.txt"
        """
        create_non_overlapping_amplicon_groups.R ${amplicon_details} ${reference_sequence_index} ${amplicon_groups}
        """
}


// extract reads that correspond to groups of amplicon to create
// subset BAM files
process extract_amplicon_regions {
    tag "${id}"

    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(id), path(bam), path(amplicon_groups)

    output:
        tuple val(id), path("${id}.*.bam"), path("${id}.*.bai"), path(amplicon_groups), emit: bam
        path "${id}.amplicon_coverage.txt", emit: coverage

    shell:
        template "extract_amplicon_regions.sh"
}


// generate allele counts for each position within each amplicon
process pileup_counts {
    tag "${id}"

    cpus 2
    memory { 4.GB * task.attempt }
    time { 1.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(id), path(bam), path(bai), path(amplicon_groups), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        path "${id}.pileup.txt"

    shell:
        template "pileup_counts.sh"
}


// call variants
process call_variants {
    tag "${id}"
    publishDir "${params.outputDir}/vcf", mode: "copy", pattern: "${id}.vcf"

    cpus 2
    memory { 4.GB * task.attempt }
    time { 1.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(id), path(bam), path(bai), path(amplicon_groups), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        path "${id}.vcf", emit: vcf
        path "${id}.variants.txt", emit: variants

    shell:
        template "vardict.sh"
}


// Picard metrics
process picard_metrics {
    tag "${id}"

    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(id), path(bam), path(amplicon_groups), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        path "${id}.alignment_metrics.txt", emit: alignment_metrics
        path "${id}.targeted_pcr_metrics.txt", emit: targeted_pcr_metrics

    shell:
        template "picard_metrics.sh"
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
    samples = check_samples(sample_sheet)

    // create groups of non-overlapping amplicons
    amplicon_groups = create_non_overlapping_amplicon_groups(
        amplicon_details,
        reference_sequence_index
    )

    // BAM file channel created by reading the validated sample sheet
    bam = check_samples.out
        .splitCsv(header: true, quote: '"')
        .map { row -> tuple("${row.ID}", file("${params.bamDir}/${row.ID}.bam", checkIfExists: true)) }

    // Picard alignment summary metrics
    picard_metrics(bam.combine(amplicon_groups).combine(reference_sequence))

    // collect Picard metrics for all samples
    collected_alignment_metrics = picard_metrics.out.alignment_metrics
        .collectFile(name: "alignment_metrics.txt", keepHeader: true)
    collected_alignment_metrics = picard_metrics.out.targeted_pcr_metrics
        .collectFile(name: "targeted_pcr_metrics.txt", keepHeader: true)

    // extract reads matching amplicons into subset BAM files for
    // each group of non-overlapping amplicons
    extract_amplicon_regions(bam.combine(amplicon_groups))

    // collect amplicon coverage data for all samples
    collected_amplicon_coverage = extract_amplicon_regions.out.coverage
        .collectFile(name: "amplicon_coverage.txt", keepHeader: true)

    // generate pileup counts
    pileup_counts(extract_amplicon_regions.out.bam.combine(reference_sequence))

    // collect pileup counts for all samples
    collected_pileup_counts = pileup_counts.out
        .collectFile(name: "pileup_counts.txt", keepHeader: true)

    // call variants with VarDict
    call_variants(extract_amplicon_regions.out.bam.combine(reference_sequence))

    // collect variant calls for all samples
    collected_variants = call_variants.out.variants
        .collectFile(name: "variants.txt", keepHeader: true)
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
