#!/usr/bin/env nextflow

// enable DSL 2 syntax
nextflow.enable.dsl = 2


// -----------------------------------------------------------------------------
// show settings and/or help
// -----------------------------------------------------------------------------

printParameterSummary()

if (params.help) {
    helpMessage()
    exit 0
}


// -----------------------------------------------------------------------------
// functions
// -----------------------------------------------------------------------------

// returns the heap size to be given to a Java task based on the task memory,
// allowing for some overhead for the JVM
def javaMemMB(task)
{
    return task.memory.toMega() - params.jvmOverhead
}


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
        check_samples_file.R \
            --samples=${sample_sheet} \
            --output=${checked_sample_sheet}
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
        create_non_overlapping_amplicon_groups.R \
            --amplicons=${amplicon_details} \
            --reference-sequence-index=${reference_sequence_index} \
            --output=${amplicon_groups}
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
        java_mem = javaMemMB(task)
        template "extract_amplicon_regions.sh"
}


// generate allele counts for each position within each amplicon
process pileup_counts {
    tag "${id}"

    memory { 4.GB * task.attempt }
    time { 1.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(id), path(bam), path(bai), path(amplicon_groups), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        path "${id}.pileup.txt"

    shell:
        java_mem = javaMemMB(task)
        template "pileup_counts.sh"
}


// call variants
process call_variants {
    tag "${id}"
    publishDir "${params.outputDir}/vcf", mode: "copy", pattern: "${id}.vcf"

    memory { 4.GB * task.attempt }
    time { 1.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(id), path(bam), path(bai), path(amplicon_groups), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        path "${id}.vcf", emit: vcf
        path "${id}.variants.txt", emit: variants

    shell:
        java_mem = javaMemMB(task)
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
        java_mem = javaMemMB(task)
        template "picard_metrics.sh"
}


// create alignment coverage report and combined alignment/coverage metrics table
process alignment_coverage_report {
    executor "local"
    publishDir "${params.outputDir}", mode: 'copy'

    input:
        path samples
        path alignment_metrics
        path targeted_pcr_metrics
        path amplicon_coverage

    output:
        path alignment_coverage_metrics
        path alignment_coverage_report

    script:
        alignment_coverage_metrics = "alignment_coverage_metrics.csv"
        alignment_coverage_report = "alignment_coverage_report.html"
        """
        alignment_coverage_report.R \
            --samples=${samples} \
            --alignment-metrics=${alignment_metrics} \
            --targeted-pcr-metrics=${targeted_pcr_metrics} \
            --amplicon-coverage=${amplicon_coverage} \
            --output-metrics=${alignment_coverage_metrics} \
            --output-report=${alignment_coverage_report}
        """
}

// fit distributions for substitution allele fractions from pileup counts and
// compute background noise thresholds
process compute_background_noise_thresholds {
    publishDir "${params.outputDir}", mode: 'copy'

    input:
        path pileup_counts

    output:
        path position_thresholds
        path dataset_thresholds

    script:
        position_thresholds = "position_noise_thresholds.txt"
        dataset_thresholds = "dataset_noise_thresholds.txt"
        """
        compute_background_noise_thresholds.R \
            --pileup-counts=${pileup_counts} \
            --position-thresholds=${position_thresholds} \
            --dataset-thresholds=${dataset_thresholds} \
            --minimum-depth=100 \
            --exclude-highest-fraction=0.1 \
            --maximum-allele-fraction=0.03 \
            --minimum-number-for-fitting=10 \
            --chunk-size=500000 \
            --read-chunk-size=100000
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
    samples = check_samples(sample_sheet)

    // create groups of non-overlapping amplicons
    amplicon_groups = create_non_overlapping_amplicon_groups(
        amplicon_details,
        reference_sequence_index
    )

    // BAM file channel created by reading the validated sample sheet
    bam = samples
        .splitCsv(header: true, quote: '"')
        .map { row -> tuple("${row.ID}", file("${params.bamDir}/${row.ID}.bam", checkIfExists: true)) }

    // Picard alignment summary metrics
    picard_metrics(bam.combine(amplicon_groups).combine(reference_sequence))

    // collect Picard metrics for all samples
    alignment_metrics = picard_metrics.out.alignment_metrics
        .collectFile(name: "alignment_metrics.txt", keepHeader: true)
    targeted_pcr_metrics = picard_metrics.out.targeted_pcr_metrics
        .collectFile(name: "targeted_pcr_metrics.txt", keepHeader: true)

    // extract reads matching amplicons into subset BAM files for
    // each group of non-overlapping amplicons
    extract_amplicon_regions(bam.combine(amplicon_groups))

    // collect amplicon coverage data for all samples
    amplicon_coverage = extract_amplicon_regions.out.coverage
        .collectFile(name: "amplicon_coverage.txt", keepHeader: true)

    // alignment coverage report
    alignment_coverage_report(samples, alignment_metrics, targeted_pcr_metrics, amplicon_coverage)

    // call variants with VarDict
    call_variants(extract_amplicon_regions.out.bam.combine(reference_sequence))

    // collect variant calls for all samples
    collected_variants = call_variants.out.variants
        .collectFile(name: "variants.txt", keepHeader: true)

    // generate pileup counts
    pileup_counts(extract_amplicon_regions.out.bam.combine(reference_sequence))

    // collect pileup counts for all samples
    collected_pileup_counts = pileup_counts.out
        .collectFile(name: "pileup_counts.txt", keepHeader: true)

    // fit distributions for substitution allele fractions from pileup counts
    // and compute background noise thresholds
    compute_background_noise_thresholds(collected_pileup_counts)
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
        VEP cache directory       : ${params.vepCacheDir}
        Species                   : ${params.vepSpecies}
        Assembly                  : ${params.vepAssembly}
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
            --help                     Show this message and exit
            --sampleSheet              CSV/TSV file containing details of sample datasets (ID and Sample columns required)
            --bamDir                   Directory in which BAM files are located
            --ampliconDetails          CSV/TSV file containing details of the amplicons (ID, Chromosome, AmpliconStart, AmpliconEnd, TargetStart, TargetEnd, Gene columns required)
            --referenceGenomeFasta     FASTA file containing the reference genome sequence (must be indexed, i.e. have an accompanying .fai file)
            --minimumAlleleFraction    Lower allele fraction limit for detection of variants (for variant callers that provide this option only)
            --vepCacheDir              Directory in which to install Ensembl VEP cache files
            --vepSpecies               The species name, e.g. homo_sapiens
            --vepAssembly              The genome assembly, e.g. GRCh37
            --outputDir                Directory to which output files are written
            --outputPrefix             Prefix for output file names

        Alternatively, override settings using a configuration file such as the
        following:

        params {
            sampleSheet           = "samples.csv"
            bamDir                = "bam"
            ampliconDetails       = "amplicons.csv"
            referenceGenomeFasta  = "/reference_data/GRCh37.fa"
            minimumAlleleFraction = 0.01
            vepCacheDir           = "vep_cache"
            vepSpecies            = "homo_sapiens"
            vepAssembly           = "GRCh37"
            outputDir             = "results"
            outputPrefix          = ""
        }

        and run as follows:
            nextflow run crukci-bioinformatics/ampliconseq -c ampliconseq.config
    """.stripIndent()
    log.info ""
}
