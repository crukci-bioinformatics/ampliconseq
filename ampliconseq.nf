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
            --amplicons ${amplicon_details} \
            --reference-sequence-index ${reference_sequence_index} \
            --output ${amplicon_groups}
        """
}


// check input files are valid and create validated versions in CSV format
// (inputs can be either TSV or CSV) for use in subsequent processes
process check_inputs {
    executor "local"

    input:
        path sample_sheet
        path specific_variants
        path blacklisted_variants
        path amplicons

    output:
        path checked_sample_sheet, emit: samples
        path checked_specific_variants, emit: specific_variants
        path checked_blacklisted_variants, emit: blacklisted_variants

    shell:
        checked_sample_sheet = "samples.checked.txt"
        checked_specific_variants = "specific_variants.checked.txt"
        checked_blacklisted_variants = "blacklisted_variants.checked.txt"
        template "check_inputs.sh"
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
        path amplicon_coverage, emit: coverage

    shell:
        java_mem = javaMemMB(task)
        amplicon_coverage = "${id}.amplicon_coverage.txt"
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
        path pileup

    shell:
        java_mem = javaMemMB(task)
        pileup = "${id}.pileup.txt"
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
        path vcf, emit: vcf
        path variants, emit: variants

    shell:
        java_mem = javaMemMB(task)
        vcf = "${id}.vcf"
        variants = "${id}.variants.txt"
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
        path alignment_metrics, emit: alignment_metrics
        path targeted_pcr_metrics, emit: targeted_pcr_metrics

    shell:
        java_mem = javaMemMB(task)
        alignment_metrics = "${id}.alignment_metrics.txt"
        targeted_pcr_metrics = "${id}.targeted_pcr_metrics.txt"
        template "picard_metrics.sh"
}


// create alignment coverage report and combined alignment/coverage metrics table
process alignment_coverage_report {
    executor "local"
    publishDir "${params.outputDir}", mode: "copy"

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
            --samples ${samples} \
            --alignment-metrics ${alignment_metrics} \
            --targeted-pcr-metrics ${targeted_pcr_metrics} \
            --amplicon-coverage ${amplicon_coverage} \
            --output-metrics ${alignment_coverage_metrics} \
            --output-report ${alignment_coverage_report}
        """
}


// fit distributions for substitution allele fractions from pileup counts and
// compute background noise thresholds
process compute_background_noise_thresholds {
    publishDir "${params.outputDir}", mode: "copy"

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
            --pileup-counts ${pileup_counts} \
            --position-thresholds ${position_thresholds} \
            --dataset-thresholds ${dataset_thresholds} \
            --minimum-depth ${params.minimumDepthForBackgroundNoise} \
            --exclude-highest-fraction ${params.excludeHighestFractionForBackgroundNoise} \
            --maximum-allele-fraction ${params.maximumAlleleFractionForBackgroundNoise} \
            --minimum-number-for-fitting ${params.minimumNumberForFittingBackgroundNoise} \
            --chunk-size ${params.chunkSizeForFittingBackgroundNoise} \
            --read-chunk-size ${params.readChunkSizeForFittingBackgroundNoise}
        """
}


// expand variant table to include specific variants and missing calls within
// sample replicates
process add_specific_variants {
    executor "local"

    input:
        path samples
        path called_variants
        path specific_variants

    output:
        path all_variants

    script:
        all_variants = "all_variants.txt"
        """
        add_specific_variants.R \
            --samples ${samples} \
            --called-variants ${called_variants} \
            --specific-variants ${specific_variants} \
            --output ${all_variants}
        """
}


// annotate variants using Ensembl VEP
process variant_effect_predictor {

    input:
        path variants
        path vep_cache_dir

    output:
        path variant_annotations

    shell:
        variant_annotations = "vep_annotations.txt"
        template "variant_effect_predictor.sh"
}


// additional variant annotations (sequence context, indel length, etc.)
process annotate_variants {
    executor "local"

    input:
        tuple path(variants), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        path variant_annotations

    shell:
        java_mem = javaMemMB(task)
        variant_annotations = "variant_annotations.txt"
        template "annotate_variants.sh"
}


// -----------------------------------------------------------------------------
// workflow
// -----------------------------------------------------------------------------

workflow {

    sample_sheet = channel.fromPath(params.sampleSheet, checkIfExists: true)
    amplicon_details = channel.fromPath(params.ampliconDetails, checkIfExists: true)
    specific_variants = channel.fromPath(params.specificVariants, checkIfExists: true)
    blacklisted_variants = channel.fromPath(params.blacklistedVariants, checkIfExists: true)

    reference_sequence_fasta_file = file("${params.referenceGenomeFasta}")
    reference_sequence_fasta = channel.fromPath(params.referenceGenomeFasta, checkIfExists: true)
    reference_sequence_index = channel.fromPath("${params.referenceGenomeFasta}.fai", checkIfExists: true)
    reference_sequence_dictionary = channel.fromPath("${params.referenceGenomeFasta}".replaceFirst("${reference_sequence_fasta_file.extension}\$", "dict"), checkIfExists: true)
    reference_sequence = reference_sequence_fasta
        .combine(reference_sequence_index)
        .combine(reference_sequence_dictionary)

    vep_cache_dir = channel.fromPath(params.vepCacheDir, checkIfExists: true)

    // create groups of non-overlapping amplicons
    amplicon_groups = create_non_overlapping_amplicon_groups(
        amplicon_details,
        reference_sequence_index
    )

    // check sample sheet
    check_inputs(sample_sheet, specific_variants, blacklisted_variants, amplicon_groups)

    samples = check_inputs.out.samples
    specific_variants = check_inputs.out.specific_variants
    blacklisted_variants = check_inputs.out.blacklisted_variants

    // BAM file channel created by reading the validated sample sheet
    bam = samples
        .splitCsv(header: true, sep: "\t")
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
    called_variants = call_variants.out.variants
        .collectFile(name: "variants.txt", keepHeader: true)

    // combine called variants with known/expected variants for specific calling
    variants = add_specific_variants(samples, called_variants, specific_variants)

    // generate pileup counts
    pileup_counts(extract_amplicon_regions.out.bam.combine(reference_sequence))

    // collect pileup counts for all samples
    collected_pileup_counts = pileup_counts.out
        .collectFile(name: "pileup_counts.txt", keepHeader: true)

    // fit distributions for substitution allele fractions from pileup counts
    // and compute background noise thresholds
    compute_background_noise_thresholds(collected_pileup_counts)

    // annotate variants using Ensembl VEP
    variant_effect_predictor(variants, vep_cache_dir)

    // additional annotations (sequence context, indel length, etc.)
    annotate_variants(variants.combine(reference_sequence))
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
        specificVariants          : ${params.specificVariants}
        blacklistedVariants       : ${params.blacklistedVariants}
        Reference genome sequence : ${params.referenceGenomeFasta}
        VEP cache directory       : ${params.vepCacheDir}
        Species                   : ${params.vepSpecies}
        Assembly                  : ${params.vepAssembly}
        Output directory          : ${params.outputDir}
        Output prefix             : ${params.outputPrefix}
        Minimum allele fraction   : ${params.minimumAlleleFraction}
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
            --specificVariants         CSV/TSV file containing specific (or known) variants that are included in the summary regardless of whether these are called or not (Sample, Chromosome, Position, Ref, Alt columns required)
            --blacklistedVariants      CSV/TSV file containing blacklisted variants that will be filtered (Chromosome, Position, Ref, Alt columns required)
            --referenceGenomeFasta     FASTA file containing the reference genome sequence (must be indexed, i.e. have an accompanying .fai file)
            --vepCacheDir              Directory in which to install Ensembl VEP cache files
            --vepSpecies               The species name, e.g. homo_sapiens
            --vepAssembly              The genome assembly, e.g. GRCh37
            --outputDir                Directory to which output files are written
            --outputPrefix             Prefix for output file names
            --minimumAlleleFraction    Lower allele fraction limit for detection of variants (for variant callers that provide this option only)

        Alternatively, override settings using a configuration file such as the
        following:

        params {
            sampleSheet           = "samples.csv"
            bamDir                = "bam"
            ampliconDetails       = "amplicons.csv"
            referenceGenomeFasta  = "/reference_data/GRCh37.fa"
            vepCacheDir           = "/reference_data/vep_cache"
            vepSpecies            = "homo_sapiens"
            vepAssembly           = "GRCh37"
            outputDir             = "results"
            outputPrefix          = ""
            minimumAlleleFraction = 0.01
        }

        and run as follows:
            nextflow run crukci-bioinformatics/ampliconseq -c ampliconseq.config
    """.stripIndent()
    log.info ""
}
