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
// check parameters
// -----------------------------------------------------------------------------

variant_caller = params.variantCaller.toLowerCase().replaceAll(/[ _\\-]/, "")
def supported_variant_callers = [ "vardict", "haplotypecaller", "mutect2" ]
if (!supported_variant_callers.contains(variant_caller)) {
    exit 1, "Unsupported variant caller - should be one of " + supported_variant_callers
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

    input:
        path amplicon_details
        path reference_sequence_index

    output:
        path amplicon_groups, emit: amplicon_groups
        path amplicon_bed_files, emit: amplicon_bed_files
        path target_bed_files, emit: target_bed_files

    script:
        amplicon_groups = "amplicon_groups.txt"
        amplicon_bed_prefix = "amplicons"
        amplicon_bed_files = "${amplicon_bed_prefix}.*.bed"
        target_bed_prefix = "targets"
        target_bed_files = "${target_bed_prefix}.*.bed"
        """
        create_non_overlapping_amplicon_groups.R \
            --amplicons ${amplicon_details} \
            --reference-sequence-index ${reference_sequence_index} \
            --output ${amplicon_groups} \
            --amplicon-bed-prefix ${amplicon_bed_prefix} \
            --target-bed-prefix ${target_bed_prefix}
        """
}


// check input files are valid and create validated versions in CSV format
// (inputs can be either TSV or CSV) for use in subsequent processes
process check_inputs {

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


// Picard metrics
process picard_metrics {
    tag "${id}"

    memory { 4.GB * task.attempt }
    time { 2.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(id), val(prefix), path(bam), path(amplicon_groups), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        path alignment_metrics, emit: alignment_metrics
        path targeted_pcr_metrics, emit: targeted_pcr_metrics

    shell:
        java_mem = javaMemMB(task)
        alignment_metrics = "${prefix}.alignment_metrics.txt"
        targeted_pcr_metrics = "${prefix}.targeted_pcr_metrics.txt"
        template "picard_metrics.sh"
}


// extract reads that correspond to groups of amplicon to create
// subset BAM files
process extract_amplicon_regions {
    tag "${id}"

    memory { 2.GB * task.attempt }
    time { 2.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(amplicon_group), path(amplicon_bed), path(target_bed), val(id), val(prefix), path(bam)

    output:
        tuple val(amplicon_group), path(amplicon_bed), path(target_bed), val(id), val(prefix), path(amplicon_bam), path(amplicon_bai), emit: bam
        path(amplicon_coverage), emit: coverage

    shell:
        java_mem = javaMemMB(task)
        amplicon_bam = "${prefix}.${amplicon_group}.bam"
        amplicon_bai = "${prefix}.${amplicon_group}.bai"
        amplicon_coverage = "${prefix}.${amplicon_group}.amplicon_coverage.txt"
        template "extract_amplicon_regions.sh"
}


// call variants and annotate with amplicon id
process call_variants {
    tag "${id}"

    memory { 4.GB * task.attempt }
    time { 8.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(amplicon_group), path(amplicon_bed), path(target_bed), val(id), val(prefix), path(amplicon_bam), path(amplicon_bai), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        tuple val(id), val(prefix), path(vcf), emit: vcf

    shell:
        java_mem = javaMemMB(task)
        vcf = "${prefix}.${amplicon_group}.vcf"
        template "call_variants.sh"
}


// merge amplicon group VCF files for each library and convert to tabular format
process collate_variants {

    publishDir "${params.outputDir}/vcf", mode: "copy", pattern: "${prefix}.vcf"

    input:
        tuple val(id), val(prefix), path(amplicon_vcfs), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        tuple val(id), path(vcf), emit: vcf
        path variants, emit: variants

    shell:
        java_mem = javaMemMB(task)
        vcf = "${prefix}.vcf"
        variants = "${prefix}.variants.txt"
        template "collate_variants.sh"
}


// generate allele counts for each position within each amplicon
process pileup_counts {
    tag "${id}"

    memory { 8.GB * task.attempt }
    time { 4.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(amplicon_group), path(amplicon_bed), path(target_bed), val(id), val(prefix), path(amplicon_bam), path(amplicon_bai), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        tuple val(id), val(prefix), path(pileup_counts)

    shell:
        java_mem = javaMemMB(task)
        pileup_counts = "${prefix}.${amplicon_group}.pileup_counts.txt"
        template "pileup_counts.sh"
}


// annotate and sort pileup counts
process annotate_and_sort_pileup_counts {

    input:
        tuple path(pileup_counts), path(samples), path(amplicons)

    output:
        path annotated_and_sorted_pileup_counts

    script:
        prefix = pileup_counts.getBaseName()
        annotated_and_sorted_pileup_counts = "${prefix}.annotated_and_sorted.txt"
        """
        annotate_and_sort_pileup_counts.R \
            --pileup-counts ${pileup_counts} \
            --samples ${samples} \
            --amplicons ${amplicons} \
            --output "${annotated_and_sorted_pileup_counts}" \
        """
}


// collate alignment and amplicon/target coverage metrics
process collate_alignment_coverage_metrics {

    publishDir "${params.outputDir}", mode: "copy", pattern: "${alignment_coverage_metrics}"
    publishDir "${params.outputDir}", mode: "copy", pattern: "${amplicon_coverage_metrics}"
    publishDir "${params.outputDir}/qc", mode: "copy", pattern: "${sorted_amplicon_coverage}"

    input:
        path samples
        path amplicons
        path alignment_metrics
        path targeted_pcr_metrics
        path amplicon_coverage
        path pileup_counts

    output:
        path alignment_coverage_metrics, emit: alignment_coverage_metrics
        path amplicon_coverage_metrics
        path sorted_amplicon_coverage

    script:
        alignment_coverage_metrics = "alignment_coverage_metrics.csv"
        amplicon_coverage_metrics = "amplicon_coverage_metrics.csv"
        sorted_amplicon_coverage = "amplicon_coverage_metrics.txt"
        """
        collate_alignment_coverage_metrics.R \
            --samples ${samples} \
            --amplicons ${amplicons} \
            --alignment-metrics ${alignment_metrics} \
            --targeted-pcr-metrics ${targeted_pcr_metrics} \
            --amplicon-coverage ${amplicon_coverage} \
            --pileup-counts ${pileup_counts} \
            --alignment-coverage-metrics ${alignment_coverage_metrics} \
            --sorted-amplicon-coverage ${sorted_amplicon_coverage} \
            --amplicon-coverage-metrics ${amplicon_coverage_metrics}
        """
}


// create coverage plots including on/off target/amplicon yield stacked bar plot
// and amplicon coverage box plot
process create_coverage_plots {

    publishDir "${params.outputDir}/qc", mode: "copy"

    input:
        path alignment_coverage_metrics
        path amplicon_coverage

    output:
        path yield_plot_svg, emit: yield_plot
        path yield_plot_pdf
        path amplicon_coverage_plot_svg, emit: amplicon_coverage_plot
        path amplicon_coverage_plot_pdf

    script:
        yield_plot_svg = "yield.svg"
        yield_plot_pdf = "yield.pdf"
        amplicon_coverage_plot_svg = "amplicon_coverage.svg"
        amplicon_coverage_plot_pdf = "amplicon_coverage.pdf"
        """
        create_coverage_plots.R \
            --alignment-metrics ${alignment_coverage_metrics} \
            --amplicon-coverage ${amplicon_coverage}
        """
}


// assess sample replicates based on correlation of SNV allele fractions
process assess_replicate_vaf {

    publishDir "${params.outputDir}/qc", mode: "copy"

    input:
        path samples
        path pileup_counts

    output:
        path vaf_table
        path vaf_heatmap_png, emit: vaf_heatmap
        path vaf_heatmap_svg
        path vaf_heatmap_pdf
        path vaf_correlation_heatmap_png, emit: vaf_correlation_heatmap
        path vaf_correlation_heatmap_svg
        path vaf_correlation_heatmap_pdf
        path mismatched_replicates, emit: mismatched_replicates

    script:
        vaf_table = "allele_fractions.txt"
        vaf_heatmap_png = "vaf_heatmap.png"
        vaf_heatmap_svg = "vaf_heatmap.svg"
        vaf_heatmap_pdf = "vaf_heatmap.pdf"
        vaf_correlation_heatmap_png = "vaf_correlation_heatmap.png"
        vaf_correlation_heatmap_svg = "vaf_correlation_heatmap.svg"
        vaf_correlation_heatmap_pdf = "vaf_correlation_heatmap.pdf"
        mismatched_replicates = "vaf_mismatched_replicates.txt"
        """
        assess_replicate_vaf.R \
            --samples ${samples} \
            --pileup-counts ${pileup_counts}
        """
}


// create QC report from collated alignment/coverage metrics and plots and the
// replicate library allele fraction correlation/clustering
process create_qc_report {

    publishDir "${params.outputDir}", mode: "copy"

    input:
        path alignment_metrics
        path yield_plot
        path amplicon_coverage_plot
        path vaf_heatmap
        path vaf_correlation_heatmap
        path mismatched_replicates

    output:
        path qc_report

    script:
        qc_report = "ampliconseq_qc_report.html"
        """
        create_qc_report.R \
            --alignment-metrics ${alignment_metrics} \
            --yield-plot ${yield_plot} \
            --amplicon-coverage-plot ${amplicon_coverage_plot} \
            --vaf-heatmap ${vaf_heatmap} \
            --vaf-correlation-heatmap ${vaf_correlation_heatmap} \
            --replicate-mismatches ${mismatched_replicates} \
            --output-report ${qc_report}
        """
}


// expand variant table to include specific variants and missing calls within
// sample replicates
process add_specific_variants {

    input:
        path samples
        path called_variants
        path specific_variants
        path reference_sequence_index

    output:
        path all_variants

    script:
        all_variants = "all_variants.txt"
        """
        add_specific_variants.R \
            --samples ${samples} \
            --called-variants ${called_variants} \
            --specific-variants ${specific_variants} \
            --reference-sequence-index ${reference_sequence_index} \
            --output ${all_variants}
        """
}


// add depth and allele fraction from pileup counts to variants
process add_pileup_allele_fractions {

    memory { 8.GB * task.attempt }
    time { 2.hour * task.attempt }
    maxRetries 1

    input:
        path variants
        path pileup_counts

    output:
        path variants_with_pileup_af

    script:
        variants_with_pileup_af = "variants_with_pileup_af.txt"
        """
        add_pileup_allele_fractions.R \
            --variants ${variants} \
            --pileup-counts ${pileup_counts} \
            --output ${variants_with_pileup_af}
        """
}


// fit distributions for substitution allele fractions from pileup counts and
// compute background noise thresholds
process compute_background_noise_thresholds {

    time { 4.hour * task.attempt }
    maxRetries 1

    publishDir "${params.outputDir}", mode: "copy"

    input:
        path pileup_counts
        path variants

    output:
        path position_noise_thresholds, emit: position_noise_thresholds
        path library_noise_thresholds, emit: library_noise_thresholds

    script:
        position_noise_thresholds = "position_noise_thresholds.txt"
        library_noise_thresholds = "library_noise_thresholds.txt"
        """
        compute_background_noise_thresholds.R \
            --pileup-counts ${pileup_counts} \
            --variants ${variants} \
            --position-thresholds ${position_noise_thresholds} \
            --library-thresholds ${library_noise_thresholds} \
            --minimum-depth ${params.minimumDepthForBackgroundNoise} \
            --exclude-highest-fraction ${params.excludeHighestFractionForBackgroundNoise} \
            --maximum-allele-fraction ${params.maximumAlleleFractionForBackgroundNoise} \
            --minimum-number-for-fitting ${params.minimumNumberForFittingBackgroundNoise} \
            --chunk-size ${params.chunkSizeForFittingBackgroundNoise} \
            --read-chunk-size ${params.readChunkSizeForFittingBackgroundNoise}
        """
}


// add background noise thresholds to variants table and applies background
// noise filters
process apply_background_noise_filters {

    input:
        path variants
        path position_noise_thresholds
        path library_noise_thresholds

    output:
        path filtered_variants

    script:
        filtered_variants = "variants_background_noise_filtered.txt"
        """
        apply_background_noise_filters.R \
            --variants ${variants} \
            --position-thresholds ${position_noise_thresholds} \
            --library-thresholds ${library_noise_thresholds} \
            --output ${filtered_variants}
        """
}


// annotate variants using Ensembl VEP
process variant_effect_predictor {

    memory 4.GB

    input:
        path variants
        path reference_sequence_index
        path vep_cache_dir
        val one_annotation_per_variant

    output:
        path vep_annotations

    shell:
        vep_annotations = "vep_annotations.txt"
        vep_pick_option = one_annotation_per_variant ? "--pick" : ""
        template "variant_effect_predictor.sh"
}


// additional variant annotations (sequence context, indel length, etc.)
process annotate_variants {

    input:
        tuple path(variants), path(amplicons), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)

    output:
        path offset_from_primer_end_annotations, emit: offset_from_primer_end_annotations
        path other_annotations, emit: other_annotations

    shell:
        java_mem = javaMemMB(task)
        offset_from_primer_end_annotations = "offset_from_primer_end_annotations.txt"
        other_annotations = "other_annotations.txt"
        template "annotate_variants.sh"
}


// gather variant calls/details for replicate libraries into a single row and
// add VEP and additional annotations
process summarize_variants {

    publishDir "${params.outputDir}", mode: "copy"

    input:
        path variants
        path blacklisted_variants
        path vep_annotations
        path offset_from_primer_end_annotations
        path other_annotations
        path reference_sequence_index

    output:
        path variant_summary_csv
        path variant_summary_tsv

    script:
        output_prefix = "variants"
        variant_summary_csv = "${output_prefix}.csv"
        variant_summary_tsv = "${output_prefix}.txt"
        vep_annotations_option = vep_annotations.name == "NO_FILE" ? "" : "--vep-annotations ${vep_annotations}"
        """
        summarize_variants.R \
            --variants ${variants} \
            --blacklist ${blacklisted_variants} \
            ${vep_annotations_option} \
            --offset-from-primer-end-annotations ${offset_from_primer_end_annotations} \
            --other-annotations ${other_annotations} \
            --reference-sequence-index ${reference_sequence_index} \
            --output-prefix ${output_prefix} \
            --minimum-depth ${params.minimumDepthForHighConfidenceCalls} \
            --minimum-alt-depth ${params.minimumAltDepthForHighConfidenceCalls}
        """
}


// -----------------------------------------------------------------------------
// workflow
// -----------------------------------------------------------------------------

workflow {

    sample_sheet = channel.fromPath(params.samples, checkIfExists: true)
    amplicon_details = channel.fromPath(params.amplicons, checkIfExists: true)
    specific_variants = channel.fromPath(params.specificVariants, checkIfExists: true)
    blacklisted_variants = channel.fromPath(params.blacklistedVariants, checkIfExists: true)

    reference_sequence_fasta_file = file("${params.referenceGenomeFasta}")
    reference_sequence_fasta = channel.fromPath(params.referenceGenomeFasta, checkIfExists: true)
    reference_sequence_index = channel.fromPath("${params.referenceGenomeFasta}.fai", checkIfExists: true)
    reference_sequence_dictionary = channel.fromPath("${params.referenceGenomeFasta}".replaceFirst("${reference_sequence_fasta_file.extension}\$", "dict"), checkIfExists: true)
    reference_sequence = reference_sequence_fasta
        .combine(reference_sequence_index)
        .combine(reference_sequence_dictionary)

    vep_cache_dir = (params.vepAnnotation ? channel.fromPath(params.vepCacheDir, checkIfExists: true) : Channel.empty())

    // create groups of non-overlapping amplicons
    create_non_overlapping_amplicon_groups(amplicon_details, reference_sequence_index)
    amplicon_groups = create_non_overlapping_amplicon_groups.out.amplicon_groups

    // check sample sheet
    check_inputs(sample_sheet, specific_variants, blacklisted_variants, amplicon_groups)

    samples = check_inputs.out.samples
    specific_variants = check_inputs.out.specific_variants
    blacklisted_variants = check_inputs.out.blacklisted_variants

    // BAM file channel created by reading the validated sample sheet
    bams = samples
        .splitCsv(header: true, sep: "\t")
        .map { row -> tuple("${row.ID}", "${row.ID}".replaceFirst(/ /, "_"), file(!params.bamDir ? "${row.BAM}" : "${params.bamDir}/${row.BAM}", checkIfExists: true)) }

    amplicon_bed_files = create_non_overlapping_amplicon_groups.out.amplicon_bed_files
        .flatten()
        .map { tuple((it =~ /.*\.(\d+)\.bed$/)[0][1], it) }

    target_bed_files = create_non_overlapping_amplicon_groups.out.target_bed_files
        .flatten()
        .map { tuple((it =~ /.*\.(\d+)\.bed$/)[0][1], it) }

    bed_files = amplicon_bed_files.join(target_bed_files)

    extract_amplicon_regions(bed_files.combine(bams))

    // collect amplicon coverage data for all samples
    amplicon_coverage = extract_amplicon_regions.out.coverage
        .collectFile(name: "amplicon_coverage.txt", keepHeader: true, sort: { it.name })

    amplicon_bams = extract_amplicon_regions.out.bam.combine(reference_sequence)

    // Picard alignment summary metrics
    picard_metrics(bams.combine(amplicon_groups).combine(reference_sequence))

    // collect Picard metrics for all samples
    alignment_metrics = picard_metrics.out.alignment_metrics
        .collectFile(name: "alignment_metrics.txt", keepHeader: true, sort: { it.name }, storeDir: "${params.outputDir}/qc")
    targeted_pcr_metrics = picard_metrics.out.targeted_pcr_metrics
        .collectFile(name: "targeted_pcr_metrics.txt", keepHeader: true, sort: { it.name }, storeDir: "${params.outputDir}/qc")

    // generate pileup counts
    pileup_counts(amplicon_bams)

    // collate pileup counts for each library
    collected_pileup_counts = pileup_counts.out.collectFile(keepHeader: true) { item -> [ "${item[1]}.pileup_counts.txt", item[2] ] }

    // annotate and sort pileup counts
    annotate_and_sort_pileup_counts(collected_pileup_counts.combine(samples).combine(amplicon_groups))

    // collect pileup counts for all libraries
    collected_pileup_counts = annotate_and_sort_pileup_counts.out
        .collectFile(name: "pileup_counts.txt", keepHeader: true, sort: { it.name }, storeDir: "${params.outputDir}")

    // call variants
    call_variants(amplicon_bams)

    // merge amplicon group VCF files for each library and convert to tabular format
    collate_variants(call_variants.out.groupTuple(by: [0, 1]).combine(reference_sequence))

    // collect variant calls for all samples
    called_variants = collate_variants.out.variants
        .collectFile(name: "variants.txt", keepHeader: true)

    // collate alignment and target coverage metrics
    collate_alignment_coverage_metrics(
        samples,
        amplicon_groups,
        alignment_metrics,
        targeted_pcr_metrics,
        amplicon_coverage,
        collected_pileup_counts
    )

    // create coverage plots including on/off target/amplicon yield stacked bar plot
    // and amplicon coverage box plot
    create_coverage_plots(collate_alignment_coverage_metrics.out.alignment_coverage_metrics, amplicon_coverage)

    // assess sample replicates based on correlation of SNV allele fractions
    assess_replicate_vaf(samples, collected_pileup_counts)

    // create QC report
    create_qc_report(
        collate_alignment_coverage_metrics.out.alignment_coverage_metrics,
        create_coverage_plots.out.yield_plot,
        create_coverage_plots.out.amplicon_coverage_plot,
        assess_replicate_vaf.out.vaf_heatmap,
        assess_replicate_vaf.out.vaf_correlation_heatmap,
        assess_replicate_vaf.out.mismatched_replicates,
    )

    // combine called variants with known/expected variants for specific calling
    add_specific_variants(samples, called_variants, specific_variants, reference_sequence_index)

    // add depth and allele fraction from pileup counts to variants
    add_pileup_allele_fractions(add_specific_variants.out, collected_pileup_counts)

    // annotate variants using Ensembl VEP
    one_annotation_per_variant = Channel.value(params.vepPickOneAnnotationPerVariant)
    variant_effect_predictor(
        add_specific_variants.out,
        reference_sequence_index,
        vep_cache_dir,
        one_annotation_per_variant
    )

    vep_annotations = ( params.vepAnnotation ? variant_effect_predictor.out : Channel.fromPath("NO_FILE") )

    // additional annotations (sequence context, indel length, etc.)
    annotate_variants(add_specific_variants.out.combine(amplicon_groups).combine(reference_sequence))

    // fit distributions for substitution allele fractions from pileup counts
    // and compute background noise thresholds
    compute_background_noise_thresholds(collected_pileup_counts, add_specific_variants.out)

    // apply background noise filters
    apply_background_noise_filters(
        add_pileup_allele_fractions.out,
        compute_background_noise_thresholds.out.position_noise_thresholds,
        compute_background_noise_thresholds.out.library_noise_thresholds
    )

    // create variant summary
    summarize_variants(
        apply_background_noise_filters.out,
        blacklisted_variants,
        vep_annotations,
        annotate_variants.out.offset_from_primer_end_annotations,
        annotate_variants.out.other_annotations,
        reference_sequence_index
    )
}


// -----------------------------------------------------------------------------
// summary of configuration parameters
// -----------------------------------------------------------------------------

def printParameterSummary() {
    log.info ""
    log.info """
        Variant calling pipeline for amplicon sequencing data
        =====================================================

        Sample sheet               : ${params.samples}
        BAM directory              : ${params.bamDir}
        Amplicon details           : ${params.amplicons}
        specificVariants           : ${params.specificVariants}
        blacklistedVariants        : ${params.blacklistedVariants}
        Reference genome sequence  : ${params.referenceGenomeFasta}
        VEP annotation             : ${params.vepAnnotation}
        VEP cache directory        : ${params.vepCacheDir}
        Species                    : ${params.vepSpecies}
        Assembly                   : ${params.vepAssembly}
        One annotation per variant : ${params.vepPickOneAnnotationPerVariant}
        Output directory           : ${params.outputDir}
        Variant caller             : ${params.variantCaller}
        Minimum allele fraction    : ${params.minimumAlleleFraction}
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
            --samples                  CSV/TSV file containing details of samples and libraries (ID and Sample columns required)
            --bamDir                   Directory in which BAM files are located
            --amplicons                CSV/TSV file containing amplicon coordinates (ID, Chromosome, AmpliconStart, AmpliconEnd, TargetStart, TargetEnd columns required)
            --specificVariants         CSV/TSV file containing specific (or known) variants that are included in the summary regardless of whether these are called or not (Sample, Chromosome, Position, Ref, Alt columns required)
            --blacklistedVariants      CSV/TSV file containing blacklisted variants that will be filtered (Chromosome, Position, Ref, Alt columns required)
            --referenceGenomeFasta     FASTA file containing the reference genome sequence (must be indexed and have an accompanying sequence dictionary)
            --vepAnnotation            Annotate variants with Ensembl VEP
            --vepCacheDir              Directory in which Ensembl VEP cache files are installed
            --vepSpecies               The species name, e.g. homo_sapiens
            --vepAssembly              The genome assembly, e.g. GRCh37
            --outputDir                Directory to which output files are written
            --variantCaller            The variant caller (VarDict, HaplotypeCaller or Mutect2)
            --minimumAlleleFraction    Lower allele fraction limit for detection of variants (for variant callers that provide this option only)

        Alternatively, override settings using a configuration file such as the
        following:

        params {
            samples               = "samples.csv"
            bamDir                = "bam"
            amplicons             = "amplicons.csv"
            referenceGenomeFasta  = "/reference_data/GRCh37.fa"
            vepAnnotation         = true
            vepCacheDir           = "/reference_data/vep_cache"
            vepSpecies            = "homo_sapiens"
            vepAssembly           = "GRCh37"
            outputDir             = "results"
            variantCaller         = "VarDict"
            minimumAlleleFraction = 0.01
        }

        and run as follows:
            nextflow run crukci-bioinformatics/ampliconseq -c ampliconseq.config
    """.stripIndent()
    log.info ""
}

