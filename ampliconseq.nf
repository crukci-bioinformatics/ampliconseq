#!/usr/bin/env nextflow

// enable DSL 2 syntax
nextflow.enable.dsl = 2


// -----------------------------------------------------------------------------
// functions
// -----------------------------------------------------------------------------

// returns the heap size to be given to a Java task based on the task memory,
// allowing for some overhead for the JVM
def javaMemMB(task, jvmOverhead)
{
    return task.memory.toMega() - jvmOverhead
}

// normalized name of the configured variant caller (lower case, no spaces,
// underscores or hyphens)
def normalizedVariantCaller()
{
    return params.variantCaller.toLowerCase().replaceAll(/[ _\\-]/, "")
}

// the set of supported variant callers
def supportedVariantCallers()
{
    return [ "vardict", "haplotypecaller", "mutect2" ]
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

    script:
        checked_sample_sheet = "samples.checked.txt"
        checked_specific_variants = "specific_variants.checked.txt"
        checked_blacklisted_variants = "blacklisted_variants.checked.txt"
        """
        set -e -o pipefail

        check_sample_sheet.R \\
            --input ${sample_sheet} \\
            --output ${checked_sample_sheet}

        check_specific_variants.R \\
            --input ${specific_variants} \\
            --samples ${checked_sample_sheet} \\
            --amplicons ${amplicons} \\
            --output ${checked_specific_variants}

        check_blacklisted_variants.R \\
            --input ${blacklisted_variants} \\
            --output ${checked_blacklisted_variants}
        """
}


// Picard metrics
process picard_metrics {
    tag "${id}"

    memory { 4.GB * task.attempt }
    time { 2.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(id), val(prefix), path(bam), path(amplicon_groups), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)
        val jvm_overhead

    output:
        path alignment_metrics, emit: alignment_metrics
        path targeted_pcr_metrics, emit: targeted_pcr_metrics

    script:
        java_mem = javaMemMB(task, jvm_overhead)
        alignment_metrics = "${prefix}.alignment_metrics.txt"
        targeted_pcr_metrics = "${prefix}.targeted_pcr_metrics.txt"
        """
        set -e -o pipefail

        # Picard CollectAlignmentSummaryMetrics bundled with GATK
        gatk --java-options "-Xmx${java_mem}m" CollectAlignmentSummaryMetrics \\
            --INPUT ${bam} \\
            --REFERENCE_SEQUENCE ${reference_sequence} \\
            --OUTPUT alignment_metrics.txt

        extract_picard_metrics.R \\
            --id "${id}" \\
            --metrics alignment_metrics.txt \\
            --output "${alignment_metrics}"

        # extract amplicon and target intervals in BED format and convert to Picard
        # interval list format
        awk 'BEGIN { FS = "\\t"; OFS = "\\t" } FNR > 1 { print \$2, \$3, \$4, \$1 }' ${amplicon_groups} > amplicons.bed
        gatk --java-options "-Xmx${java_mem}m" BedToIntervalList \\
            --INPUT amplicons.bed \\
            --SEQUENCE_DICTIONARY ${reference_sequence_dictionary} \\
            --OUTPUT amplicons.interval_list.txt

        awk 'BEGIN { FS = "\\t"; OFS = "\\t" } FNR > 1 { print \$2, \$5, \$6, \$1 }' ${amplicon_groups} > targets.bed
        gatk --java-options "-Xmx${java_mem}m" BedToIntervalList \\
            --INPUT targets.bed \\
            --SEQUENCE_DICTIONARY ${reference_sequence_dictionary} \\
            --OUTPUT targets.interval_list.txt

        # Picard CollectTargetedPcrMetrics bundled with GATK
        gatk --java-options "-Xmx${java_mem}m" CollectTargetedPcrMetrics \\
            --INPUT ${bam} \\
            --REFERENCE_SEQUENCE ${reference_sequence} \\
            --AMPLICON_INTERVALS amplicons.interval_list.txt \\
            --TARGET_INTERVALS targets.interval_list.txt \\
            --OUTPUT targeted_pcr_metrics.txt

        extract_picard_metrics.R \\
            --id "${id}" \\
            --metrics targeted_pcr_metrics.txt \\
            --output "${targeted_pcr_metrics}"
        """
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
        val maximum_distance
        val require_both_ends_anchored
        val jvm_overhead

    output:
        tuple val(amplicon_group), path(amplicon_bed), path(target_bed), val(id), val(prefix), path(amplicon_bam), path(amplicon_bai), emit: bam
        path(amplicon_coverage), emit: coverage

    script:
        java_mem = javaMemMB(task, jvm_overhead)
        amplicon_bam = "${prefix}.${amplicon_group}.bam"
        amplicon_bai = "${prefix}.${amplicon_group}.bai"
        amplicon_coverage = "${prefix}.${amplicon_group}.amplicon_coverage.txt"
        """
        set -e -o pipefail

        JAVA_OPTS="-Xmx${java_mem}m" extract-amplicon-regions \\
            --id "${id}" \\
            --input ${bam} \\
            --amplicon-intervals ${amplicon_bed} \\
            --output "${amplicon_bam}" \\
            --coverage "${amplicon_coverage}" \\
            --maximum-distance ${maximum_distance} \\
            --require-both-ends-anchored=${require_both_ends_anchored} \\
            --unmark-duplicate-reads
        """
}


// call variants and annotate with amplicon id
process call_variants {
    tag "${id}"

    memory { 4.GB * task.attempt }
    time { 8.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(amplicon_group), path(amplicon_bed), path(target_bed), val(id), val(prefix), path(amplicon_bam), path(amplicon_bai), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)
        val minimum_allele_fraction
        val maximum_reads_per_alignment_start
        val jvm_overhead

    output:
        tuple val(id), val(prefix), path(vcf), emit: vcf

    script:
        java_mem = javaMemMB(task, jvm_overhead)
        variant_caller = normalizedVariantCaller()
        vcf = "${prefix}.${amplicon_group}.vcf"
        """
        set -e -o pipefail

        if [ "${variant_caller}" == "vardict" ]; then

            JAVA_OPTS="-Xmx${java_mem}m -XX:-UsePerfData" vardict-java \\
                -b ${amplicon_bam} \\
                -G ${reference_sequence} \\
                -N "${id}" \\
                -f ${minimum_allele_fraction} \\
                -z -c 1 -S 2 -E 3 -g 4 ${target_bed} \\
                > vardict.txt

            cat vardict.txt \\
                | teststrandbias.R \\
                > vardict.teststrandbias.txt

            cat vardict.teststrandbias.txt \\
                | var2vcf_valid.pl -N "${id}" -E -P 0 -f ${minimum_allele_fraction} \\
                > variants.vcf

        elif [ "${variant_caller}" == "haplotypecaller" ]; then

            gatk --java-options "-Xmx${java_mem}m" HaplotypeCaller \\
                --input ${amplicon_bam} \\
                --intervals ${target_bed} \\
                --reference ${reference_sequence} \\
                --output haplotypecaller.vcf \\
                --max-reads-per-alignment-start ${maximum_reads_per_alignment_start} \\
                --native-pair-hmm-threads 1 \\
                --force-active

            gatk --java-options "-Xmx${java_mem}m" VariantFiltration \\
                --variant haplotypecaller.vcf \\
                --reference ${reference_sequence} \\
                --filter-name "QualByDepth" --filter-expression "QD < 2.0" \\
                --filter-name "StrandOddsRatio" --filter-expression "SOR > 3.0" \\
                --filter-name "RMSMappingQuality" --filter-expression "MQ < 40.0" \\
                --filter-name "MappingQualityRankSumTest" --filter-expression "MQRankSum < -12.5" \\
                --output variants.vcf

        elif [ "${variant_caller}" == "mutect2" ]; then

            gatk --java-options "-Xmx${java_mem}m" Mutect2 \\
                --input ${amplicon_bam} \\
                --intervals ${target_bed} \\
                --reference ${reference_sequence} \\
                --output mutect.vcf \\
                --max-reads-per-alignment-start ${maximum_reads_per_alignment_start} \\
                --minimum-allele-fraction ${minimum_allele_fraction} \\
                --native-pair-hmm-threads 1 \\
                --force-active

            gatk --java-options "-Xmx${java_mem}m" FilterMutectCalls \\
                --variant mutect.vcf \\
                --reference ${reference_sequence} \\
                --output variants.vcf \\
                --min-allele-fraction ${minimum_allele_fraction}

        else
            echo "Unrecognized variant caller: ${variant_caller}" >&2
            exit 1
        fi

        JAVA_OPTS="-Xmx${java_mem}m" annotate-vcf-with-amplicon-ids \\
            --input variants.vcf \\
            --target-intervals ${target_bed} \\
            --output "${vcf}"
        """
}


// merge amplicon group VCF files for each library and convert to tabular format
process collate_variants {

    input:
        tuple val(id), val(prefix), path(amplicon_vcfs), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)
        val jvm_overhead

    output:
        tuple val(id), path(vcf), emit: vcf
        path variants, emit: variants

    script:
        java_mem = javaMemMB(task, jvm_overhead)
        vcf = "${prefix}.vcf"
        variants = "${prefix}.variants.txt"
        """
        set -e -o pipefail

        for amplicon_vcf in ${amplicon_vcfs}
        do
            echo \${amplicon_vcf} >> vcf_list.txt
        done

        gatk --java-options "-Xmx${java_mem}m" MergeVcfs \\
            --INPUT vcf_list.txt \\
            --SEQUENCE_DICTIONARY ${reference_sequence_dictionary} \\
            --OUTPUT "${vcf}"

        add-assorted-annotations-to-vcf \\
            --input "${vcf}" \\
            --output merged_annotated.vcf \\
            --reference-sequence ${reference_sequence} \\
            --sequence-context-length 1

        # left-align and normalize indels and split multi-allelic variants
        # note that this can still result in deletions that have an '*' for the
        # alt allele which will need fixing
        bcftools norm \\
            --multiallelics -both \\
            --fasta-ref ${reference_sequence} \\
            merged_annotated.vcf \\
            --output merged_annotated_normalized.vcf

        gatk --java-options "-Xmx${java_mem}m" VariantsToTable \\
            --variant merged_annotated_normalized.vcf \\
            --output merged_annotated_normalized.txt \\
            --show-filtered \\
            --fields AMPLICON \\
            --fields CHROM \\
            --fields POS \\
            --fields REF \\
            --fields ALT \\
            --fields MULTIALLELIC \\
            --fields FILTER \\
            --fields QUAL \\
            --fields FivePrimeContext \\
            --genotype-fields DP \\
            --asGenotypeFieldsToTake AD \\
            --asGenotypeFieldsToTake AF

        tidy_variant_table.R \\
            --input merged_annotated_normalized.txt \\
            --id "${id}" \\
            --output "${variants}"
        """
}


// generate allele counts for each position within each amplicon
process pileup_counts {
    tag "${id}"

    memory { 8.GB * task.attempt }
    time { 4.hour * task.attempt }
    maxRetries 2

    input:
        tuple val(amplicon_group), path(amplicon_bed), path(target_bed), val(id), val(prefix), path(amplicon_bam), path(amplicon_bai), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)
        val minimum_mapping_quality
        val minimum_base_quality
        val jvm_overhead

    output:
        tuple val(id), val(prefix), path(pileup_counts)

    script:
        java_mem = javaMemMB(task, jvm_overhead)
        pileup_counts = "${prefix}.${amplicon_group}.pileup_counts.txt"
        """
        JAVA_OPTS="-Xmx${java_mem}m" pileup-counts \\
            --id "${id}" \\
            --input ${amplicon_bam} \\
            --amplicon-intervals ${target_bed} \\
            --reference-sequence ${reference_sequence} \\
            --output "${pileup_counts}" \\
            --minimum-mapping-quality ${minimum_mapping_quality} \\
            --minimum-base-quality ${minimum_base_quality}
        """
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

    input:
        path samples
        path amplicons
        path alignment_metrics
        path targeted_pcr_metrics
        path amplicon_coverage
        path pileup_counts

    output:
        path alignment_coverage_metrics, emit: alignment_coverage_metrics
        path amplicon_coverage_metrics, emit: amplicon_coverage_metrics
        path sorted_amplicon_coverage, emit: sorted_amplicon_coverage

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

    input:
        path alignment_coverage_metrics
        path amplicon_coverage

    output:
        path yield_plot_svg, emit: yield_plot
        path yield_plot_pdf, emit: yield_plot_pdf
        path amplicon_coverage_plot_svg, emit: amplicon_coverage_plot
        path amplicon_coverage_plot_pdf, emit: amplicon_coverage_plot_pdf

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

    input:
        path samples
        path pileup_counts

    output:
        path vaf_table, emit: vaf_table
        path vaf_heatmap_png, emit: vaf_heatmap
        path vaf_heatmap_svg, emit: vaf_heatmap_svg
        path vaf_heatmap_pdf, emit: vaf_heatmap_pdf
        path vaf_correlation_heatmap_png, emit: vaf_correlation_heatmap
        path vaf_correlation_heatmap_svg, emit: vaf_correlation_heatmap_svg
        path vaf_correlation_heatmap_pdf, emit: vaf_correlation_heatmap_pdf
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

    input:
        path alignment_metrics
        path yield_plot
        path amplicon_coverage_plot
        path vaf_heatmap
        path vaf_correlation_heatmap
        path mismatched_replicates

    output:
        path qc_report, emit: qc_report

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

    input:
        path pileup_counts
        path variants
        val minimum_depth
        val exclude_highest_fraction
        val maximum_allele_fraction
        val minimum_number_for_fitting
        val chunk_size
        val read_chunk_size

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
            --minimum-depth ${minimum_depth} \
            --exclude-highest-fraction ${exclude_highest_fraction} \
            --maximum-allele-fraction ${maximum_allele_fraction} \
            --minimum-number-for-fitting ${minimum_number_for_fitting} \
            --chunk-size ${chunk_size} \
            --read-chunk-size ${read_chunk_size}
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
        val vep_species
        val vep_assembly

    output:
        path vep_annotations

    script:
        vep_annotations = "vep_annotations.txt"
        vep_pick_option = one_annotation_per_variant ? "--pick" : ""
        """
        set -e -o pipefail

        create_distinct_vcf.R \\
            --input ${variants} \\
            --reference-sequence-index ${reference_sequence_index} \\
            --output distinct_variants.vcf

        vep \\
            --input_file distinct_variants.vcf \\
            --format vcf \\
            --output_file distinct_variants.vep.vcf \\
            --vcf \\
            --offline \\
            --dir_cache ${vep_cache_dir} \\
            --species ${vep_species} \\
            --assembly ${vep_assembly} \\
            --buffer_size 100 \\
            --no_stats \\
            --dont_skip \\
            ${vep_pick_option} \\
            --everything \\
            --no_escape

        vep_vcf_to_tabular.R \\
            --vcf distinct_variants.vep.vcf \\
            --output ${vep_annotations}
        """
}


// additional variant annotations (sequence context, indel length, etc.)
process annotate_variants {

    input:
        tuple path(variants), path(amplicons), path(reference_sequence), path(reference_sequence_index), path(reference_sequence_dictionary)
        val sequence_context_length
        val jvm_overhead

    output:
        path offset_from_primer_end_annotations, emit: offset_from_primer_end_annotations
        path other_annotations, emit: other_annotations

    script:
        java_mem = javaMemMB(task, jvm_overhead)
        offset_from_primer_end_annotations = "offset_from_primer_end_annotations.txt"
        other_annotations = "other_annotations.txt"
        """
        set -e -o pipefail

        create_distinct_vcf.R \\
            --input ${variants} \\
            --reference-sequence-index ${reference_sequence_index} \\
            --output distinct_variants.vcf

        add-assorted-annotations-to-vcf \\
            --input distinct_variants.vcf \\
            --output distinct_variants.annotated.vcf \\
            --reference-sequence ${reference_sequence} \\
            --sequence-context-length ${sequence_context_length}

        gatk --java-options "-Xmx${java_mem}m" VariantsToTable \\
            --variant distinct_variants.annotated.vcf \\
            --output distinct_variants.annotated.txt \\
            --fields CHROM \\
            --fields POS \\
            --fields REF \\
            --fields ALT \\
            --fields FivePrimeContext \\
            --fields ThreePrimeContext \\
            --fields IndelLength

        echo -e "Chromosome\\tPosition\\tRef\\tAlt\\t5' context\\tAlleles\\t3' context\\tIndel length" > ${other_annotations}
        awk 'BEGIN { FS = "\\t"; OFS = "\\t" } FNR > 1 { print \$1, \$2, \$3, \$4, \$5, \$3"/"\$4, \$6, \$7 }' distinct_variants.annotated.txt >> ${other_annotations}

        add_offset_from_primer_end.R \\
            --variants ${variants} \\
            --amplicons ${amplicons} \\
            --output ${offset_from_primer_end_annotations}
        """
}


// gather variant calls/details for replicate libraries into a single row and
// add VEP and additional annotations
process summarize_variants {

    input:
        path variants
        path blacklisted_variants
        path vep_annotations
        path offset_from_primer_end_annotations
        path other_annotations
        path reference_sequence_index
        val minimum_depth
        val minimum_alt_depth

    output:
        path variant_summary_csv, emit: variant_summary_csv
        path variant_summary_tsv, emit: variant_summary_tsv

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
            --minimum-depth ${minimum_depth} \
            --minimum-alt-depth ${minimum_alt_depth}
        """
}


// -----------------------------------------------------------------------------
// workflow
// -----------------------------------------------------------------------------

workflow {

    main:

    // show settings and/or help
    printParameterSummary()

    if (params.help) {
        helpMessage()
        System.exit(0)
    }

    // check parameters
    if (!supportedVariantCallers().contains(normalizedVariantCaller())) {
        error "Unsupported variant caller - should be one of " + supportedVariantCallers()
    }

    if (params.containsKey('outputDir')) {
        error "The outputDir parameter is no longer supported - set the output directory using the -output-dir command line option or the outputDir setting in a configuration file (outside the params block)."
    }

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

    vep_cache_dir = (params.vepAnnotation ? channel.fromPath(params.vepCacheDir, checkIfExists: true) : channel.empty())

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
        .map { bed -> tuple((bed =~ /.*\.(\d+)\.bed$/)[0][1], bed) }

    target_bed_files = create_non_overlapping_amplicon_groups.out.target_bed_files
        .flatten()
        .map { bed -> tuple((bed =~ /.*\.(\d+)\.bed$/)[0][1], bed) }

    bed_files = amplicon_bed_files.join(target_bed_files)

    extract_amplicon_regions(bed_files.combine(bams), params.maxDistanceFromAmpliconEnd, params.requireBothEndsAnchored, params.jvmOverhead)

    // collect amplicon coverage data for all samples
    amplicon_coverage = extract_amplicon_regions.out.coverage
        .collectFile(name: "amplicon_coverage.txt", keepHeader: true, sort: { file -> file.name })

    amplicon_bams = extract_amplicon_regions.out.bam.combine(reference_sequence)

    // Picard alignment summary metrics
    picard_metrics(bams.combine(amplicon_groups).combine(reference_sequence), params.jvmOverhead)

    // collect Picard metrics for all samples
    alignment_metrics = picard_metrics.out.alignment_metrics
        .collectFile(name: "alignment_metrics.txt", keepHeader: true, sort: { file -> file.name })
    targeted_pcr_metrics = picard_metrics.out.targeted_pcr_metrics
        .collectFile(name: "targeted_pcr_metrics.txt", keepHeader: true, sort: { file -> file.name })

    // generate pileup counts
    pileup_counts(amplicon_bams, params.minimumMappingQualityForPileup, params.minimumBaseQualityForPileup, params.jvmOverhead)

    // collate pileup counts for each library
    collected_pileup_counts = pileup_counts.out.collectFile(keepHeader: true) { item -> [ "${item[1]}.pileup_counts.txt", item[2] ] }

    // annotate and sort pileup counts
    annotate_and_sort_pileup_counts(collected_pileup_counts.combine(samples).combine(amplicon_groups))

    // collect pileup counts for all libraries
    collected_pileup_counts = annotate_and_sort_pileup_counts.out
        .collectFile(name: "pileup_counts.txt", keepHeader: true, sort: { file -> file.name })

    // call variants
    call_variants(amplicon_bams, params.minimumAlleleFraction, params.maximumReadsPerAlignmentStart, params.jvmOverhead)

    // merge amplicon group VCF files for each library and convert to tabular format
    collate_variants(call_variants.out.groupTuple(by: [0, 1]).combine(reference_sequence), params.jvmOverhead)

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
    variant_effect_predictor(
        add_specific_variants.out,
        reference_sequence_index,
        vep_cache_dir,
        params.vepPickOneAnnotationPerVariant,
        params.vepSpecies,
        params.vepAssembly
    )

    vep_annotations = ( params.vepAnnotation ? variant_effect_predictor.out : channel.fromPath("NO_FILE") )

    // additional annotations (sequence context, indel length, etc.)
    annotate_variants(add_specific_variants.out.combine(amplicon_groups).combine(reference_sequence), params.sequenceContextLength, params.jvmOverhead)

    // fit distributions for substitution allele fractions from pileup counts
    // and compute background noise thresholds
    compute_background_noise_thresholds(
        collected_pileup_counts,
        add_specific_variants.out,
        params.minimumDepthForBackgroundNoise,
        params.excludeHighestFractionForBackgroundNoise,
        params.maximumAlleleFractionForBackgroundNoise,
        params.minimumNumberForFittingBackgroundNoise,
        params.chunkSizeForFittingBackgroundNoise,
        params.readChunkSizeForFittingBackgroundNoise
    )

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
        reference_sequence_index,
        params.minimumDepthForHighConfidenceCalls,
        params.minimumAltDepthForHighConfidenceCalls
    )

    publish:

    // per-sample merged VCF files
    variant_vcfs = collate_variants.out.vcf

    // alignment and coverage metrics summary tables
    coverage_metrics = collate_alignment_coverage_metrics.out.alignment_coverage_metrics
        .mix(collate_alignment_coverage_metrics.out.amplicon_coverage_metrics)
    sorted_amplicon_coverage = collate_alignment_coverage_metrics.out.sorted_amplicon_coverage

    // aggregated per-sample metrics and pileup counts
    alignment_metrics_table = alignment_metrics
    targeted_pcr_metrics_table = targeted_pcr_metrics
    pileup_counts_table = collected_pileup_counts

    // coverage plots
    coverage_plots = create_coverage_plots.out.yield_plot
        .mix(
            create_coverage_plots.out.yield_plot_pdf,
            create_coverage_plots.out.amplicon_coverage_plot,
            create_coverage_plots.out.amplicon_coverage_plot_pdf
        )

    // replicate VAF assessment
    replicate_vaf = assess_replicate_vaf.out.vaf_table
        .mix(
            assess_replicate_vaf.out.vaf_heatmap,
            assess_replicate_vaf.out.vaf_heatmap_svg,
            assess_replicate_vaf.out.vaf_heatmap_pdf,
            assess_replicate_vaf.out.vaf_correlation_heatmap,
            assess_replicate_vaf.out.vaf_correlation_heatmap_svg,
            assess_replicate_vaf.out.vaf_correlation_heatmap_pdf,
            assess_replicate_vaf.out.mismatched_replicates
        )

    // QC report
    qc_report = create_qc_report.out.qc_report

    // background noise thresholds
    background_noise_thresholds = compute_background_noise_thresholds.out.position_noise_thresholds
        .mix(compute_background_noise_thresholds.out.library_noise_thresholds)

    // variant summary tables
    variant_summary = summarize_variants.out.variant_summary_csv
        .mix(summarize_variants.out.variant_summary_tsv)
}


output {
    variant_vcfs               { path 'vcf' }
    coverage_metrics           { path '.'  }
    sorted_amplicon_coverage   { path 'qc' }
    alignment_metrics_table    { path 'qc' }
    targeted_pcr_metrics_table { path 'qc' }
    pileup_counts_table        { path '.'  }
    coverage_plots             { path 'qc' }
    replicate_vaf              { path 'qc' }
    qc_report                  { path '.'  }
    background_noise_thresholds { path '.'  }
    variant_summary            { path '.'  }
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
        Output directory           : ${workflow.outputDir}
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
            variantCaller         = "VarDict"
            minimumAlleleFraction = 0.01
        }

        // the output directory is set outside the params block (or with -output-dir)
        outputDir = "results"

        and run as follows:
            nextflow run crukci-bioinformatics/ampliconseq -c ampliconseq.config
    """.stripIndent()
    log.info ""
}

