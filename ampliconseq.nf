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

bamDir = params.bamDir
if (bamDir && !bamDir.endsWith("/")) {
    bamDir = "${params.bamDir}/"
}


// -----------------------------------------------------------------------------
// processes
// -----------------------------------------------------------------------------

// check sample sheet is valid a create a validated version to be used in
// subsequent processes
process check_samples {
    input:
        path install_dir
        path sample_sheet

    output:
        path "samples.checked.csv"

    script:
        """
        Rscript ${install_dir}/R/check_samples_file.R ${sample_sheet} samples.checked.csv
        """
}


// create non-overlapping amplicon groups where none of the amplicons overlap
// with another amplicon from the same group
process create_non_overlapping_amplicon_groups {
    input:
        path install_dir
        path amplicon_details
        path reference_genome_index

    output:
        path "amplicon_groups.txt", emit: amplicons
        path "amplicon_groups.*.txt", emit: amplicon_groups

    script:
        """
        Rscript ${install_dir}/R/create_non_overlapping_amplicon_groups.R ${amplicon_details} ${reference_genome_index}
        """
}


// extract reads that correspond to a set of amplicon intervals to create a
// subset BAM file
process extract_amplicon_regions {
    tag "${id} ${group}"

    input:
        tuple val(id), path(bam), val(group), path(amplicons)

    output:
        tuple val(id), val(group), path(output_bam), path(output_bai), emit: bam
        tuple val(id), path(coverage), emit: coverage

    script:
        output_bam = "${id}.${group}.bam"
        output_bai = "${id}.${group}.bai"
        coverage = "${id}.amplicon_coverage.${group}.txt"
        """
        awk 'BEGIN { FS = "\t"; OFS = "\t" } NR > 1 { print \$2, \$3 - 1, \$4 }' ${amplicons} > amplicons.bed

        extract-amplicon-regions \
            --bam ${bam} \
            --intervals amplicons.bed \
            --amplicon-bam ${output_bam} \
            --coverage ${coverage} \
            --maximum-distance 0 \
            --require-both-ends-anchored \
            --unmark-duplicate-reads
        """
}


// -----------------------------------------------------------------------------
// workflow
// -----------------------------------------------------------------------------

workflow {

    install_dir = channel.fromPath(projectDir, checkIfExists: true)
    sample_sheet = channel.fromPath(params.sampleSheet, checkIfExists: true)
    amplicon_details = channel.fromPath(params.ampliconDetails, checkIfExists: true)
    reference_genome = channel.fromPath(params.referenceGenomeFasta, checkIfExists: true)
    reference_genome_index = channel.fromPath("${params.referenceGenomeFasta}.fai", checkIfExists: true)

    check_samples(install_dir, sample_sheet)

    // create groups of non-overlapping amplicons
    create_non_overlapping_amplicon_groups(
        install_dir,
        amplicon_details,
        reference_genome_index
    )

    // amplicon group channel in which each value is a tuple that includes
    // the amplicon file for that group and the group number extracted from the
    // amplicon file name
    amplicon_groups = create_non_overlapping_amplicon_groups.out.amplicon_groups
        .flatten()
        .map { file ->
            def group = (file.name =~ /(\d+)/)[0][1] as Integer
            return tuple(group, file)
        }

    // BAM file channel created by reading the validated sample sheet
    bam = check_samples.out
        .splitCsv(header: true, quote: '"')
        .map { row -> tuple("${row.ID}", file("${bamDir}${row.ID}.bam", checkIfExists: true)) }

    // combine BAM file and amplicon group channels
    // results in every combination of BAM file and amplicon group
    bam_amplicon_groups = bam
        .combine(amplicon_groups)

    // extract reads matching amplicons into a subset BAM file
    extract_amplicon_regions(bam_amplicon_groups)

    // concatenate coverage files for each BAM file
    amplicon_coverage = extract_amplicon_regions.out.coverage
        .collectFile(keepHeader: true) { id, coverage -> [ "${id}.amplicon_coverage.txt", coverage ] }

}


// -----------------------------------------------------------------------------
// summary of configuration parameters
// -----------------------------------------------------------------------------

def printParameterSummary() {
    log.info ""
    log.info """
        Variant calling pipeline for amplicon sequencing data
        =====================================================

        Sample sheet             : ${params.sampleSheet}
        BAM directory            : ${params.bamDir}
        Reference sequence FASTA : ${params.referenceGenomeFasta}
        Output directory         : ${params.outputDir}
        Output prefix            : ${params.outputPrefix}
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
            --output-dir                  Directory to which output files are written
            --output-prefix               Prefix for output file names

        Alternatively, override settings using a configuration file such as the
        following, in which parameter names used are the camelCase equivalent of the
        above options:

        params {
            sampleSheet          = "samples.csv"
            bamDir               = "bam"
            ampliconDetails      = "amplicons.csv"
            referenceGenomeFasta = "/reference_data/GRCh37.fa"
            outputDir            = "results"
            outputPrefix         = ""
        }

        and run as follows:
            nextflow run crukci-bioinformatics/ampliconseq -c ampliconseq.config
    """.stripIndent()
    log.info ""
}

