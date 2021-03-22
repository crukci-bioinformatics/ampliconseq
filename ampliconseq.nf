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
// processes
// -----------------------------------------------------------------------------

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


process create_non_overlapping_interval_groups {
    input:
        path install_dir
        path amplicon_details
        path reference_genome_index

    output:
        path "amplicon_groups.csv", emit: amplicons
        path "amplicons.*.bed", emit: amplicon_bed_files
        path "targets.*.bed", emit: target_bed_files

    script:
        """
        Rscript ${install_dir}/R/create_non_overlapping_interval_groups.R ${amplicon_details} ${reference_genome_index}
        """
}


// -----------------------------------------------------------------------------
// workflow
// -----------------------------------------------------------------------------

workflow {

    install_dir = channel.fromPath(projectDir, checkIfExists: true)
    sample_sheet = channel.fromPath(params.sampleSheet, checkIfExists: true)
    amplicon_details = channel.fromPath(params.ampliconDetails, checkIfExists: true)
    reference_genome_fasta = channel.fromPath(params.referenceGenomeFasta, checkIfExists: true)
    reference_genome_index = channel.fromPath("${params.referenceGenomeFasta}.fai", checkIfExists: true)

    check_samples(install_dir, sample_sheet)

    create_non_overlapping_interval_groups(
        install_dir,
        amplicon_details,
        reference_genome_index
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

