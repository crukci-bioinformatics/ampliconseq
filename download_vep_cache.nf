#!/usr/bin/env nextflow

// enable DSL 2 syntax
nextflow.enable.dsl = 2

// -----------------------------------------------------------------------------
// default parameter settings
// -----------------------------------------------------------------------------

params.help     = false
params.cacheDir = "${launchDir}/vep_cache"
params.species  = "homo_sapiens"
params.assembly = "GRCh37"


printParameterSummary()

if (params.help) {
    helpMessage()
    exit 0
}


process download_vep_cache {

    input:
        path cache_dir

    script:
        """
        vep_install \
            --AUTO cf \
            --CACHEDIR ${cache_dir} \
            --SPECIES ${params.species} \
            --ASSEMBLY ${params.assembly} \
            --CONVERT \
            --NO_BIOPERL \
            --NO_HTSLIB \
            --NO_TEST \
            --NO_UPDATE
        """
}


workflow {
    cache_dir = channel.fromPath(params.cacheDir, checkIfExists: true)
    download_vep_cache(cache_dir)
}


// -----------------------------------------------------------------------------
// summary of configuration parameters
// -----------------------------------------------------------------------------

def printParameterSummary() {
    log.info ""
    log.info """
        Variant calling pipeline for amplicon sequencing data
        =====================================================

        Ensembl VEP cache download

        VEP cache directory : ${params.cacheDir}
        Species             : ${params.species}
        Assembly            : ${params.assembly}
    """.stripIndent()
    log.info ""
}


// ----------------------------------------------------------------------------
// help/usage
// ----------

def helpMessage() {
    log.info """
        Usage:
            nextflow run crukci-bioinformatics/ampliconseq --main-script download_vep_cache.nf

        Options:
            --help            Show this message and exit
            --cache-dir       Directory in which to install Ensembl VEP cache files
            --species         The species name, e.g. homo_sapiens
            --assembly        The genome assembly, e.g. GRCh37

    """.stripIndent()
    log.info ""
}


