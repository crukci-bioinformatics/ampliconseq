#!/usr/bin/env nextflow

// enable DSL 2 syntax
nextflow.enable.dsl = 2


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
            --SPECIES ${params.vepSpecies} \
            --ASSEMBLY ${params.vepAssembly} \
            --CONVERT \
            --NO_BIOPERL \
            --NO_HTSLIB \
            --NO_TEST \
            --NO_UPDATE
        """
}


workflow {
    cache_dir = channel.fromPath(params.vepCacheDir, checkIfExists: true)
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

        VEP cache directory : ${params.vepCacheDir}
        Species             : ${params.vepSpecies}
        Assembly            : ${params.vepAssembly}
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
            --help           Show this message and exit
            --vepCacheDir    Directory in which to install Ensembl VEP cache files
            --vepSpecies     The species name, e.g. homo_sapiens
            --vepAssembly    The genome assembly, e.g. GRCh37

        Alternatively, override settings using a configuration file such as the
        following:

        params {
            vepCacheDir = "vep_cache"
            vepSpecies  = "homo_sapiens"
            vepAssembly = "GRCh37"
        }

        and run as follows:
            nextflow run crukci-bioinformatics/ampliconseq -c ampliconseq.config

    """.stripIndent()
    log.info ""
}


