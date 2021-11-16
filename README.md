# ampliconseq

Variant calling pipeline for amplicon sequencing data.

## Introduction

`ampliconseq` is an analysis pipeline for calling single nucleotide variants
(SNVs) and indels in tagged amplicon sequencing data. Variants are called using
[GATK HaplotypeCaller](https://gatk.broadinstitute.org), preferred for germline
or clonal somatic mutations especially in FFPE samples, or
[VarDict](https://github.com/AstraZeneca-NGS/VarDictJava) which can identify low
allele fraction (AF) SNVs in circulating tumour DNA in plasma samples. In
addition to caller-specific filters, the pipeline models the background
substitution noise at each amplicon position to identify and filter erroneous
SNV calls. Alignment and target coverage metrics are computed and compiled into
a QC report. Variants are annotated using
[Ensembl Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html)
(VEP).

The `ampliconseq` pipeline is run using the [Nextflow](https://www.nextflow.io)
scientic workflow system and all dependencies and tools are packaged in a
[Docker container](https://hub.docker.com/r/crukcibioinformatics/ampliconseq).
The inputs to the pipeline are BAM files containing aligned sequence reads.

The `ampliconseq` pipeline has the following features:

* Choice of variant callers: [GATK HaplotypeCaller](https://gatk.broadinstitute.org) and [VarDict](https://github.com/AstraZeneca-NGS/VarDictJava)
* Alignment and coverage QC report using metrics calculated by [Picard](https://broadinstitute.github.io/picard)
* Annotation of variants using [Ensembl Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html)
* Support for overlapping amplicon targets by partitioning reads prior to variant calling
* Support for calling and filtering low allele fraction SNVs, e.g. for circulating tumour DNA in plasma samples with allele fractions down to 0.1%, by fitting probability distributions to model background noise
* *Specific calling* of known mutations
* Minimal barrier to installation with the only requirements being a Java runtime, Nextflow, and either [Docker](https://www.docker.com) or [Singularity](https://sylabs.io/docs) to run a container in which all other dependencies and tools are packaged
* Simple installation, only Java, Nextflow and either Docker or Singularity required
* Scale easily from deployment on multi-core workstation to high-performance compute cluster or cloud with simple configuration change
* Accompanying visualization tool for viewing and assessing SNV calls

The `ampliconseq` pipeline was developed by the
[Bioinformatics Core](https://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core)
in collaboration with James Brenton's research group at the
[Cancer Research UK Cambridge Institute](https://www.cruk.cam.ac.uk).

---

## Quickstart

1. Install Nextflow (Java 8 or later required).

        curl -s https://get.nextflow.io | bash

2. Download Ensembl VEP cache.

    This step can be skipped if the VEP cache for the relevant species and genome
assembly is already installed.

        nextflow -bg -q run crukci-bioinformatics/ampliconseq \
            -main-script download_vep_cache.nf \
            -with-singularity \
            --vepCacheDir /path_to/vep_cache \
            --vepSpecies homo_sapiens \
            --vepAssembly GRCh37

    This will download the `ampliconseq` pipeline from
[GitHub](https://github.com/crukci-bioinformatics/ampliconseq) including the
single step workflow for downloading the VEP cache (`download_vep_cache.nf`).
It will also download the Docker container in which Ensembl VEP is installed
from [Docker Hub](https://hub.docker.com/r/crukcibioinformatics/ampliconseq)
from which a Singularity container is built. Use `-with-docker` to use Docker
instead of Singularity.

    Substitute the top-level VEP cache directory as required; note that this
directory must exist.

    The Ensembl VEP cache is quite large (around 15G for homo sapiens) and the
download and unpacking carried out in this step can take several minutes.

3. Create a configuration file (`ampliconseq.config`) specifying the sample sheet, the amplicon definition file, the reference genome, the VEP cache directory, the variant caller and various other parameters.

4. Create a sample sheet.

5. Create an amplicon definition file, providing amplicon and target region coordinates for each amplicon.

6. Run the `ampliconseq` pipeline specifying the configuration file and execution profile.

        nextflow \
            -c ampliconseq.config \
            -log logs/ampliconseq.log \
            run crukci-bioinformatics/ampliconseq \
            -with-singularity \
            -profile bigserver \
            -with-report logs/ampliconseq_report.html \
            -with-timeline logs/ampliconseq_timeline.html

---

