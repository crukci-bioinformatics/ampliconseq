# ampliconseq

Variant calling pipeline for amplicon sequencing data.

## Table of Contents

* [Introduction](#introduction)
* [Quickstart guide](#quickstart)
* [Installing ampliconseq](#installation)
    * [Installing ampliconseq](#installation)
    * [Installing a specific version/release](#install_specific_release)
    * [Requirements](#requirements)

---

## <a name="introduction">Introduction</a>

*ampliconseq* is an analysis pipeline for calling single nucleotide variants
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

The ampliconseq pipeline is run using the [Nextflow](https://www.nextflow.io)
scientic workflow system and all dependencies and tools are packaged in a
[Docker](https://www.docker.com) container that can be run either using Docker
or [Singularity](https://sylabs.io/docs). The inputs to the pipeline are BAM
files containing aligned sequence reads.

The ampliconseq pipeline has the following features:

* Choice of variant callers: [GATK HaplotypeCaller](https://gatk.broadinstitute.org) and [VarDict](https://github.com/AstraZeneca-NGS/VarDictJava)
* Alignment and coverage QC report using metrics calculated by [Picard](https://broadinstitute.github.io/picard)
* Annotation of variants using [Ensembl Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html)
* Support for overlapping amplicon targets by partitioning reads prior to variant calling
* Support for calling and filtering low allele fraction SNVs, e.g. for circulating tumour DNA in plasma samples with allele fractions down to 0.1%, by fitting probability distributions to model background noise
* *Specific calling* of known mutations
* Combines variant call outputs for replicate amplicon libraries from the same sample and assigns confidence levels based on whether calls are made or filtered in each of the replicates
* Minimal barrier to installation with the only requirements being a Java runtime, [Nextflow](https://www.nextflow.io), and either [Docker](https://www.docker.com) or [Singularity](https://sylabs.io/docs) to run a container in which all other dependencies and tools are packaged
* Simple installation, only Java, Nextflow and either Docker or Singularity required
* Scale easily from deployment on multi-core workstation to high-performance compute cluster or cloud with simple configuration change
* Accompanying visualization tool for viewing and assessing SNV calls

The ampliconseq pipeline was developed by the
[Bioinformatics Core](https://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core)
in collaboration with James Brenton's research group at the
[Cancer Research UK Cambridge Institute](https://www.cruk.cam.ac.uk).

---

## <a name="quickstart">Quickstart guide</a>

1. Install Nextflow (Java 8 or later required).

        curl -s https://get.nextflow.io | bash

    This creates a file named `nextflow` in the current directory. For
    convenience, move this to some location on your `PATH`. For more details on
    installing Nextflow, see the
    [Nextflow documentation](https://www.nextflow.io).

2. Download Ensembl VEP cache.

    This step can be skipped if the VEP cache for the relevant species and genome
assembly is already installed or if variant annotation is not required.

        nextflow run crukci-bioinformatics/ampliconseq \
            -main-script download_vep_cache.nf \
            -with-singularity \
            --vepCacheDir /path_to/vep_cache \
            --vepSpecies homo_sapiens \
            --vepAssembly GRCh37

    This will download the ampliconseq pipeline from
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

3. Create a samples file (`samples.txt`) containing library and sample identifiers.

4. Create an amplicon definition file (`amplicons.csv`) providing amplicon and target region coordinates for each amplicon.

5. Create a configuration file (`ampliconseq.config`) specifying the sample sheet, the amplicon definition file, the reference genome, the VEP cache directory, the variant caller and various other parameters.

6. Run the ampliconseq pipeline specifying the configuration file and execution profile.

        nextflow \
            run crukci-bioinformatics/ampliconseq \
            -c ampliconseq.config \
            -with-singularity \
            -profile bigserver \
            -with-report logs/ampliconseq_report.html \
            -with-timeline logs/ampliconseq_timeline.html

---

## <a name="installation">Installing ampliconseq</a>

The ampliconseq pipeline is downloaded and run using the Nextflow workflow
system. Dependencies, including GATK, VarDict, Picard, Ensembl Variant Effect
Predictor, R and various R packages are packaged as a
[Docker container](https://hub.docker.com/r/crukcibioinformatics/ampliconseq)
that can be run with either [Docker](https://www.docker.com) or
[Singularity](https://sylabs.io/docs). The container is also downloaded by
Nextflow. The only requirements are a recent version of Nextflow and either
Docker or Singularity. Nextflow requires Java 8 or above and can be installed as
shown in the Quickstart section above. See the
[Nextflow documentation](https://www.nextflow.io/docs/latest/index.html) for
more details.

### <a name="install_specific_release">Installing a specific release of ampliconseq</a>

Using the latest stable
[release](https://github.com/crukci-bioinformatics/ampliconseq/releases)
of ampliconseq is recommended. A specific version of ampliconseq can be
installed using `nextflow pull` as follows:

    nextflow pull crukci-bioinformatics/ampliconseq -r 1.0.0

When a specific version of ampliconseq is installed in this way the revision
also needs to be specified when running the pipeline using `nextflow run`.

    nextflow run crukci-bioinformatics/ampliconseq -r 1.0.0 -c ampliconseq.config

Run `nextflow info` to view details about the currently installed version.

    nextflow info crukci-bioinformatics/ampliconseq

### <a name="updating">Updating ampliconseq</a>

The latest snapshot of ampliconseq will be downloaded and run if no revision is
specified using the `-r` or `-revision` command line option when running
ampliconseq for the first time. Subsequent runs will use this snapshot version
but Nextflow detects if there have been revisions to the pipeline since then and
displays a message such as the following:

    NOTE: Your local project version looks outdated - a different revision is available in the remote repository [961d1d72a2]

Run the following command to update ampliconseq to the latest revision on the
master branch:

    nextflow pull crukci-bioinformatics/ampliconseq -r master

### <a name="requirements">Requirements</a>

* [Nextflow](https://www.nextflow.io) 20.10.0 or above
* [Singularity](https://sylabs.io/docs) or [Docker](https://www.docker.com)

Dependencies, including GATK, VarDict, Ensembl VEP, R and various R packages,
are packaged in a
[Docker container](https://hub.docker.com/r/crukcibioinformatics/ampliconseq)
that will be downloaded automatically by Nextflow.

ampliconseq can be run without a container by installing the following tools
and packages:

* R 4.1.0 or above and the following packages:
    * tidyverse
    * optparse
    * fitdistrplus
    * nozzle.r1
    * base64
    * svglite
    * rsvg
    * ComplexHeatmap (from Bioconductor)
* GATK 4.2.0.0 or above
* VarDict (Java version)
* Ensembl Variant Effect Predictor

These can be installed manually or, more straightforwardly, using Conda. The
pipeline assumes that the executables, `R`, `gatk`, `vardict-java` and `vep`,
are available on your `PATH`. The Docker container recipe (Dockerfile) uses
Conda to install the dependencies and the Conda environment file located in
the [GitHub repository](https://github.com/crukci-bioinformatics/ampliconseq)
within the `docker` subdirectory can be used to install these dependencies such
that the pipeline can be run without using the container.



