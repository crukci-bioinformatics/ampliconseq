# ampliconseq

Variant calling pipeline for amplicon sequencing data.

## Table of Contents

* [Introduction](#introduction)
* [Quickstart guide](#quickstart)
* [Installing ampliconseq](#installation)
    * [Installing a specific release or version](#install_specific_release)
    * [Updating ampliconseq](#updating)
    * [Requirements](#requirements)
    * [Running the pipeline without a container](#running_without_container)
* [Configuration](#configuration)
    * [Sample sheet](#sample_sheet)
    * [Amplicon coordinates file](#amplicon_coordinates_file)
    * [Configuration file](#configuration_file)
    * [Configuration using command line arguments](#configuration_using_command_line_args)

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

4. Create an amplicon coordinates file (`amplicons.csv`).

5. Create a configuration file (`ampliconseq.config`) specifying the sample sheet, the amplicon coordinates file, the reference genome, the VEP cache directory, the variant caller and various other parameters.

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

### <a name="running_without_container">Running the pipeline without a container</a>

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

Additionally, the pipeline contains a number of custom Java tools written
using the [HTSJDK](https://github.com/samtools/htsjdk) library. These are
available from the
[releases](https://github.com/crukci-bioinformatics/ampliconseq/releases)
page on the GitHub repository; download and unpack the tarball file named
ampliconseq-1.0.0-tools.tar.gz, substituting the version number as appropriate,
and ensure that the `bin` subdirectory is available on the `PATH`.

---

## <a name="configuration">Configuring the pipeline</a>

The ampliconseq pipeline requires a sample sheet file, an amplicon coordinates
file and, optionally, a configuration file in which the input files and
parameter settings are specified.

The input files are aligned sequence BAM files, in which there is a single BAM
file for each library and where each library contains amplified DNA for all
amplicons within the panel. The reference genome sequence FASTA file to which
the sequence reads were aligned must be specified in the configuration file or
as a command line argument; this needs to be indexed and have an accompanying
sequence dictionary.

### <a name="sample_sheet">Sample sheet</a>

The samples sheet provides details about each of the amplicon libraries. It can
be either a tab-delimited (TSV) or comma-separated value (CSV) file. The
following columns are expected.

Column  | Required | Description
--------|----------|------------
ID      | yes      | The library identifier or barcode
Sample  | yes      | The name or identifier of the sample from which the library was created
BAM     | no       | The BAM file name or path (can be a relative or absolute path)

Replicate libraries created from the same sample will share the same `Sample`
name. The pipeline creates a variant call summary table in which variants called
within replicates of the same sample are gathered and reported together with an
assigned confidence level. The QC report clusters libraries based on variant
allele fractions and identifies possible sample library mispairings where
replicate libraries are quite dissimilar to each other but more similar to other
libraries in the run.

An example sample sheet containing duplicate libraries for each of two samples
is given below (note that a single run can contain hundreds of libraries).

```
ID                  Sample        BAM
SLX-12850.FLD0011   JBLAB-2493    FLD0011.bam
SLX-12850.FLD0012   JBLAB-2493    FLD0012.bam
SLX-12850.FLD0013   JBLAB-3401    FLD0013.bam
SLX-12850.FLD0014   JBLAB-3401    FLD0014.bam
```

If the sample sheet does not contain a BAM column the pipeline will assume that
BAM files follow a file naming convention in which the ID is the prefix to which
the '.bam' extension is added, e.g. `SLX-12850.FLD0011.bam` for the first
library in the example sample sheet given above. There is a `bamDir`
configuration parameter that can be set in order to avoid having to specify the
full path for each BAM file in the sample sheet; it is prepended to the BAM file
name given in the sample sheet or to the default file name based on the ID if
the sample sheet does not contain the BAM column.

### <a name="amplicon_coordinates_file">Amplicon coordinates file</a>

### <a name="configuration_file">Configuration file</a>

### <a name="configuration_using_command_line_args">Configuration using command line arguments</a>


