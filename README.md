# ampliconseq

Variant calling pipeline for amplicon sequencing data.

## Table of Contents

* [Introduction](#introduction)
* [Quickstart guide](#quickstart)
* [Installing ampliconseq](#installation)
    * [Installing a specific release or version](#install_specific_release)
    * [Updating ampliconseq](#updating)
    * [Requirements](#requirements)
    * [Installing dependencies to run the pipeline without a container](#running_without_container)
* [Configuration](#configuration)
    * [Sample sheet](#sample_sheet)
    * [Amplicon coordinates file](#amplicon_coordinates_file)
    * [Configuration parameters](#configuration_parameters)
        * [Command line arguments](#configuration_using_command_line_args)
        * [Configuration file](#configuration_file)
    * [Reference genome sequence FASTA file](#reference_genome_fasta)
* [Running ampliconseq](#running)
    * [Running with a container](#running_with_container)
    * [Execution profiles](#execution_profiles)

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
[Cancer Research UK Cambridge Institute](https://www.cruk.cam.ac.uk) (CRUK CI).

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
            -config ampliconseq.config \
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
installed using `nextflow pull` as follows using the `-r` or `-revision` option:

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

### <a name="running_without_container">Installing dependencies to run the pipeline without a container</a>

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
Conda to install the dependencies and the Conda environment file, `conda.yml`,
located in the
[GitHub repository](https://github.com/crukci-bioinformatics/ampliconseq)
within the `docker` subdirectory can be used to install these dependencies such
that the pipeline can be run without using the container.

        conda env create -f conda.yml

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
be either a tab-delimited (TSV) or comma-separated value (CSV) file. By default,
the ampliconseq pipeline expects a file named `samples.csv` in the directory in
which the pipeline is run, but this can be can be changed within the
configuration file or using a command line argument (see
[below](#configuration_file)).

The sample sheet is expected to have the following columns.

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

The amplicon coordinates file provides the start and end coordinates for each
amplicon and the start and end position for the target region that excludes the
primer sequences. By default, the ampliconseq pipeline expects a file named
`amplicons.csv` in the directory in which the pipeline is run, but this can be
changed within the configuration file or using a command line argument (see
[below](#configuration_file)).

The following columns are all required.

Column        | Description
--------------|------------
ID            | The amplicon identifier
Chromosome    | The chromosome
AmpliconStart | The start coordinate of the amplicon (includes primer sequence)
AmpliconEnd   | The end coordinate of the amplicon (includes primer sequence)
TargetStart   | The start coordinate of the target region (excludes primer)
TargetEnd     | The end coordinate of the target region (excludes primer)

The following snippet from an amplicon coordinates file contains a set of
amplicons targeting the TP53 gene.

```
ID              Chromosome  AmpliconStart  AmpliconEnd  TargetStart  TargetEnd
TP53_D0008_001  chr17       7572903        7573031      7572924      7573010
TP53_D0008_002  chr17       7573904        7574019      7573922      7574000
TP53_D0008_003  chr17       7573975        7574077      7573997      7574049
TP53_D0008_004  chr17       7576789        7576917      7576812      7576898
TP53_D0008_005  chr17       7576873        7576961      7576895      7576936
TP53_D0008_006  chr17       7576996        7577115      7577015      7577097
TP53_D0008_007  chr17       7577074        7577182      7577094      7577157
TP53_D0008_008  chr17       7577434        7577528      7577453      7577506
TP53_D0008_009  chr17       7577484        7577612      7577503      7577590
TP53_D0008_010  chr17       7577561        7577667      7577587      7577649
```

### <a name="configuration_parameters">Configuration parameters</a>

The ampliconseq pipeline has a number of configuration parameters. Use the
`--help` option to see usage instructions and details of each.

    nextflow run crukci-bioinformatics/ampliconseq --help

These can be set either as command-line options or using a configuration file.

#### <a name="configuration_using_command_line_args">Configuration using command line arguments</a>

It is possible to set configuration parameters using command line arguments.
This can become unwieldy when changing a large number of settings and the
alternative use of a configuration file is generally preferred.

As an example, the path to the reference genome sequence FASTA file, to which
the sequence data were aligned, can be specified using the
`referenceGenomeFasta` paramter as follows:

    nextflow run crukci-bioinformatics/ampliconseq --referenceGenomeFasta /data/reference_data reference_genomes/homo_sapiens/GRCh37/fasta/GRCh37.fa 

The following example specifies the sample sheet and amplicon coordinates file,
and instructs the pipeline to annotate variants using Ensembl VEP for the given
species and assembly. It also specifies the variant caller to use (VarDict) and
the minimum allele fraction of variants to be called.

    nextflow run crukci-bioinformatics/ampliconseq \
        --sampleSheet samples.txt \
        --ampliconDetails /data/reference_data/ampliconseq/tp53_panel/amplicons.csv \
        --referenceGenomeFasta /data/reference_data/reference_genomes/homo_sapiens/GRCh37/fasta/GRCh37.fa \
        --vepAnnotation \
        --vepCacheDir vep_cache \
        --vepSpecies homo_sapiens \
        --vepAssembly GRCh37 \
        --variantCaller vardict \
        --minimumAlleleFraction 0.01

The default parameter settings can be found in the
[`nextflow.config`](nextflow.config) file in the
[GitHub repository](https://github.com/crukci-bioinformatics/ampliconseq).

#### <a name="configuration_file">Configuration file</a>

A more convenient way of setting pipeline parameters makes use of a
configuration file, e.g. `ampliconseq.config`. This is specified when running
the pipeline using the `-config` (or `-c`) option.

    nextflow run crukci-bioinformatics/ampliconseq -c ampliconseq.config

The following is a sample configuration file that sets the same sample sheet,
amplicon coordinates file, VEP annotation and variant calling settings as in
the above example using command line arguments.

```
params {
    sampleSheet           = "samples.txt"
    ampliconDetails       = "/data/reference_data/ampliconseq/tp53_panel/amplicons.csv"
    referenceGenomeFasta  = "/data/reference_data/reference_genomes/homo_sapiens/GRCh37/fasta/GRCh37.fa"
    vepAnnotation         = true
    vepSpecies            = "homo_sapiens"
    vepAssembly           = "GRCh37"
    outputDir             = "results"
    variantCaller         = "vardict"
    minimumAlleleFraction = 0.01
}
```

This takes the form of a `name = value` syntax with a separate line for each
parameter within curly braces bounding a `params` block. Note that file names
and paths and other character or string values need to be quoted while numeric
values do not, and boolean parameters such as `vepAnnotation` can be set to
`true` or `false`.

See the [Nextflow documentation](https://www.nextflow.io/docs/latest) for more
details about the configuration syntax.

The configuration file will normally contain a subset of the parameters
specified in the [`nextflow.config`](nextflow.config) found at the top level of
the [GitHub repository](https://github.com/crukci-bioinformatics/ampliconseq).
[`nextflow.config`](nextflow.config) contains the default settings that are
overridden by the configuration file specified with the `-config` option when
running the pipeline.

### <a name="reference_genome_fasta">Reference genome sequence FASTA file</a>

The reference genome FASTA file using in aligning the sequence data to generate
the BAM files, the primary input to the pipeline, needs to be specified using
with the `referenceGenomeFasta` parameter. This FASTA file needs to be indexed,
e.g. using `samtools faidx`, and have an accompanying sequence dictionary that
can be created with `samtools dict` or the GATK/Picard CreateSequenceDictionary
tool.

### <a name="ensembl_vep_cache">Ensembl Variant Effect Predictor (VEP) cache</a>

The pipeline can annotate variants using Ensembl VEP. It runs VEP in offline
mode using a pre-downloaded annotation cache. The cache for a particular species
and genome assembly can be downloaded using a single step supplementary workflow
(`download_vep_cache`) as shown below:

    nextflow run crukci-bioinformatics/ampliconseq \
        -main-script download_vep_cache.nf \
        -with-singularity \
        --vepCacheDir /path_to/vep_cache \
        --vepSpecies homo_sapiens \
        --vepAssembly GRCh37

This assumes use of VEP installed in the container and run using Singularity.
Use `-with-docker` to use Docker instead of Singularity or remove the
`-with-singularity` argument if not using a container (the `vep_install`
installed with VEP will need to be available on the `PATH`).

Substitute the top-level VEP cache directory as required; note that this
directory must already exist.

The Ensembl VEP cache is quite large (around 15G for homo sapiens) and the
download and unpacking carried out in this step can take several minutes.

## <a name="running">Running ampliconseq</a>

### <a name="running_with_container">Running the pipeline using a container</a>

The most straightforward way to run the ampliconseq pipeline is to use the
pre-packaged container either using Docker or Singularity by specifying the
`-with-docker` or `-with-singularity` flag.

    nextflow run crukci-bioinformatics/ampliconseq -with-docker

    nextflow run crukci-bioinformatics/ampliconseq -with-singularity

Nextflow will automatically fetch the container from
[Docker Hub](https://hub.docker.com/r/crukcibioinformatics/ampliconseq) and,
if using Singularity, will build the Singularity image from the Docker
container. Singularity is more likely than Docker to be available on
high-performance cluster computing platforms.

When using Singularity, the pipeline assumes that the user bind control feature
is enabled and sets `singularity.autoMounts = true` in the Nextflow
configuration file. See the
[Nextflow documentation](https://www.nextflow.io/docs/latest/index.html) for
more details on this.

Alternatively, the use of the Docker container can be specified in the
configuration file by adding the following line:

    docker.enabled = true

Similarly, to enable Singularity, instead add the following line:

    singularity.enabled = true

These can also be added as part of an execution profile (see next section).

### <a name="execution_profiles">Execution profiles</a>

Resource settings are configured using Nextflow profiles. The ampliconseq
pipeline provides three profiles - `standard`, `bigserver` and `cluster`
configured for running on servers and the high-performance compute cluster at
CRUK CI. These specifythe maximum number of CPUs or memory that can be used at
any one time during the pipeline run or the maximum number of jobs that can be
submitted to the cluster to be run in parallel.

A custom profile can be created in the configuration file, e.g.
`ampliconseq.config`, an example of which is shown below.

    myprofile {
        process {
            executor = 'slurm'
            queue = 'long'
        }
        executor {
            queueSize = 25
            pollInterval = 30.sec
            jobName = { "'$task.name'" }
        }
        singularity.enabled = true
    }

This profile can be specified using the `-profile` command line option.

    nextflow run crukci-bioinformatics/ampliconseq -c ampliconseq.config -profile myprofile

With this profile, Nextflow will submit jobs to cluster nodes using the SLURM
resource manager. A maximum number of 25 jobs that will be submitted to the
'long' queue for running in parallel and Nextflow will poll every 30 seconds
to check for completed jobs. Use of Singularity for running jobs using the
container is enabled so it is not necessary to specify this separately with the
`-with-singularity` option.
