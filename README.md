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
* [Reference data](#reference_data)
    * [Reference genome sequence FASTA file](#reference_genome_fasta)
    * [Ensembl Variant Effect Predictor (VEP) cache](#ensembl_vep_cache)
* [Running ampliconseq](#running)
    * [Running with a container](#running_with_container)
    * [Execution profiles](#execution_profiles)
    * [Nextflow reports](#nextflow_reports)
    * [Nextflow log files and work directories](#nextflow_log_files_and_work_directories)

---

## <a name="introduction">Introduction</a>

*ampliconseq* is an analysis pipeline for calling single nucleotide variants
(SNVs) and indels in targeted amplicon sequencing data. Variants are called
using [GATK HaplotypeCaller](https://gatk.broadinstitute.org), preferred for
germline or clonal somatic mutations, especially in FFPE samples, or
[VarDict](https://github.com/AstraZeneca-NGS/VarDictJava) which can identify low
allele fraction SNVs in circulating tumour DNA from plasma samples. In addition
to caller-specific filters, the pipeline models the background substitution
noise at each amplicon position to identify and filter SNV calls with very low
allele fractions that are not distinguishable from noise. Alignment and target
coverage metrics are computed and compiled into a QC report. Variants are
annotated using
[Ensembl Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html)
(VEP).

The ampliconseq pipeline is executed using the
[Nextflow](https://www.nextflow.io) scientic workflow system and all
dependencies and tools are packaged in a [Docker](https://www.docker.com)
container that can be run either using Docker or
[Singularity](https://sylabs.io/docs). The inputs to the pipeline are BAM files
containing sequence reads aligned to the reference genome.

The ampliconseq pipeline has the following features:

* Choice of variant callers: [GATK HaplotypeCaller](https://gatk.broadinstitute.org) and [VarDict](https://github.com/AstraZeneca-NGS/VarDictJava)
* Alignment and coverage QC report using metrics calculated by [Picard](https://broadinstitute.github.io/picard) CollectAlignmentSummaryMetrics and CollectTargetedPcrMetrics
* Annotation of variants using [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html)
* Support for overlapping amplicon targets by partitioning reads prior to variant calling
* Support for calling and filtering low allele fraction SNVs, e.g. for circulating tumour DNA in plasma samples with allele fractions down to 0.1%, by fitting probability distributions to model background noise
* *Specific calling* of known mutations
* Assignment of confidence level based on whether a variant is called or filtered in each of a set of replicate libraries (usually duplicate libraries)
* Minimal barrier to installation with the only requirements being a Java runtime, [Nextflow](https://www.nextflow.io), and either [Docker](https://www.docker.com) or [Singularity](https://sylabs.io/docs) to run a container in which all other dependencies and tools are packaged
* Scales easily from deployment on multi-core workstation to high-performance compute cluster or cloud with only a simple configuration change
* Accompanying visualization tool for viewing and assessing SNV calls

The ampliconseq pipeline was developed by the
[Bioinformatics Core](https://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core)
in collaboration with
[James Brenton's research group](https://www.cruk.cam.ac.uk/research-groups/brenton-group)
at the
[Cancer Research UK Cambridge Institute](https://www.cruk.cam.ac.uk) (CRUK CI).

---

## <a name="quickstart">Quickstart guide</a>

1. Install Nextflow (Java 8 or later required).

        curl -s https://get.nextflow.io | bash

    This creates a file named `nextflow` in the current directory. For
    convenience, move this to some location on your `PATH`. See the
    [Nextflow documentation](https://www.nextflow.io) for more details on
    installing Nextflow.

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
and from this build a Singularity container. Use `-with-docker` to use Docker
instead of Singularity.

    Substitute the top-level VEP cache directory as required; note that this
step will fail if the directory doesn't already exist.

    The VEP cache can be quite large (around 15G for homo sapiens) and
downloading and unpacking the cache may take several minutes.

3. Create a samples file (`samples.txt`) containing the sample identifier and BAM file for each library.

4. Create an amplicon coordinates file (`amplicons.csv`).

5. Create a configuration file (`ampliconseq.config`) specifying the sample sheet, the amplicon coordinates file, the reference genome, the VEP cache directory, the variant caller and various other parameters.

6. Run the ampliconseq pipeline specifying the configuration file and execution profile.

        nextflow run crukci-bioinformatics/ampliconseq \
            -config ampliconseq.config \
            -with-singularity \
            -profile bigserver \
            -with-report ampliconseq_report.html \
            -with-timeline ampliconseq_timeline.html

---

## <a name="installation">Installing ampliconseq</a>

The ampliconseq pipeline is downloaded and run using the Nextflow workflow
system. Dependencies, including GATK, VarDict, Picard, Ensembl Variant Effect
Predictor, R and various R packages, are packaged as a
[Docker container](https://hub.docker.com/r/crukcibioinformatics/ampliconseq)
that can be run with either [Docker](https://www.docker.com) or
[Singularity](https://sylabs.io/docs). The container is also downloaded by
Nextflow. The only requirements are a recent version of Nextflow and either
Docker or Singularity. Nextflow requires Java 8 or above and can be installed as
shown in the Quickstart section above (see the
[Nextflow documentation](https://www.nextflow.io/docs/latest/index.html) for
more details).

### <a name="install_specific_release">Installing a specific release of ampliconseq</a>

Using the latest stable
[release](https://github.com/crukci-bioinformatics/ampliconseq/releases)
of ampliconseq is recommended. A specific version of ampliconseq can be
installed using `nextflow pull` with the `-revision` (or `-r`) option:

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
but Nextflow detects if there have been revisions to the pipeline since and
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
* GATK 4.2.0.0 or above (includes the Picard tools used to calculate various metrics)
* VarDict (Java version) 1.8.2 or above
* Ensembl Variant Effect Predictor release 104 or later

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
BAM     | no       | The BAM file name (where the directory can be specified as a configuration parameter) or path (can be a relative or absolute path)

Replicate libraries created from the same sample will share the same `Sample`
name or identifier. This sample-based grouping of libraries is used in the
pipeline when creating the variant call summary table, in which variants called
within replicates of the same sample are gathered and reported together with a
confidence level. Is is also used when identifying possible sample library
mispairings as part of the QC report; libraries are clustered based on variant
allele fractions and replicate libraries from the same sample are expected to
cluster together.

An example sample sheet containing duplicate libraries for each of two samples
is given below. This is a small snippet of a sample sheet; runs typically
contain tens or hundreds of libraries.

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
the sample sheet does not contain a BAM column.

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


The following parameters can be configured. These can be set either as
[command line options](#configuration_using_command_line_args) or using a
[configuration file](#configuration_file).


parameter                                | default value   | description
-----------------------------------------|-----------------|-----------------------------------------
samples                                  | samples.csv     | CSV or TSV file giving the sample name and BAM file for each library (ID and Sample columns required, optional BAM column).
bamDir                                   |                 | Directory in which BAM files are located; paths to BAM files specified in the sample sheet are relative to this directory or to the launch directory if not specified. Alternatively, this parameter can be left unset and full paths given in the BAM column within the samples file.
amplicons                                | amplicons.csv   | CSV/TSV file containing amplicon coordinates (ID, Chromosome, AmpliconStart, AmpliconEnd, TargetStart, TargetEnd columns required).
specificVariants                         |                 | CSV/TSV file containing specific (or known) variants that are included in the summary regardless of whether these are called or not (Sample, Chromosome, Position, Ref, Alt columns required).
blacklistedVariants                      |                 | CSV/TSV file containing blacklisted variants that will be filtered (Chromosome, Position, Ref, Alt columns required).
referenceGenomeFasta                     | /reference_data/GRCh37.fa | FASTA file containing the reference genome sequence (must be indexed and have an accompanying sequence dictionary).
vepAnnotation                            | false           | Annotate variants with Ensembl VEP.
vepCacheDir                              | /reference_data/vep_cache | Directory in which Ensembl VEP cache files are installed.
vepSpecies                               | homo_sapiens    | The species name of the VEP annotation cache.
vepAssembly                              | GRCh37          | The genome assembly of the VEP annotation cache.
outputDir                                |                 | Directory to which output files are written or the launch directory if not specified.
variantCaller                            | VarDict         | The variant caller (VarDict or HaplotypeCaller).
minimumAlleleFraction                    | 0.01            | Lower allele fraction limit for detection of variants (for variant callers that provide this option only).
maximumReadsPerAlignmentStart            | 2500            | Maximum number of reads to retain per alignment start position; reads above this threshold will be downsampled (specific to GATK HaplotypeCaller).
minimumMappingQualityForPileup           | 1               | Minimum mapping quality of reads to include in the pileup, i.e. when computing depths and allele fractions.
minimumBaseQualityForPileup              | 10              | Minimum base quality at a given locus for reads to include in the pileup, i.e. when computing depths and allele fractions.
minimumDepthForBackgroundNoise           | 100             | Minimum depth of coverage at a given locus for a library to be included when computing background noise.
excludeHighestFractionForBackgroundNoise | 0.1             | Fraction of measurements with the highest allele fraction to exclude from fitting a distribution to the background noise (assumes these are not due to error/noise).
maximumAlleleFractionForBackgroundNoise  | 0.03            | Maximum allele fraction to include in fitting a distribution to the background noise (assumes anything above this is not due to error/noise).
minimumNumberForFittingBackgroundNoise   | 10              | Minimum number of libraries required to fit a background noise distribution.
chunkSizeForFittingBackgroundNoise       | 500000          | Maximum number of pileup count rows to process in a chunk when fitting background noise distributions.
readChunkSizeForFittingBackgroundNoise   | 100000          | Chunk size for reading pileup count records prior to chunking for fitting background noise distributions.
sequenceContextLength                    | 5               | The length of the sequence context bordering the variant on the 5' and 3' ends to be included in the output table.
minimumDepthForHighConfidenceCalls       | 100             | Minimum depth for high-confidence variant calls.
jvmOverhead                              | 192             | The memory overhead to allow for the Java Virtual Machine in addition to the memory specified for each Java process.

#### <a name="configuration_using_command_line_args">Configuration using command line arguments</a>

It is possible to set configuration parameters using command line arguments.
This can become unwieldy when changing a large number of settings and the
alternative use of a configuration file is generally preferred.

As an example, the path to the reference genome sequence FASTA file, to which
the sequence data were aligned, can be specified using the
`referenceGenomeFasta` paramter as follows:

    nextflow run crukci-bioinformatics/ampliconseq --referenceGenomeFasta /data/reference_data/reference_genomes/homo_sapiens/GRCh37/fasta/GRCh37.fa

The following example specifies the sample sheet and amplicon coordinates file,
and instructs the pipeline to annotate variants using Ensembl VEP for the given
species and assembly. It also specifies the variant caller to use (VarDict) and
the minimum allele fraction of variants that it can attempt to identify.

    nextflow run crukci-bioinformatics/ampliconseq \
        --samples samples.txt \
        --amplicons /data/reference_data/ampliconseq/tp53_panel/amplicons.csv \
        --referenceGenomeFasta /data/reference_data/reference_genomes/homo_sapiens/GRCh37/fasta/GRCh37.fa \
        --vepAnnotation \
        --vepCacheDir /data/reference_data/vep_cache \
        --vepSpecies homo_sapiens \
        --vepAssembly GRCh37 \
        --variantCaller vardict \
        --minimumAlleleFraction 0.01

The default parameter settings can be found in the
[`nextflow.config`](nextflow.config) file that is installed as part of the
pipeline.

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
    samples               = "samples.txt"
    amplicons             = "/data/reference_data/ampliconseq/tp53_panel/amplicons.csv"
    referenceGenomeFasta  = "/data/reference_data/reference_genomes/homo_sapiens/GRCh37/fasta/GRCh37.fa"
    vepAnnotation         = true
    vepCacheDir           = "/data/reference_data/vep_cache"
    vepSpecies            = "homo_sapiens"
    vepAssembly           = "GRCh37"
    outputDir             = "results"
    variantCaller         = "vardict"
    minimumAlleleFraction = 0.01
}
```

This takes the form of a `name = value` syntax with a separate line for each
parameter within curly braces bounding a `params` block. Note that file names
and paths and other character or string values need to be in quotation marks
while numeric values do not, and boolean parameters such as `vepAnnotation` can
be set to `true` or `false`.

See the [Nextflow documentation](https://www.nextflow.io/docs/latest) for more
details about the configuration syntax.

The configuration file will normally contain a subset of the parameters
specified in the [`nextflow.config`](nextflow.config) found at the top level of
the [GitHub repository](https://github.com/crukci-bioinformatics/ampliconseq).
[`nextflow.config`](nextflow.config) contains the default settings, some or all
of which are overridden by the configuration file specified with the `-config`
option when running the pipeline.

## <a name="reference_data">Reference data</a>

The pipeline requires a reference genome FASTA file and, optionally, an
annotation database or cache for Ensembl Variant Effect Predictor (VEP).

### <a name="reference_genome_fasta">Reference genome sequence FASTA file</a>

The primary inputs to the pipeline are BAM files containing alignments for
sequence reads mapped to a reference genome. The reference genome sequence FASTA
file used in the alignment process must be specified using the
`referenceGenomeFasta` parameter. This FASTA file needs to be indexed,
e.g. using `samtools faidx`, and have an accompanying sequence dictionary that
can be created with `samtools dict` or the GATK/Picard CreateSequenceDictionary
tool.

### <a name="ensembl_vep_cache">Ensembl Variant Effect Predictor (VEP) cache</a>

The pipeline can annotate variants using Ensembl VEP. It runs VEP in an offline
mode using a pre-downloaded annotation cache. The cache for a particular species
and genome assembly can be downloaded using a single step supplementary workflow
(`download_vep_cache`) as shown below:

    nextflow run crukci-bioinformatics/ampliconseq \
        -main-script download_vep_cache.nf \
        -with-singularity \
        --vepCacheDir /path_to/vep_cache \
        --vepSpecies homo_sapiens \
        --vepAssembly GRCh37

The `-with-singularity` argument indicates that VEP will be run in a container
using Singularity. Use `-with-docker` to use Docker instead of Singularity or
remove the `-with-singularity` argument if not using a container, in which case
the `vep_install` tool that is installed with VEP will need to be available on
the `PATH`.

Substitute the top-level VEP cache directory as required; note that the download
will fail if the directory doesn't already exist.

The VEP cache can be quite large (around 15G for homo sapiens) and downloading
and unpacking the cache may take several minutes.

## <a name="running">Running ampliconseq</a>

### <a name="running_with_container">Running the pipeline using a container</a>

The most straightforward way to run the ampliconseq pipeline is to use the
pre-packaged container with either Docker or Singularity by specifying the
`-with-docker` or `-with-singularity` flag.

    nextflow run crukci-bioinformatics/ampliconseq -config ampliconseq.config -with-docker

    nextflow run crukci-bioinformatics/ampliconseq -config ampliconseq.config -with-singularity

Nextflow will automatically fetch the container from
[Docker Hub](https://hub.docker.com/r/crukcibioinformatics/ampliconseq) and will
build the Singularity image from the Docker container when running with
Singularity. Singularity is more likely than Docker to be available on
high-performance cluster computing platforms.

When using Singularity, the pipeline assumes that the user bind control feature
is enabled and sets `singularity.autoMounts = true` in the Nextflow
configuration file. See the
[Nextflow documentation](https://www.nextflow.io/docs/latest/index.html) for
more details on this.

Alternatively, the use of the Docker container can be specified in the
configuration file by adding the following line:

    docker.enabled = true

Similarly, to enable Singularity, instead add the following:

    singularity.enabled = true

These can also be added as part of an execution profile (see next section).

### <a name="execution_profiles">Execution profiles</a>

Resource settings are configured using Nextflow profiles. The ampliconseq
pipeline provides three profiles - `standard`, `bigserver` and `cluster`
configured for running on servers and the high-performance compute cluster at
CRUK CI. These specify the maximum number of CPUs or memory that can be used at
any one time during the pipeline run or the maximum number of jobs that can be
submitted to the cluster to be run in parallel.

A custom profile can be created in the configuration file, e.g.
`ampliconseq.config`, an example of which is shown below.

    profiles {
        myprofile {
            process.executor = 'local'
            executor {
                cpus = 8
                memory = 32.GB
            }
            singularity.enabled = true
        }
    }

The new profile, `myserver`, allows for up to 8 CPUs to be used at any one time
and a total of 32G of memory. The pipeline specifies how many CPU cores and how
much memory each process requires and Nextflow ensures that the overall resource
allocation does not exceed that specified in the profile.

The *local* executor (the default) runs processes on the computer on which
Nextflow is launched. The processes are parallelized by spawning multiple
threads and taking advantage of the multi-core architecture provided by the CPU.

Setting `singularity.enabled = true` in the profile tells Nextflow to use the
container with Singularity; it is not necessary to specify this separately with
the `-with-singularity` option.

This profile can be selected by using the `-profile` command line option.

    nextflow run crukci-bioinformatics/ampliconseq -config ampliconseq.config -profile myprofile

The following profile tells Nextflow to submit jobs to nodes on a compute
cluster using the Slurm resource manager. It allows for a maximum of 25 jobs to
be submitted to the 'long' queue for running in parallel and tells Nextflow to
poll every 30 seconds to check for completed jobs.

    profiles {
        mycluster {
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
    }

### <a name="nextflow_reports">Nextflow reports</a>

Nextflow can provide a useful summary report detailing the completion status,
execution time and memory used by each task, and a timeline chart.

Use the `-with-report` and `-with-timeline` command line options to produce
these reports, e.g.

     nextflow run crukci-bioinformatics/ampliconseq \
        -config ampliconseq.config \
        -with-report ampliconseq.report.html \
        -with-timeline ampliconseq.timeline.html

### <a name="nextflow_log_files_and_work_directories">Nextflow log files and work directories</a>

Nextflow logs information to a hidden file named `.nextflow.log` in the launch
directory in which ampliconseq is run. This contains logging information that
can help with debugging problems with a pipeline run. It will, for example, show
which task(s) failed and the directory in which that task was run. An
alternative log file name can be specified using the `-log` command line
argument.

    nextflow -log ampliconseq.log run crukci-bioinformatics/ampliconseq -config ampliconseq.config

`nextflow help` gives more details on command line options for Nextflow.

Nextflow runs each task within its own directory. These directories are created
under a work directory, by default a subdirectory of the launch directory named
`work` but which is configurable with the `-work-dir` command line option. Each
task directory contains hidden files with names such as `.command.sh` and
`.command.out`, inspection of which can be helpful when debugging pipeline runs.

Intermediate files created during a pipeline execution are written to the work
directories. The final outputs are written either to the launch directory or the
directory specified using the `--outputDir` command line option or the
`outputDir` parameter. The `work` directory (and all its subdirectories) can be
deleted on successful completion of the pipeline unless other Nextflow pipeline
runs are also making use of the same top-level work directory.

