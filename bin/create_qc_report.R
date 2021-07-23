#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
#
# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Create QC report from alignment coverage metrics and plots and results from
# the sample replicate correlation analysis based on variant allele fractions.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--alignment-metrics"), dest = "alignment_metrics_file",
              help = "Alignment and coverage metrics file"),

  make_option(c("--yield-plot"), dest = "yield_plot_file",
              help = "Yield stacked bar plot file"),

  make_option(c("--amplicon-coverage-plot"), dest = "amplicon_coverage_plot_file",
              help = "Amplicon coverage boxplot file"),

  make_option(c("--vaf-heatmap"), dest = "vaf_heatmap_file",
              help = "Variant allele fraction heatmap file"),

  make_option(c("--vaf-correlation-heatmap"), dest = "vaf_correlation_heatmap_file",
              help = "Variant allele fraction correlation heatmap file"),

  make_option(c("--replicate-mismatches"), dest = "replicate_mismatch_file",
              help = "Table of correlations between mismatched replicates"),

  make_option(c("--output-report"), dest = "output_report_file",
              help = "Output HTML report")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

alignment_metrics_file <- opt$alignment_metrics_file
yield_plot_file <- opt$yield_plot_file
amplicon_coverage_plot_file <- opt$amplicon_coverage_plot_file
vaf_heatmap_file <- opt$vaf_heatmap_file
vaf_correlation_heatmap_file <- opt$vaf_correlation_heatmap_file
replicate_mismatch_file <- opt$replicate_mismatch_file
output_report_file <- opt$output_report_file

if (is.null(alignment_metrics_file)) stop("Alignment metrics file must be specified")
if (is.null(yield_plot_file)) stop("Yield plot file must be specified")
if (is.null(amplicon_coverage_plot_file)) stop("Amplicon coverage plot file must be specified")
if (is.null(vaf_heatmap_file)) stop("VAF heatmap file must be specified")
if (is.null(vaf_correlation_heatmap_file)) stop("VAF correlation heatmap file must be specified")
if (is.null(replicate_mismatch_file)) stop("Replicate mismatches file must be specified")
if (is.null(output_report_file)) stop("Output report file must be specified")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Nozzle.R1)
  library(base64)
})

# function for creating base64-encoded image that can be embedded within the report
to_embedded_image <- function(file, mimetype)
{
  tempFile <- tempfile(pattern = "base64_", fileext = ".txt")
  encode(file, tempFile)
  encoded <- readChar(tempFile, file.info(tempFile)$size)
  file.remove(tempFile)
  encoded <- str_c("data:", mimetype, ";base64,", encoded)
  encoded
}

# Nozzle report

report <- newCustomReport("Amplicon Sequencing QC Report")
# report <- setReportSubTitle(report, title)

alignment_metrics_section <- newSection("Alignment Metrics")

alignment_metrics_section <- addTo(alignment_metrics_section, newParagraph("
Alignment metrics were computed using Picard CollectAlignmentSummaryMetrics and
percentages of bases on-, near- and off-amplicon were computed using Picard
CollectTargetedPcrMetrics.
"))

alignment_metrics_section <- addTo(alignment_metrics_section, newParagraph("
Pileup counts at target loci were generated with a custom Java utility written
using htsjdk (http://samtools.github.io/htsjdk). This tool is packaged as part
of the ampliconseq pipeline
(https://github.com/crukci-bioinformatics/ampliconseq).
"))

alignment_metrics_section <- addTo(alignment_metrics_section, newParagraph("
Usable reads exclude those that are unmapped or have a mapping quality of zero,
indicating that they align equally well to more than one location.
"))

alignment_metrics_section <- addTo(alignment_metrics_section, newParagraph("
Overlapping segments of reads from a read pair are not double counted in the
numbers of usable bases and mean target coverage values.
"))

alignment_metrics <- read_csv(alignment_metrics_file, col_types = cols(ID = "c", Sample = "c", .default = "d"))

alignment_metrics_table <- alignment_metrics %>%
  select(
    ID,
    Sample,
    Reads,
    `% reads mapped` = `% reads aligned`,
    `% bases mapped` = `% bases aligned`,
    `% mismatches`,
    `% bases on amplicon`,
    `% bases off amplicon`,
    `% reads assigned`,
    `% bases assigned`,
    `% bases on target`,
    `% usable bases on target` = `% bases on target (usable)`,
    `Mean target coverage`
  ) %>%
  as.data.frame() %>%
  newTable("Alignment Metrics", significantDigits = 3)

alignment_metrics_section <- addTo(alignment_metrics_section, alignment_metrics_table)

yield_section <- newSection("Sample Yields")

yield_section <- addTo(yield_section, newParagraph("
The yield is the total number of bases sequenced per sample.
"))

yield_section <- addTo(yield_section, newParagraph("
In the stacked bar chart below, the yield is divided into subsets depending on
where reads are mapped to the genome and whether these match and are assigned to
amplicon target regions.
"))

yield_section <- addTo(yield_section, newParagraph("
Reads considered to be near amplicon regions are those overlapping 250bp regions
flanking amplicon targets. Bases that are aligned within amplicon regions are
considered to be 'on target' if the read maps to the expected amplicon start or
end; how well a read or read pair must match an amplicon is configurable (see
the ampliconseq pipeline documentation for details).
"))

yield_section <- addTo(yield_section, newParagraph("
On-target bases that are from reads with zero mapping quality, i.e. could
align equally well elsewhere on the genome, or those with low base qualities are
considered unusable. Similarly, only one of two bases from overlapping segments
of reads from a read pair is counted within the usable on-target subset.
"))

yield_figure <- newFigure(
  to_embedded_image(yield_plot_file, "image/svg+xml"),
  "Yield for each sample/library."
)
yield_section <- addTo(yield_section, yield_figure)

amplicon_coverage_section <- newSection("Amplicon Coverage")

amplicon_coverage_section <- addTo(amplicon_coverage_section, newParagraph("
The mean coverage for each target amplicon is computed only for reads that have
been assigned to the amplicon based on one or both of the read ends aligning to
expected genomic locations for that amplicon. Overlapping segments of reads from
a read pair are not double counted.
"))

amplicon_coverage_figure <- newFigure(
  to_embedded_image(amplicon_coverage_plot_file, "image/svg+xml"),
  "Mean coverage for each target amplicon across all samples based on reads assigned to amplicons."
)
amplicon_coverage_section <- addTo(amplicon_coverage_section, amplicon_coverage_figure)

replicates_section <- newSection("Replicate libraries")

replicates_section <- addTo(replicates_section, newParagraph("
Variant allele fractions for single nucleotide substitutions were used to
cluster sample libraries. Variant loci were identified based on a configurable
minimum allele fraction within at least one library and the initial set of
variants were filtered to include those that could potentially help in
distinguishing between libraries - those with a narrow range of allele fractions
across all libraries or too many missing values (usually where the depth of
coverage is too low) were excluded.
"))

replicates_section <- addTo(replicates_section, newParagraph("
Some libraries may also have been excluded from the clustering because they have
too many missing allele fraction values among the set of substitutions used due
to low coverage.
"))

vaf_heatmap_figure <- newFigure(
  to_embedded_image(vaf_heatmap_file, "image/png"),
  "Heatmap of allele fractions for a set of single nucleotide substitutions. used to cluster sample libraries. Sample replicates are displayed with the same colour in the annotation bar."
)
replicates_section <- addTo(replicates_section, vaf_heatmap_figure)

vaf_correlation_heatmap_figure <- newFigure(
  to_embedded_image(vaf_correlation_heatmap_file, "image/png"),
  "Heatmap showing correlations between libraries of allele fractions for a set of single nucleotide substitutions. Sample replicates are displayed with the same colour in the annotation bar."
)
replicates_section <- addTo(replicates_section, vaf_correlation_heatmap_figure)

replicate_mismatches <- read_tsv(replicate_mismatch_file, col_types = cols(Correlation = "d", .default = "c"))

if (nrow(replicate_mismatches) == 0) {
  replicates_section <- addTo(replicates_section, newParagraph("
No mismatched replicate libraries detected.
"))
} else {
  replicate_mismatches_table <- replicate_mismatches %>%
    as.data.frame() %>%
    newTable("Correlations between mismatched libraries of allele fractions for single nucleotide substitutions.", significantDigits = 3)

  replicates_section <- addTo(replicates_section, replicate_mismatches_table)
}

report <- addTo(report, alignment_metrics_section)
report <- addTo(report, yield_section)
report <- addTo(report, amplicon_coverage_section)
report <- addTo(report, replicates_section)

temp_report_file <- tempfile(pattern = "report_", tmpdir = getwd(), fileext = "")
writeReport(report, filename = temp_report_file, output = HTML.REPORT)
file.rename(str_c(temp_report_file, ".html"), output_report_file)

