#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
#
# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Extract and collate most useful alignment and coverage metrics from Picard
# CollectAlignmentSummaryMetrics and CollectTargetedPcrMetrics, the summary from
# the ExtractAmpliconRegions utility and the pileup counts.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--alignment-metrics"), dest = "alignment_metrics_file",
              help = "Alignment metrics from Picard CollectAlignmentSummaryMetrics in collated tabular format"),

  make_option(c("--targeted-pcr-metrics"), dest = "targeted_pcr_metrics_file",
              help = "Targeted PCR metrics from Picard CollectTargetedPcrMetrics in collated tabular format"),

  make_option(c("--amplicon-coverage"), dest = "amplicon_coverage_file",
              help = "Amplicon coverage metrics file"),

  make_option(c("--pileup-counts"), dest = "pileup_counts_file",
              help = "Pileup counts file"),

  make_option(c("--output-metrics"), dest = "output_metrics_file",
              help = "Amplicon coverage metrics file")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

alignment_metrics_file <- opt$alignment_metrics_file
targeted_pcr_metrics_file <- opt$targeted_pcr_metrics_file
amplicon_coverage_file <- opt$amplicon_coverage_file
pileup_counts_file <- opt$pileup_counts_file
output_metrics_file <- opt$output_metrics_file

if (is.null(alignment_metrics_file)) stop("Alignment metrics file must be specified")
if (is.null(targeted_pcr_metrics_file)) stop("Targeted PCR metrics file must be specified")
if (is.null(amplicon_coverage_file)) stop("Alignment coverage file must be specified")
if (is.null(pileup_counts_file)) stop("Pileup counts file must be specified")
if (is.null(output_metrics_file)) stop("Output merged metrics file must be specified")

suppressPackageStartupMessages(library(tidyverse))

# Picard alignment summary metrics
alignment_metrics <- read_tsv(alignment_metrics_file, col_types = cols(.default = "c")) %>%
  filter(CATEGORY %in% c("PAIR", "UNPAIRED"))

alignment_metrics <- alignment_metrics %>%
  select(
    ID,
    Sample,
    Reads = PF_READS,
    MeanReadLength = MEAN_READ_LENGTH,
    ReadsAligned = PF_READS_ALIGNED,
    MismatchRate = PF_MISMATCH_RATE
  )

# Picard targeted PCR metrics
targeted_pcr_metrics <- read_tsv(targeted_pcr_metrics_file, col_types = cols(.default = "c")) %>%
  select(-any_of(c("SAMPLE", "LIBRARY", "READ_GROUP")))

if (nrow(targeted_pcr_metrics) != nrow(alignment_metrics)) {
  stop("unexpected number of entries in targeted PCR metrics file")
}

if (any(alignment_metrics$ID != targeted_pcr_metrics$ID)) {
  stop("inconsistent IDs in alignment metrics and targeted PCR metrics files")
}

# Picard TargetedPcrMetrics are affected by duplicate marking (unlike
# AlignmentSummaryMetrics which are not)
# Only using metrics that are unaffected, e.g. the number of bases aligned and
# the division of those between off-, near- and on-amplicon loci
targeted_pcr_metrics <- targeted_pcr_metrics %>%
  select(
    ID,
    Bases = PF_BASES,
    BasesAligned = PF_BASES_ALIGNED,
    BasesOnAmplicon = ON_AMPLICON_BASES,
    BasesNearAmplicon = NEAR_AMPLICON_BASES,
    BasesOffAmplicon = OFF_AMPLICON_BASES
  )

# PF_BASES_ALIGNED is slightly higher than PF_ALIGNED_BASES from
# CollectAlignmentSummaryMetrics - not sure why

# PF_BASES_ALIGNED excludes reads that fail the vendor's filter and secondary
# alignments, includes supplementary records

# PF_BASES_ALIGNED = ON_AMPLICON_BASES + NEAR_AMPLICON_BASES + OFF_AMPLICON_BASES

# NEAR_AMPLICON_BASES are bases that align outside an amplicon for a read
# that at least partially aligns to a wider region enclosing an amplicon (Picard
# uses a 'near distance' of 250 by default)

# Picard only includes really 'off-amplicon' bases (not near-amplicon bases)
# when computing the off amplicon rate
#   PCT_OFF_AMPLICON = OFF_AMPLICON_BASES / PF_BASES_ALIGNED

# amplicon coverage
amplicon_coverage <- read_tsv(amplicon_coverage_file, col_types = cols(Reads = "i", Bases = "i", .default = "c"))

# assigned reads and bases per library
assignment_metrics <- amplicon_coverage %>%
  group_by(ID) %>%
  summarize(ReadsAssigned = sum(Reads), BasesAssigned = sum(Bases))

if (nrow(assignment_metrics) != nrow(assignment_metrics)) {
  stop("unexpected number of entries in amplicon coverage file")
}

if (any(assignment_metrics$ID != alignment_metrics$ID)) {
  stop("inconsistent IDs in alignment metrics and amplicon coverage files")
}

# read pileup counts in chunks and collect numbers of bases on target and the
# usable (filtered) bases on target
pileup_metrics <- NULL

collect_pileup_metrics <- function(pileup_counts, pos) {
  pileup_metrics <<- pileup_metrics %>%
    bind_rows(
      pileup_counts %>%
        group_by(ID) %>%
        summarize(
          TargetPositions = n(),
          BasesOnTarget = sum(`Depth unfiltered`),
          BasesOnTargetUsable = sum(Depth)
        )
    )
}

result <- read_tsv_chunked(
  pileup_counts_file,
  SideEffectChunkCallback$new(collect_pileup_metrics),
  chunk_size = 100000,
  col_types = cols(Depth = "i", `Depth unfiltered` = "i", .default = "c"),
  progress = TRUE)

# sum up base and target position counts for libraries that were split across
# more than one chunk and compute the mean target coverage
pileup_metrics <- pileup_metrics %>%
  group_by(ID) %>%
  summarize(across(everything(), sum)) %>%
  mutate(MeanTargetCoverage = round(BasesOnTargetUsable / TargetPositions))

# merge metrics into single data frame
metrics <- alignment_metrics %>%
  left_join(targeted_pcr_metrics, by = "ID") %>%
  left_join(assignment_metrics, by = "ID") %>%
  left_join(pileup_metrics, by = "ID")

# write extracted metrics to output file
write_tsv(metrics, output_metrics_file, na = "")

