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

  make_option(c("--samples"), dest = "samples_file",
              help = "CSV/TSV file containing sample identifiers (ID and Sample columns required)"),

  make_option(c("--amplicons"), dest = "amplicons_file",
              help = "Tab-separated file containing amplicon details the order of which is used in sorting the amplicon coverage metrics (ID column required)"),

  make_option(c("--alignment-metrics"), dest = "alignment_metrics_file",
              help = "Alignment metrics from Picard CollectAlignmentSummaryMetrics in collated tabular format"),

  make_option(c("--targeted-pcr-metrics"), dest = "targeted_pcr_metrics_file",
              help = "Targeted PCR metrics from Picard CollectTargetedPcrMetrics in collated tabular format"),

  make_option(c("--amplicon-coverage"), dest = "amplicon_coverage_file",
              help = "Amplicon coverage metrics file"),

  make_option(c("--pileup-counts"), dest = "pileup_counts_file",
              help = "Pileup counts file"),

  make_option(c("--alignment-coverage-metrics"), dest = "alignment_coverage_metrics_file",
              help = "Output alignment coverage metrics CSV file"),

  make_option(c("--sorted-amplicon-coverage"), dest = "sorted_amplicon_coverage_file",
              help = "Sorted amplicon coverage table"),

  make_option(c("--amplicon-coverage-metrics"), dest = "amplicon_coverage_metrics_file",
              help = "Output amplicon coverage metrics CSV file")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

samples_file <- opt$samples_file
amplicons_file <- opt$amplicons_file
alignment_metrics_file <- opt$alignment_metrics_file
targeted_pcr_metrics_file <- opt$targeted_pcr_metrics_file
amplicon_coverage_file <- opt$amplicon_coverage_file
pileup_counts_file <- opt$pileup_counts_file
alignment_coverage_metrics_file <- opt$alignment_coverage_metrics_file
sorted_amplicon_coverage_file <- opt$sorted_amplicon_coverage_file
amplicon_coverage_metrics_file <- opt$amplicon_coverage_metrics_file

if (is.null(samples_file)) stop("Sample sheet file must be specified")
if (is.null(amplicons_file)) stop("Amplicon coordinates file must be specified")
if (is.null(alignment_metrics_file)) stop("Alignment metrics file must be specified")
if (is.null(targeted_pcr_metrics_file)) stop("Targeted PCR metrics file must be specified")
if (is.null(amplicon_coverage_file)) stop("Amplicon coverage file must be specified")
if (is.null(pileup_counts_file)) stop("Pileup counts file must be specified")
if (is.null(alignment_coverage_metrics_file)) stop("Output alignment coverage metrics file must be specified")
if (is.null(sorted_amplicon_coverage_file)) stop("Sorted amplicon coverage output file must be specified")
if (is.null(amplicon_coverage_metrics_file)) stop("Output amplicon coverage metrics file must be specified")

suppressPackageStartupMessages(library(tidyverse))


# sample sheet
samples <- read_tsv(samples_file, col_types = cols(ID = "f", .default = "c"))
samples <- select(samples, ID, Sample)
ids <- levels(samples$ID)


# Picard alignment summary metrics
alignment_metrics <- read_tsv(
  alignment_metrics_file,
  col_types = cols(
    PF_READS = "i",
    MEAN_READ_LENGTH = "d",
    PF_READS_ALIGNED = "i",
    PF_MISMATCH_RATE = "d",
    .default = "c")
  ) %>%
  filter(CATEGORY %in% c("PAIR", "UNPAIRED")) %>%
  transmute(
    ID,
    Reads = PF_READS,
    `Mean read length` = round(MEAN_READ_LENGTH, digits = 1),
    `Reads aligned` = PF_READS_ALIGNED,
    `% reads aligned` = round(100.0 * PF_READS_ALIGNED / PF_READS, digits = 2),
    `% mismatches` = round(100.0 * PF_MISMATCH_RATE, digits = 3)
  )


# Picard targeted PCR metrics
# Note reading number of bases in as doubles as can exceed maximum size for an
# integer
targeted_pcr_metrics <- read_tsv(
  targeted_pcr_metrics_file,
  col_types = cols(
    PF_BASES = "d",
    PF_BASES_ALIGNED = "d",
    ON_AMPLICON_BASES = "d",
    NEAR_AMPLICON_BASES = "d",
    OFF_AMPLICON_BASES = "d",
    .default = "c")
  ) %>%
  transmute(
    ID,
    Bases = PF_BASES,
    `Bases aligned` = PF_BASES_ALIGNED,
    `% bases aligned` = round(100.0 * PF_BASES_ALIGNED / PF_BASES, digits = 2),
    `Bases on amplicon` = ON_AMPLICON_BASES,
    `% bases on amplicon` = round(100.0 * ON_AMPLICON_BASES / PF_BASES, digits = 2),
    `Bases near amplicon` = NEAR_AMPLICON_BASES,
    `% bases near amplicon` = round(100.0 * NEAR_AMPLICON_BASES / PF_BASES, digits = 2),
    `Bases off amplicon` = OFF_AMPLICON_BASES,
    `% bases off amplicon` = round(100.0 * OFF_AMPLICON_BASES / PF_BASES, digits = 2)
  )

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


# merge metrics into single data frame
metrics <- left_join(alignment_metrics, targeted_pcr_metrics, by = "ID")

# add sample identifiers and sort into order given in sample sheet
metrics <- metrics %>%
  mutate(ID = factor(ID, levels = ids)) %>%
  left_join(samples, by = "ID") %>%
  select(ID, Sample, everything()) %>%
  arrange(ID)


# amplicon coverage
amplicon_coverage <- read_tsv(amplicon_coverage_file, col_types = cols(Reads = "i", Bases = "i", .default = "c"))

# amplicon coordinates files
amplicons <- read_tsv(amplicons_file, col_types = cols(.default = "c"))

# sort by ID and amplicon
amplicon_coverage <- amplicon_coverage %>%
  mutate(ID = factor(ID, levels = ids)) %>%
  mutate(Amplicon = factor(Amplicon, levels = amplicons$ID)) %>%
  arrange(ID, Amplicon)

# write sorted amplicon coverage table
write_tsv(amplicon_coverage, sorted_amplicon_coverage_file, na = "")

# add sample identifiers
amplicon_coverage <- amplicon_coverage %>%
  left_join(samples, by = "ID") %>%
  select(ID, Sample, everything())

# assigned reads and bases per library
assignment_metrics <- amplicon_coverage %>%
  group_by(ID) %>%
  summarize(`Reads assigned` = sum(Reads), `Bases assigned` = sum(Bases)) %>%
  arrange(ID)

if (nrow(assignment_metrics) != nrow(assignment_metrics)) {
  stop("unexpected number of entries in amplicon coverage file")
}

if (any(assignment_metrics$ID != metrics$ID)) {
  stop("inconsistent IDs in alignment metrics and amplicon coverage files")
}

# compute percentages
assignment_metrics <- assignment_metrics %>%
  left_join(select(metrics, ID, Reads, Bases), by = "ID") %>%
  mutate(
    `% reads assigned` = round(100.0 * `Reads assigned` / Reads, digits = 2),
    `% bases assigned` = round(100.0 * `Bases assigned` / Bases, digits = 2)
  ) %>%
  select(
    ID,
    `Reads assigned`,
    `% reads assigned`,
    `Bases assigned`,
    `% bases assigned`
  )

# merge metrics into single data frame
metrics <- left_join(metrics, assignment_metrics, by = "ID")


# read pileup counts in chunks and collect numbers of bases on target and the
# usable (filtered) bases on target
pileup_metrics <- NULL

collect_pileup_metrics <- function(pileup_counts, pos) {
  pileup_metrics <<- pileup_metrics %>%
    bind_rows(
      pileup_counts %>%
        group_by(ID) %>%
        summarize(
          `Target positions` = n(),
          `Bases on target` = sum(`Depth unfiltered`),
          `Bases on target (usable)` = sum(Depth)
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
# more than one chunk
pileup_metrics <- pileup_metrics %>%
  group_by(ID) %>%
  summarize(across(everything(), sum))

# compute percentages and mean target coverage
pileup_metrics <- pileup_metrics %>%
  left_join(select(metrics, ID, Bases), by = "ID") %>%
  mutate(
    `% bases on target` = round(100.0 * `Bases on target` / Bases, digits = 2),
    `% bases on target (usable)` = round(100.0 * `Bases on target (usable)` / Bases, digits = 2),
    `Mean target coverage` = round(`Bases on target (usable)` / `Target positions`)
  ) %>%
  select(
    ID,
    `Bases on target`,
    `% bases on target`,
    `Bases on target (usable)`,
    `% bases on target (usable)`,
    `Mean target coverage`
  )

# merge metrics into single data frame
metrics <- left_join(metrics, pileup_metrics, by = "ID")


# write merged alignment coverage metrics to output file
write_csv(metrics, alignment_coverage_metrics_file, na = "")


# write amplicon coverage metrics output file
amplicon_coverage %>%
  select(ID, Sample, Amplicon, Reads, `Read pairs`, `Mean coverage`) %>%
  write_csv(amplicon_coverage_metrics_file, na = "")

