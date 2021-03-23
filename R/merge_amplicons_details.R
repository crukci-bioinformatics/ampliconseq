#!/usr/bin/env Rscript

# Creates the amplicon details file required by the new version of the
# ampliconseq pipeline from the amplicon and target intervals files used in
# previous versions and the previous amplion details file (containing ID and
# gene columns).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4)
{
  stop("Usage: Rscript merge_amplicons_details.R amplicon_intervals_file target_intervals_file amplicon_details_file output_csv_file")
}

amplicon_intervals_file <- args[1]
target_intervals_file <- args[2]
amplicon_details_file <- args[3]
output_file <- args[4]

suppressPackageStartupMessages(library(tidyverse))


# function to read and check contents of an intervals file
read_intervals <- function(intervals_file) {
  intervals <- read_tsv(intervals_file, comment = "@", col_names = c("Chromosome", "Start", "End", "Strand", "ID"), col_types = "ciicc") %>%
    select(-Strand)

  if (nrow(intervals) == 0) {
    stop("no intervals found in ", intervals_file)
  }

  if (nrow(filter(intervals, is.na(intervals$ID))) > 0) {
    stop("missing IDs in ", intervals_file)
  }

  duplicates <- intervals %>%
    count(ID) %>%
    filter(n > 1)
  if (nrow(duplicates) > 0) {
    stop("duplicate IDs found in ", intervals_file, ": '", str_c(duplicates$ID, collapse = "', '"), "'")
  }

  incomplete <- filter(intervals, if_any(everything(), ~ is.na(.x)))
  if (nrow(incomplete) > 0) {
    stop("missing values found in ", intervals_file, " for following amplicons: '", str_c(incomplete$ID, collapse = "', '"), "'")
  }

  incorrect_start_and_end <- filter(intervals, Start < 1 | Start > End)
  if (nrow(incorrect_start_and_end)) {
    stop("invalid start and/or end coordinates in ", intervals_file, " for following amplicons: '", str_c(incorrect_start_and_end$ID, collapse = "', '"), "'")
  }

  duplicates <- intervals %>%
    count(Chromosome, Start, End) %>%
    filter(n > 1) %>%
    left_join(intervals, by = c("Chromosome", "Start", "End"))
  if (nrow(duplicates) > 0) {
    stop("amplicons with same genomic coordinates in ", intervals_file, ": '", str_c(duplicates$ID, collapse = "', '"), "'")
  }

  intervals
}


# read amplicon and target intervals files

amplicons <- read_intervals(amplicon_intervals_file) %>%
  rename(AmpliconStart = Start, AmpliconEnd = End)

targets <- read_intervals(target_intervals_file) %>%
  rename(TargetChromosome = Chromosome, TargetStart = Start, TargetEnd = End)


# check for inconsistencies between target and amplicon intervals

if(nrow(amplicons) != nrow(targets)) {
  stop(amplicon_intervals_file, " and ", target_intervals_file, " have different numbers of entries")
}

missing <- anti_join(targets, amplicons, by = "ID")
if (nrow(missing) > 0) {
  stop("IDs found in ", target_intervals_file, " without corresponding entry in ", amplicon_intervals_file, ": '", str_c(missing$ID, collapse = "', '"), "'")
}

missing <- anti_join(amplicons, targets, by = "ID")
if (nrow(missing) > 0) {
  stop("IDs found in ", amplicon_intervals_file, " without corresponding entry in ", target_intervals_file, ": '", str_c(missing$ID, collapse = "', '"), "'")
}


# join amplicon and target intervals into single table

amplicons <- amplicons %>%
  left_join(targets, by = "ID") %>%
  select(ID, Chromosome, AmpliconStart, AmpliconEnd, TargetChromosome, TargetStart, TargetEnd)


# further checks for inconsistencies in amplicon coordinates

mismatched_chromosomes <- filter(amplicons, Chromosome != TargetChromosome)
if (nrow(mismatched_chromosomes) > 0) {
  stop("mismatched chromosomes in ", amplicon_intervals_file, " and ", target_intervals_file, " for the following amplicons: '", str_c(mismatched_chromosomes$ID, collapse = "', '"), "'")
}

targets_outside_amplicons <- filter(amplicons, TargetStart < AmpliconStart | TargetEnd > AmpliconEnd)
if (nrow(targets_outside_amplicons) > 0) {
  stop("Target intervals not within amplicon intervals: '", str_c(targets_outside_amplicons$ID, collapse = "', '"), "'")
}


# read and check amplicon details file

amplicon_details <- read_tsv(amplicon_details_file, col_types = cols(.default = col_character()))

expected_columns <- c("ID", "Gene")
missing_columns <- setdiff(expected_columns, colnames(amplicon_details))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", amplicon_details_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

amplicon_details <- select(amplicon_details, all_of(expected_columns))

unexpected_ids <- anti_join(amplicon_details, amplicons, by = "ID")
if (nrow(unexpected_ids) > 0) {
  stop("IDs found in ", amplicon_details_file, " without corresponding entry in ", amplicon_intervals_file, ": '", str_c(missing$ID, collapse = "', '"), "'")
}


# add amplicon details and write to output file

amplicons %>%
  select(-TargetChromosome) %>%
  left_join(amplicon_details, by = "ID") %>%
  write_csv(output_file, na = "")


