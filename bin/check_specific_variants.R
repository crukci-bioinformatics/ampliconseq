#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
#
# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Checks the specific variants file is valid and and outputs a cut-down
# version with only the expected columns.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--input"), dest = "input_file",
              help = "CSV/TSV file containing specific/known variants (Sample, Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--samples"), dest = "samples_file",
              help = "CSV/TSV file containing sample identifiers (ID and Sample columns required)"),

  make_option(c("--amplicons"), dest = "amplicons_file",
              help = "CSV/TSV file containing details of the amplicons (ID, Chromosome, TargetStart and TargetEnd columns required)"),

  make_option(c("--output"), dest = "output_file",
              help = "Output specific variants file in the format required for subsequent pipeline processes")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

input_file <- opt$input_file
samples_file <- opt$samples_file
amplicons_file <- opt$amplicons_file
output_file <- opt$output_file

if (is.null(input_file)) stop("Input file must be specified")
if (is.null(samples_file)) stop("Sample sheet file must be specified")
if (is.null(amplicons_file)) stop("Amplicon details file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

suppressPackageStartupMessages(library(tidyverse))

# read specific variants
if (str_ends(str_to_lower(input_file), "\\.csv")) {
  variants <- read_csv(input_file, col_types = cols(.default = col_character()))
} else {
  variants <- read_tsv(input_file, col_types = cols(.default = col_character()))
}

# check for expected columns
expected_columns <- c("Sample", "Chromosome", "Position", "Ref", "Alt")
missing_columns <- setdiff(expected_columns, colnames(variants))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", input_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

variants <- select(variants, all_of(expected_columns))

# check for missing values in any of the expected columns
missing_values <- filter(variants, if_any(everything(), is.na))
if (nrow(missing_values) > 0) {
  stop("missing values found in ", input_file)
}

# convert positions to integer values
variants <- mutate(variants, Position = parse_integer(Position))
missing_values <- filter(variants, is.na(Position))
if (nrow(missing_values) > 0) {
  stop("variants with non-integer coordinates found in ", variants_file)
}

# check for multi-allelic variants
multiallelic_variants <- filter(variants, if_any(Ref:Alt, ~ str_detect(.x, ",")))
if (nrow(multiallelic_variants) > 0) {
  stop(input_file, " should not contain multi-allelic variants")
}

# read sample sheet and check for unmatched sample identifiers
samples <- read_tsv(samples_file, col_types = cols(.default = "c"))
samples <- select(samples, ID, Sample)

unmatched_samples <- variants %>%
  anti_join(samples, by = "Sample") %>%
  distinct(Sample)
if (nrow(unmatched_samples) > 0) {
  stop("found samples in ", input_file, " that don't appear in the sample sheet: ", str_c(unmatched_samples$Sample, collapse = ", "))
}

# read and check amplicon details file
if (str_ends(str_to_lower(amplicons_file), "\\.csv")) {
  amplicons <- read_csv(amplicons_file, col_types = cols(.default = col_character()))
} else {
  amplicons <- read_tsv(amplicons_file, col_types = cols(.default = col_character()))
}

expected_columns <- c("ID", "Chromosome", "TargetStart", "TargetEnd")
missing_columns <- setdiff(expected_columns, colnames(amplicons))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", amplicons_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}
amplicons <- select(amplicons, all_of(expected_columns))

missing_values <- filter(amplicons, if_any(everything(), is.na))
if (nrow(missing_values) > 0) {
  stop("amplicons with missing values found in ", amplicons_file, ": '", str_c(missing_values$ID, collapse = "', '"), "'")
}

# convert target coordinates to integer values
amplicons <- mutate(amplicons, across(TargetStart:TargetEnd, parse_integer))
missing_values <- filter(amplicons, if_any(TargetStart:TargetEnd, is.na))
if (nrow(missing_values) > 0) {
  stop("amplicons with non-integer coordinates found in ", amplicons_file, ": '", str_c(missing_values$ID, collapse = "', '"), "'")
}

# assign amplicons to specific variants
amplicons <- rename(amplicons, Amplicon = ID)

variants_within_amplicons <- variants %>%
  left_join(amplicons, by = "Chromosome") %>%
  filter(Position >= TargetStart, Position <= TargetEnd) %>%
  select(!TargetStart:TargetEnd)

unmatched_amplicons <- anti_join(variants, variants_within_amplicons, by = colnames(variants))
if (nrow(unmatched_amplicons) > 0) {
  stop("found variants in ", input_file, " that don't aren't within amplicon target regions")
}

# add library/dataset ID and write to TSV file
samples %>%
  inner_join(variants_within_amplicons, by = "Sample") %>%
  select(Sample, ID, Amplicon, everything()) %>%
  distinct() %>%
  write_tsv(output_file)

