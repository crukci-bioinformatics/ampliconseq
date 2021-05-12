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
              help = "CSV/TSV file containing sample identifiers (Sample column required)"),

  make_option(c("--output"), dest = "output_file",
              help = "Output specific variants file in the format required for subsequent pipeline processes")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

input_file <- opt$input_file
samples_file <- opt$samples_file
output_file <- opt$output_file

if (is.null(input_file)) stop("Input file must be specified")
if (is.null(samples_file)) stop("Sample sheet file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

suppressPackageStartupMessages(library(tidyverse))

if (str_ends(str_to_lower(input_file), "\\.csv")) {
  specific_variants <- read_csv(input_file, col_types = cols(.default = col_character()))
} else {
  specific_variants <- read_tsv(input_file, col_types = cols(.default = col_character()))
}

# check for expected columns
expected_columns <- c("Sample", "Chromosome", "Position", "Ref", "Alt")
missing_columns <- setdiff(expected_columns, colnames(specific_variants))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", input_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

specific_variants <- select(specific_variants, all_of(expected_columns))

# check for missing values in any of the expected columns
missing_values <- filter(specific_variants, if_any(everything(), is.na))
if (nrow(missing_values) > 0) {
  stop("missing values found in ", input_file)
}

# check for multi-allelic variants
multiallelic_variants <- filter(specific_variants, if_any(Ref:Alt, ~ str_detect(.x, ",")))
if (nrow(multiallelic_variants) > 0) {
  stop(input_file, " should not contain multi-allelic variants")
}

# read sample sheet and check for unmatched sample identifiers
if (str_ends(str_to_lower(samples_file), "\\.csv")) {
  samples <- read_csv(samples_file, col_types = cols(.default = col_character()))
} else {
  samples <- read_tsv(samples_file, col_types = cols(.default = col_character()))
}

expected_columns <- "Sample"
missing_columns <- setdiff(expected_columns, colnames(samples))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", samples_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

unmatched_samples <- specific_variants %>%
  anti_join(samples, by = "Sample") %>%
  distinct(Sample)
if (nrow(unmatched_samples) > 0) {
  stop("found samples in ", input_file, " that don't appear in the sample sheet: ", str_c(unmatched_samples$Sample, collapse = ", "))
}

# write specific variants to CSV file
specific_variants %>%
  distinct() %>%
  write_csv(output_file)

