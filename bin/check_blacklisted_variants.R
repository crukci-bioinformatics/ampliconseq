#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
#
# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Checks the blacklisted variants file is valid and and outputs a cut-down
# version with only the expected columns.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--input"), dest = "input_file",
              help = "CSV/TSV file containing blacklisted variants (Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--output"), dest = "output_file",
              help = "Output blacklist file in the format required for subsequent pipeline processes")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

input_file <- opt$input_file
output_file <- opt$output_file

if (is.null(input_file)) stop("Input file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

suppressPackageStartupMessages(library(tidyverse))

if (str_ends(str_to_lower(input_file), "\\.csv")) {
  blacklisted_variants <- read_csv(input_file, col_types = cols(.default = col_character()))
} else {
  blacklisted_variants <- read_tsv(input_file, col_types = cols(.default = col_character()))
}

# check for expected columns
expected_columns <- c("Chromosome", "Position", "Ref", "Alt")
missing_columns <- setdiff(expected_columns, colnames(blacklisted_variants))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", input_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

blacklisted_variants <- select(blacklisted_variants, all_of(expected_columns))

# check for missing values in any of the expected columns
missing_values <- filter(blacklisted_variants, if_any(everything(), is.na))
if (nrow(missing_values) > 0) {
  stop("missing values found in ", input_file)
}

# check for multi-allelic variants
multiallelic_variants <- filter(blacklisted_variants, if_any(Ref:Alt, ~ str_detect(.x, ",")))
if (nrow(multiallelic_variants) > 0) {
  stop(input_file, " should not contain multi-allelic variants")
}

# write specific variants to CSV file
blacklisted_variants %>%
  distinct() %>%
  write_csv(output_file)

