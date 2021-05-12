#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
#
# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Checks the sample sheet file is valid, i.e. has the expected columns, etc.,
# and outputs a cut-down version with only the expected columns.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--input"), dest = "input_file",
              help = "CSV/TSV file containing details of sample datasets (ID and Sample columns required)"),

  make_option(c("--output"), dest = "output_file",
              help = "Output sample sheet file in the format required for subsequent pipeline processes")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

input_file <- opt$input_file
output_file <- opt$output_file

if (is.null(input_file)) stop("Sample sheet file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

suppressPackageStartupMessages(library(tidyverse))

# read and check sample sheet
if (str_ends(str_to_lower(input_file), "\\.csv")) {
  samples <- read_csv(input_file, col_types = cols(.default = col_character()))
} else {
  samples <- read_tsv(input_file, col_types = cols(.default = col_character()))
}

expected_columns <- c("ID", "Sample")
missing_columns <- setdiff(expected_columns, colnames(samples))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", input_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

if (nrow(samples) == 0) {
  stop("empty samples file: ", input_file)
}

if (nrow(filter(samples, is.na(samples$ID))) > 0) {
  stop("missing IDs found in ", input_file)
}

if (nrow(filter(samples, is.na(samples$Sample))) > 0) {
  stop("missing Sample names found in ", input_file)
}

duplicates <- samples %>%
  count(ID) %>%
  filter(n > 1) %>%
  pull(ID)
if (length(duplicates) > 0) {
  stop("duplicate IDs found in ", input_file, ": '", str_c(duplicates, collapse = "', '"), "'")
}

# write new samples file containing only the expected columns
samples %>%
  select(all_of(expected_columns)) %>%
  write_csv(output_file)
