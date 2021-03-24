#!/usr/bin/env Rscript

# Checks the sample sheet file is valid, i.e. has the expected columns, etc.,
# and outputs a cut-down version with only the expected columns.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
{
  message("Usage: Rscript check_samples_file.R samples_file output_samples_file")
  quit(status = 1)
}

samples_file <- args[1]
output_file <- args[2]

suppressPackageStartupMessages(library(tidyverse))

# read and check sample sheet
if (str_ends(str_to_lower(samples_file), "\\.csv")) {
  samples <- read_csv(samples_file, col_types = cols(.default = col_character()))
} else {
  samples <- read_tsv(samples_file, col_types = cols(.default = col_character()))
}

expected_columns <- c("ID", "Sample")
missing_columns <- setdiff(expected_columns, colnames(samples))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", samples_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

if (nrow(samples) == 0) {
  stop("empty samples file: ", samples_file)
}

if (nrow(filter(samples, is.na(samples$ID))) > 0) {
  stop("missing IDs found in ", samples_file)
}

if (nrow(filter(samples, is.na(samples$Sample))) > 0) {
  stop("missing Sample names found in ", samples_file)
}

duplicates <- samples %>%
  count(ID) %>%
  filter(n > 1) %>%
  pull(ID)
if (length(duplicates) > 0) {
  stop("duplicate IDs found in ", samples_file, ": '", str_c(duplicates, collapse = "', '"), "'")
}

# write new samples file containing only the expected columns
samples %>%
  select(all_of(expected_columns)) %>%
  write_csv(output_file)

