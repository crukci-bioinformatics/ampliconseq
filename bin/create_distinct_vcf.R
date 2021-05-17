#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Creates a VCF file containing the unique set of variants from the given input
# table, e.g. to be used as input to Ensembl Variant Effect Predictor.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--input"), dest = "input_file",
              help = "CSV/TSV file containing variants (Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--reference-sequence-index"), dest = "reference_sequence_index_file",
              help = "Index file for the reference genome sequence used for chromosome sort order (expected to have .fai extension)"),

  make_option(c("--output"), dest = "output_file",
              help = "Output VCF file containing an entry for each distinct variant")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

input_file <- opt$input_file
reference_sequence_index_file <- opt$reference_sequence_index_file
output_file <- opt$output_file

if (is.null(input_file)) stop("Input variant file must be specified")
if (is.null(reference_sequence_index_file)) stop("Reference sequence index file must be specified")
if (is.null(output_file)) stop("Output VCF file must be specified")

suppressPackageStartupMessages(library(tidyverse))

if (str_ends(str_to_lower(input_file), "\\.csv")) {
  variants <- read_csv(input_file, col_types = cols(.default = col_character()))
} else {
  variants <- read_tsv(input_file, col_types = cols(.default = col_character()))
}

expected_columns <- c("Chromosome", "Position", "Ref", "Alt")
missing_columns <- setdiff(expected_columns, colnames(variants))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", input_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

variants <- select(variants, all_of(expected_columns))

missing_values <- filter(variants, if_any(everything(), is.na))
if (nrow(missing_values) > 0) {
  stop("missing values found in ", input_file)
}

# convert positions to integer values
variants <- mutate(variants, Position = parse_integer(Position))
missing_values <- filter(variants, is.na(Position))
if (nrow(missing_values) > 0) {
  stop("variants with non-integer coordinates found in ", input_file)
}

# read reference genome index file
chromosomes <- read_tsv(reference_sequence_index_file, col_types = "ciiii", col_names = c("Chromosome", "Length", "Offset", "Linebases", "Linewidth"))

# get unique set of variants in coordinate sorted order
distinct_variants <- variants %>%
  distinct(Chromosome, Position, Ref, Alt) %>%
  mutate(Chromosome = factor(Chromosome, levels = chromosomes$Chromosome)) %>%
  arrange(Chromosome, Position, Ref, Alt)

# write output VCF
write_lines("##fileformat=VCFv4.2", output_file)

distinct_variants %>%
  transmute(`#CHROM` = Chromosome, POS = Position, ID = NA, REF = Ref, ALT = Alt, QUAL = NA, FILTER = NA, INFO = NA) %>%
  write_tsv(output_file, append = TRUE, col_names = TRUE, na = ".")

