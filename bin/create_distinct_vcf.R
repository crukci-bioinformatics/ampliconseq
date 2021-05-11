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

  make_option(c("--variants"), dest = "variants_file",
              help = "TSV file containing variants (Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--output"), dest = "output_file",
              help = "Output VCF file containing an entry for each distinct variant")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

variants_file <- opt$variants_file
output_file <- opt$output_file

if (is.null(variants_file)) stop("Input variants file must be specified")
if (is.null(output_file)) stop("Output VCF file must be specified")

suppressPackageStartupMessages(library(tidyverse))

variants <- read_tsv(variants_file)

expected_columns <- c("Chromosome", "Position", "Ref", "Alt")
missing_columns <- setdiff(expected_columns, colnames(variants))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", variants_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

# note that dplyr distinct function retains the same ordering for the subset of
# rows
distinct_variants <- variants %>%
  select(all_of(expected_columns)) %>%
  distinct()

distinct_variants %>%
  transmute(`#CHROM` = Chromosome, POS = Position, ID = NA, REF = Ref, ALT = Alt, QUAL = NA, FILTER = NA, INFO = NA) %>%
  write_tsv(output_file, na = ".")

