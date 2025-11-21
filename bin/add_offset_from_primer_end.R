#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2025 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Annotate variants with the offset from the nearest primer end.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--variants"), dest = "variants_file",
              help = "TSV file containing variants (Amplicon, Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--amplicons"), dest = "amplicons_file",
              help = "TSV file containing details of the amplicons (ID, Chromosome, TargetStart and TargetEnd columns required)"),

  make_option(c("--output"), dest = "output_file",
              help = "Output TSV file containing annotated variants")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

variants_file <- opt$variants_file
amplicons_file <- opt$amplicons_file
output_file <- opt$output_file

if (is.null(variants_file)) stop("Input variant file must be specified")
if (is.null(amplicons_file)) stop("Amplicon details file must be specified")
if (is.null(output_file)) stop("Output annotated variant file must be specified")

suppressPackageStartupMessages(library(tidyverse))

# read variants
variants <- read_tsv(variants_file, col_types = cols(
  Position = "n",
  .default = "c"))

expected_columns <- c("Amplicon", "Chromosome", "Position", "Ref", "Alt")
missing_columns <- setdiff(expected_columns, colnames(variants))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", variants_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

# read and check amplicon details file
amplicons <- read_tsv(amplicons_file, col_types = cols(
  TargetStart = "n",
  TargetEnd = "n",
  .default = "c"))

expected_columns <- c("ID", "Chromosome", "TargetStart", "TargetEnd")
missing_columns <- setdiff(expected_columns, colnames(amplicons))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", amplicons_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}
amplicons <- select(amplicons, all_of(expected_columns))

# annotate with offset/distance from nearest primer end
# assume left-aligned and normalized indels
variants <- variants %>%
  left_join(amplicons, by = c("Amplicon" = "ID", "Chromosome")) %>%
  mutate(RefLength = str_length(Ref), AltLength = str_length(Alt)) %>%
  mutate(Start = Position, End = Position + RefLength - 1) %>%
  mutate(Insertion = AltLength > RefLength & str_sub(Alt, 1, RefLength) == Ref) %>%
  mutate(Deletion = RefLength > AltLength & str_sub(Ref, 1, AltLength) == Alt) %>%
  mutate(End = ifelse(Deletion, Position + RefLength - 1, End)) %>%
  mutate(Start = ifelse(Deletion, End - RefLength + AltLength + 1, Start)) %>%
  mutate(Start = ifelse(Insertion, Position + RefLength, Start)) %>%
  mutate(End = ifelse(Insertion, Start - 1, End)) %>%
  mutate(`Offset from primer end` = pmax(pmin(Start - TargetStart, TargetEnd - End) + 1, 0)) %>%
  select(-TargetStart, -TargetEnd, -RefLength, -AltLength, -Start, -End, -Insertion, -Deletion)

# write annotated variants to TSV file
write_tsv(variants, output_file, na = "")

