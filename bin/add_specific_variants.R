#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
#
# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Adds missing calls to the variant table so that entries specifying a 'no call'
# are included for libraries where the call was made in a replicate for the same
# sample or for variants listed for specific calling in the specific variants
# file.
# Specific variants are those of particular interest and which need to be
# included in the final variants table along with allele fractions, depths, etc.
# regardless of whether these are called or not.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--samples"), dest = "samples_file",
              help = "Sample sheet file (Sample column required)"),

  make_option(c("--called-variants"), dest = "variants_file",
              help = "Variant calls file (ID, Amplicon, Chromosome, Position, Ref, Alt and Filters columns required)"),

  make_option(c("--specific-variants"), dest = "specific_variants_file",
              help = "Specific/known variants file (Sample, Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--reference-sequence-index"), dest = "reference_sequence_index_file",
              help = "Index file for the reference genome sequence used for chromosome sort order (expected to have .fai extension)"),

  make_option(c("--output"), dest = "output_file",
              help = "Output file containing all variants including rows for missing calls")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

samples_file <- opt$samples_file
variants_file <- opt$variants_file
specific_variants_file <- opt$specific_variants_file
reference_sequence_index_file <- opt$reference_sequence_index_file
output_file <- opt$output_file

if (is.null(samples_file)) stop("Sample sheet file must be specified")
if (is.null(variants_file)) stop("Variant calls file must be specified")
if (is.null(specific_variants_file)) stop("Specific variants file must be specified")
if (is.null(reference_sequence_index_file)) stop("Reference sequence index file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

suppressPackageStartupMessages(library(tidyverse))

# read input files (no further checks needed)
samples <- read_tsv(samples_file, col_types = cols(.default = col_character()))
variants <- read_tsv(variants_file, col_types = cols(.default = col_character()))
specific_variants <- read_tsv(specific_variants_file, col_types = cols(.default = col_character()))

samples <- select(samples, Sample, ID)

sample_variant_columns <- c("Sample", "Amplicon", "Chromosome", "Position", "Ref", "Alt")

specific_variants <- specific_variants %>%
  select(all_of(sample_variant_columns)) %>%
  distinct() %>%
  mutate(Specific = "true")

# variants by sample
sample_variants <- samples %>%
  inner_join(variants, by = "ID") %>%
  select(all_of(sample_variant_columns)) %>%
  distinct()

# add specific variants
sample_variants <- sample_variants %>%
  full_join(specific_variants, by = sample_variant_columns) %>%
  mutate(Specific = replace_na(Specific, "false"))

# add variant calls
variants <- sample_variants %>%
  left_join(samples, by = "Sample") %>%
  left_join(variants, by = c("ID", "Amplicon", "Chromosome", "Position", "Ref", "Alt")) %>%
  select(Sample, ID, Amplicon, Chromosome, Position, Ref, Alt, everything()) %>%
  mutate(Filters = replace_na(Filters, "NOT_CALLED"))

# read reference genome index file
chromosomes <- read_tsv(reference_sequence_index_file, col_types = "cnnnn", col_names = c("Chromosome", "Length", "Offset", "Linebases", "Linewidth"))

# sort variants and write to output file
variants %>%
  mutate(Chromosome = factor(Chromosome, levels = chromosomes$Chromosome)) %>%
  arrange(Sample, Chromosome, Position, Ref, Alt, Amplicon, ID) %>%
  write_tsv(output_file, na = "")

