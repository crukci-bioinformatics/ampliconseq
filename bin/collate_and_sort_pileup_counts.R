#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Collates and sorts pileup counts.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--sample"), dest = "sample",
              help = "Identifier for the sample used to populate a Sample column in the output table"),

  make_option(c("--amplicons"), dest = "amplicons_file",
              help = "Tab-separated file containing amplicon details the order of which is used in sorting the pileup counts (ID column required)"),

  make_option(c("--output"), dest = "output_file",
              help = "Output table")
)

option_parser <- OptionParser(usage = "usage: %prog [options] pileup_counts_files", option_list = option_list, add_help_option = TRUE)
arguments <- parse_args(option_parser, positional_arguments = TRUE)

input_files <- arguments$args
if (length(input_files) == 0) stop("Pileup counts file(s) must be specified")

opt <- arguments$options

sample <- opt$sample
amplicons_file <- opt$amplicons_file
output_file <- opt$output_file

if (is.null(amplicons_file)) stop("Amplicons file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

suppressPackageStartupMessages(library(tidyverse))

amplicons <- read_tsv(amplicons_file, col_types = cols(.default = "c"))

pileup_counts <- map_dfr(input_files, read_tsv, col_types = cols(Position = "i", .default = "c"))

pileup_counts <- pileup_counts %>%
  mutate(Sample = sample) %>%
  select(ID, Sample, everything()) %>%
  mutate(Amplicon = factor(Amplicon, levels = amplicons$ID)) %>%
  arrange(Amplicon, Position)

write_tsv(pileup_counts, output_file, na = "")

