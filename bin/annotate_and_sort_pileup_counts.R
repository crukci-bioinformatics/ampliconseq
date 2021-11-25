#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Collates and sorts pileup counts and annotates with the sample name.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--pileup-counts"), dest = "pileup_counts_file",
              help = "Pileup counts file"),

  make_option(c("--samples"), dest = "samples_file",
              help = "Sample sheet file (ID and Sample columns required)"),

  make_option(c("--amplicons"), dest = "amplicons_file",
              help = "Tab-separated file containing amplicon details the order of which is used in sorting the pileup counts (ID column required)"),

  make_option(c("--output"), dest = "output_file",
              help = "Output table")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

pileup_counts_file <- opt$pileup_counts_file
samples_file <- opt$samples_file
amplicons_file <- opt$amplicons_file
output_file <- opt$output_file

if (is.null(pileup_counts_file)) stop("Pileup counts file must be specified")
if (is.null(samples_file)) stop("Samples file must be specified")
if (is.null(amplicons_file)) stop("Amplicons file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

suppressPackageStartupMessages(library(tidyverse))

pileup_counts <- read_tsv(pileup_counts_file, col_types = cols(Position = "i", .default = "c"))

samples <- read_tsv(samples_file, col_types = cols(.default = "c"))
samples <- select(samples, ID, Sample)

amplicons <- read_tsv(amplicons_file, col_types = cols(.default = "c"))

pileup_counts <- pileup_counts %>%
  left_join(samples, by = "ID") %>%
  select(ID, Sample, everything()) %>%
  mutate(Amplicon = factor(Amplicon, levels = amplicons$ID)) %>%
  arrange(Amplicon, Position)

write_tsv(pileup_counts, output_file, na = "")

