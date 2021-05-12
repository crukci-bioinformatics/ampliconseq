#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Computes the allele fraction from the pileup counts and adds this and the
# pileup depth as independently-computed values to the variants table

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--variants"), dest = "variants_file",
              help = "TSV file containing variants (ID, Amplicon, Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--pileup-counts"), dest = "pileup_counts_file",
              help = "TSV file containing pileup counts (ID, Amplicon, Chromosome, Position, 'Reference base', Depth, 'A count', 'C count', 'G count' and 'T count' columns required)"),

  make_option(c("--output"), dest = "output_file",
              help = "Output variants file including additional columns for depth and allele fraction computed from the pileup counts")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

variants_file <- opt$variants_file
pileup_counts_file <- opt$pileup_counts_file
output_file <- opt$output_file

if (is.null(variants_file)) stop("Input variant file must be specified")
if (is.null(pileup_counts_file)) stop("Pileup counts file must be specified")
if (is.null(output_file)) stop("Output VCF file must be specified")

suppressPackageStartupMessages(library(tidyverse))

variants <- read_tsv(variants_file, col_types = cols(.default = "c"))

pileup_counts <- read_tsv(pileup_counts_file, col_types = cols(.default = "c"))

pileup_counts <- pileup_counts %>%
  semi_join(variants, by = c("ID", "Amplicon", "Chromosome", "Position")) %>%
  select(ID, Amplicon, Chromosome, Position, Ref = `Reference base`, `A count`:`T count`, Depth) %>%
  pivot_longer(`A count`:`T count`, names_to = "Alt", values_to = "Count") %>%
  mutate(Alt = str_remove(Alt, " count$")) %>%
  mutate(across(c(Count, Depth), parse_integer)) %>%
  mutate(`Allele fraction` = Count / Depth) %>%
  transmute(ID, Amplicon, Chromosome, Position, Ref, Alt, `Depth (pileup)` = Depth, `Allele fraction (pileup)` = Count / Depth)

variants %>%
  left_join(pileup_counts, by = c("ID", "Amplicon", "Chromosome", "Position", "Ref", "Alt")) %>%
  write_tsv(output_file, na = "")

