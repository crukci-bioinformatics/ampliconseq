#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Adds background noise thresholds to variants table and applies background
# noise filters.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--variants"), dest = "variants_file",
              help = "TSV file containing variants (ID, Amplicon, Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--position-thresholds"), dest = "position_thresholds_file",
              help = "TSV file containing position noise thresholds (Amplicon, Chromosome, Position, 'Reference base', 'Alternate allele' and 'Allele fraction threshold' columns required)"),

  make_option(c("--library-thresholds"), dest = "library_thresholds_file",
              help = "TSV file containing position noise thresholds (ID, 'Reference base', 'Alternate allele' and 'Allele fraction threshold' columns required)"),

  make_option(c("--output"), dest = "output_file",
              help = "Output variants file")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

variants_file <- opt$variants_file
position_thresholds_file <- opt$position_thresholds_file
library_thresholds_file <- opt$library_thresholds_file
output_file <- opt$output_file

if (is.null(variants_file)) stop("Input variant file must be specified")
if (is.null(position_thresholds_file)) stop("Position noise thresholds file must be specified")
if (is.null(library_thresholds_file)) stop("Library noise thresholds file must be specified")
if (is.null(output_file)) stop("Output VCF file must be specified")

suppressPackageStartupMessages(library(tidyverse))

variants <- read_tsv(variants_file, col_types = cols(`Allele fraction (pileup)` = "d", .default = "c"))

position_thresholds <- read_tsv(position_thresholds_file, col_types = cols(`Allele fraction threshold` = "d", .default = "c")) %>%
  rename(Ref = `Reference base`, Alt = `Alternate allele`, `Position noise threshold` = `Allele fraction threshold`)

library_thresholds <- read_tsv(library_thresholds_file, col_types = cols(`Allele fraction threshold` = "d", .default = "c")) %>%
  rename(Ref = `Reference base`, Alt = `Alternate allele`, `Library noise threshold` = `Allele fraction threshold`)

variants <- variants %>%
  left_join(position_thresholds, by = c("Amplicon", "Chromosome", "Position", "Ref", "Alt")) %>%
  mutate(Filters = ifelse(!is.na(`Allele fraction (pileup)`) & !is.na(`Position noise threshold`) & `Allele fraction (pileup)` > `Position noise threshold`, str_c(Filters, "position_noise", sep = ","), Filters)) %>%
  left_join(library_thresholds, by = c("ID", "Ref", "Alt")) %>%
  mutate(Filters = ifelse(!is.na(`Allele fraction (pileup)`) & !is.na(`Library noise threshold`) & `Allele fraction (pileup)` > `Library noise threshold`, str_c(Filters, "library_noise", sep = ","), Filters)) %>%
  mutate(Filters = str_remove(Filters, "^PASS,"))

write_tsv(variants, output_file, na = "")

