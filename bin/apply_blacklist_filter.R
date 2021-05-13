#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Apply filter for blacklisted variants

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--variants"), dest = "variants_file",
              help = "TSV file containing variants (Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--blacklist"), dest = "blacklisted_variants_file",
              help = "TSV file containing blacklisted variants (Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--output"), dest = "output_file",
              help = "Output variants file")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

variants_file <- opt$variants_file
blacklisted_variants_file <- opt$blacklisted_variants_file
output_file <- opt$output_file

if (is.null(variants_file)) stop("Input variant file must be specified")
if (is.null(blacklisted_variants_file)) stop("Blacklisted variants file must be specified")
if (is.null(output_file)) stop("Output VCF file must be specified")

suppressPackageStartupMessages(library(tidyverse))

variants <- read_tsv(variants_file, col_types = cols(.default = "c"))

blacklisted_variants <- read_tsv(blacklisted_variants_file, col_types = cols(.default = "c"))

blacklisted_variants <- mutate(blacklisted_variants, Blacklist = TRUE)

variants <- variants %>%
  left_join(blacklisted_variants, by = c("Chromosome", "Position", "Ref", "Alt")) %>%
  mutate(Blacklist = replace_na(Blacklist, FALSE)) %>%
  mutate(Filters = ifelse(Blacklist, str_c(Filters, "blacklist", sep = ","), Filters)) %>%
  mutate(Filters = str_remove(Filters, "^PASS,")) %>%
  select(!Blacklist)

write_tsv(variants, output_file, na = "")

