#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Gather variant calls/details for replicate libraries into a single row and add
# annotations.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--variants"), dest = "variants_file",
              help = "TSV file containing variants (Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--vep-annotations"), dest = "vep_file",
              help = "TSV file containing Ensembl VEP annotations (Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--other-annotations"), dest = "annotation_file",
              help = "TSV file containing additional annotations (Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--output-prefix"), dest = "output_prefix",
              help = "Prefix for output variant summary files in CSV and TSV format")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

variants_file <- opt$variants_file
vep_file <- opt$vep_file
annotation_file <- opt$annotation_file
output_prefix <- opt$output_prefix

if (is.null(variants_file)) stop("Input variant file must be specified")
if (is.null(vep_file)) stop("Ensembl VEP annotations file must be specified")
if (is.null(annotation_file)) stop("Additional annotations file must be specified")
if (is.null(output_prefix)) stop("Prefix for output files must be specified")

suppressPackageStartupMessages(library(tidyverse))

variants <- read_tsv(variants_file, col_types = cols(.default = "c"))

# rounding for allele fraction and noise thresholds
variants <- mutate(variants, across(c(`Allele fraction (pileup)`, `Position noise threshold`, `Library noise threshold`), ~ round(parse_double(.x), digits = 5)))

# collapse/condense variants for all replicate libraries for a sample into a
# single row for each distinct variant
replicates <- variants %>%
  distinct(Sample, ID) %>%
  group_by(Sample) %>%
  mutate(ReplicateNumber = row_number()) %>%
  ungroup()

variants <- variants %>%
  select(
    Sample, ID, Amplicon, Chromosome, Position, Ref, Alt, Specific,
    # Multiallelic,
    Filters, Quality, Depth, `Allele fraction`,
    `Depth (pileup)`, `Allele fraction (pileup)`,
    `Position noise threshold`, `Library noise threshold`
  ) %>%
  left_join(replicates, by = c("Sample", "ID")) %>%
  pivot_wider(names_from = ReplicateNumber, values_from = !c(Sample, Amplicon, Chromosome, Position, Ref, Alt, Specific, `Position noise threshold`, ReplicateNumber), names_sep = " ") %>%
  select(
    Sample, Amplicon, Chromosome, Position, Ref, Alt, Specific,
    starts_with("ID "),
    # starts_with("Multiallelic "),
    starts_with("Filters "), starts_with("Quality "),
    matches("^Depth [0-9]+$"), matches("^Allele fraction [0-9]+$"),
    starts_with("Depth (pileup) "), starts_with("Allele fraction (pileup) "),
    `Position noise threshold`, starts_with("Library noise threshold ")
  )

# TODO add IGV links

# read Ensembl VEP annotations and add to the variant table
vep_annotations <- read_tsv(vep_file, col_types = cols(.default = "c"))

variants <- left_join(variants, vep_annotations, by = c("Chromosome", "Position", "Ref", "Alt"))


# read additional annotations and add to the varaint table
annotations <- read_tsv(annotation_file, col_types = cols(.default = "c"))

variants <- left_join(variants, annotations, by = c("Chromosome", "Position", "Ref", "Alt"))

# write variant summary table to CSV and TSV files
write_csv(variants, str_c(output_prefix, ".csv"), na = "")
write_tsv(variants, str_c(output_prefix, ".txt"), na = "")

