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

  make_option(c("--blacklist"), dest = "blacklist_file",
              help = "TSV file containing blacklisted variants (Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--vep-annotations"), dest = "vep_file",
              help = "TSV file containing Ensembl VEP annotations (optional; Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--offset-from-primer-end-annotations"), dest = "offset_from_primer_end_annotations_file",
              help = "TSV file containing annotations for offset from primer end (Amplicon, Chromosome, Position, Ref, Alt and Offset from primer end columns required)"),

  make_option(c("--other-annotations"), dest = "annotation_file",
              help = "TSV file containing additional annotations (Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--reference-sequence-index"), dest = "reference_sequence_index_file",
              help = "Index file for the reference genome sequence used for chromosome sort order (expected to have .fai extension)"),

  make_option(c("--output-prefix"), dest = "output_prefix",
              help = "Prefix for output variant summary files in CSV and TSV format"),

  make_option(c("--minimum-depth"), dest = "minimum_depth", type = "integer", default = 100,
              help = "Minimum depth for high-confidence variant calls (default: %default)")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

variants_file <- opt$variants_file
blacklist_file <- opt$blacklist_file
vep_file <- opt$vep_file
offset_from_primer_end_annotations_file <- opt$offset_from_primer_end_annotations_file
annotation_file <- opt$annotation_file
reference_sequence_index_file <- opt$reference_sequence_index_file
output_prefix <- opt$output_prefix
minimum_depth <- opt$minimum_depth

if (is.null(variants_file)) stop("Input variant file must be specified")
if (is.null(blacklist_file)) stop("Blacklisted variants file must be specified")
#if (is.null(vep_file)) stop("Ensembl VEP annotations file must be specified")
if (is.null(offset_from_primer_end_annotations_file)) stop("Offset from primer end annotations file must be specified")
if (is.null(annotation_file)) stop("Additional annotations file must be specified")
if (is.null(reference_sequence_index_file)) stop("Reference sequence index file must be specified")
if (is.null(output_prefix)) stop("Prefix for output files must be specified")

if (!is.numeric(minimum_depth) || minimum_depth <= 0) {
  stop("Invalid minimum depth of coverage for positions to be included in fitting noise distribution")
}

suppressPackageStartupMessages(library(tidyverse))

# read variants
variants <- read_tsv(variants_file, col_types = cols(
  Multiallelic = "l",
  `Allele fraction (pileup)` = "d",
  `Position noise threshold` = "d",
  `Library noise threshold` = "d",
  .default = "c"))

# rounding for allele fraction and noise thresholds
variants <- mutate(variants, across(c(`Allele fraction (pileup)`, `Position noise threshold`, `Library noise threshold`), ~ round(.x, digits = 5)))

# assign confidence based on the number of replicates in which the variant is
# called and passes filters and for which there was sufficient depth of coverage
confidence <- variants %>%
  mutate(Confident = Filters == "PASS" & Depth >= minimum_depth) %>%
  group_by(Sample, Amplicon, Chromosome, Position, Ref, Alt) %>%
  summarize(ConfidentCount = sum(Confident), ReplicateCount = n(), .groups = "drop") %>%
  mutate(Confidence = ifelse(ConfidentCount == 0, "low", "medium")) %>%
  mutate(Confidence = ifelse(ConfidentCount == ReplicateCount, "high", Confidence)) %>%
  select(!c(ConfidentCount, ReplicateCount))

# collapse/condense variants for all replicate libraries for a sample into a
# single row for each distinct variant
replicates <- variants %>%
  distinct(Sample, ID) %>%
  group_by(Sample) %>%
  mutate(ReplicateNumber = row_number()) %>%
  ungroup()

variants <- variants %>%
  select(
    Sample, Amplicon, Chromosome, Position, Ref, Alt, Specific,
    ID, Filters, Quality, Depth, `Allele fraction`,
    `Depth (pileup)`, `Allele fraction (pileup)`,
    `Position noise threshold`, `Library noise threshold`
  ) %>%
  left_join(confidence, by = c("Sample", "Amplicon", "Chromosome", "Position", "Ref", "Alt")) %>%
  left_join(replicates, by = c("Sample", "ID")) %>%
  pivot_wider(names_from = ReplicateNumber, values_from = !c(Sample:Specific, `Position noise threshold`, Confidence, ReplicateNumber), names_sep = " ") %>%
  select(
    Sample, Amplicon, Chromosome, Position, Ref, Alt, Specific, Confidence,
    starts_with("ID "), starts_with("Filters "), starts_with("Quality "),
    matches("^Depth [0-9]+$"), matches("^Allele fraction [0-9]+$"),
    starts_with("Depth (pileup) "), starts_with("Allele fraction (pileup) "),
    `Position noise threshold`, starts_with("Library noise threshold ")
  )

# read blacklist and add column to indicate which variants are blacklisted
blacklisted_variants <- read_tsv(blacklist_file, col_types = cols(.default = "c"))
blacklisted_variants <- mutate(blacklisted_variants, Blacklist = "true")
variants <- variants %>%
  left_join(blacklisted_variants, by = c("Chromosome", "Position", "Ref", "Alt")) %>%
  mutate(Blacklist = replace_na(Blacklist, "false")) %>%
  select(Sample:Specific, Blacklist, everything())

# read Ensembl VEP annotations and add to the variant table
if (!is.null(vep_file)) {
  vep_annotations <- read_tsv(vep_file, col_types = cols(.default = "c"))
  variants <- left_join(variants, vep_annotations, by = c("Chromosome", "Position", "Ref", "Alt"))
}

# read offset from primer end annotations and add to the variant table
offset_annotations <- read_tsv(offset_from_primer_end_annotations_file, col_types = cols(.default = "c"))
variants <- left_join(variants, offset_annotations, by = c("Amplicon", "Chromosome", "Position", "Ref", "Alt"))

# read additional annotations and add to the varaint table
annotations <- read_tsv(annotation_file, col_types = cols(.default = "c"))
variants <- left_join(variants, annotations, by = c("Chromosome", "Position", "Ref", "Alt"))

# read reference genome index file and use chromosome order in sorting
chromosomes <- read_tsv(reference_sequence_index_file, col_types = "cnnnn", col_names = c("Chromosome", "Length", "Offset", "Linebases", "Linewidth"))
variants <- variants %>%
  mutate(Chromosome = factor(Chromosome, levels = chromosomes$Chromosome)) %>%
  arrange(Sample, Chromosome, Position, Ref, Alt, Amplicon)

# write variant summary table to CSV and TSV files
write_csv(variants, str_c(output_prefix, ".csv"), na = "")
write_tsv(variants, str_c(output_prefix, ".txt"), na = "")

