#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
#
# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Makes minor modifications to the variant table output by GATK VariantsToTable
# including:
# - add ID column
# - separate allele depths into separate columns
# - rounds quality scores
# - computes the allele fraction from the allelic depths if no values extracted

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--input"), dest = "input_file",
              help = "Variant table from GATK VariantsToTable (AMPLICON, CHROM, POS, REF, ALT, MULTI-ALLELIC, TYPE, FILTER, QUAL, {ID}.DP {ID}.AD {ID}.AF columns required)"),

  make_option(c("--id"), dest = "id", help = "Library ID"),

  make_option(c("--output"), dest = "output_file",
              help = "Output file containing all variants including rows for missing calls")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

id <- opt$id
input_file <- opt$input_file
output_file <- opt$output_file

if (is.null(id)) stop("Library ID must be specified")
if (is.null(input_file)) stop("Input variants file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

suppressPackageStartupMessages(library(tidyverse))

# read variants file
variants <- read_tsv(input_file, col_types = cols(.default = "c"))

# rename columns
colnames(variants) <- c("Amplicon", "Chromosome", "Position", "Ref", "Alt", "Multiallelic", "Type", "Filters", "Quality", "Depth", "AD", "Allele fraction")

# separate Ref and Alt allelic depths from AD column
variants <- separate(variants, AD, into = c("Ref count", "Alt count"), sep = ",")

# round quality scores
variants <- mutate(variants, Quality = round(parse_double(Quality)))

# compute allele fraction if all values are missing, i.e. the variant caller
# doesn't create an AF field
if (all(is.na(variants$`Allele fraction`))) {
  variants <- variants %>%
    mutate(across(c(`Ref count`, `Alt count`, Depth), parse_integer)) %>%
    mutate(`Allele fraction` = round(`Alt count` / Depth, digits = 5))
}

# add ID column
variants <- variants %>%
  mutate(ID = id) %>%
  select(ID, everything())

# sort variants and write to output file
write_tsv(variants, output_file, na = "")

