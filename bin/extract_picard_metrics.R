#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core

# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Extract the metrics table from the Picard metrics file and adds an ID column
# containing the given identifier.

# Note that if an ID column was already present in the Picard metrics file it
# would be overwritten. Currently no Picard metrics files have an ID column.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--id"), dest = "id",
              help = "Identifier for the dataset used to populate an ID column in the output table"),

  make_option(c("--metrics"), dest = "metrics_file",
              help = "Picard metrics file"),

  make_option(c("--output"), dest = "output_file",
              help = "Output table")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

id <- opt$id
metrics_file <- opt$metrics_file
output_file <- opt$output_file

if (is.null(metrics_file)) stop("Metrics file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

suppressPackageStartupMessages(library(tidyverse))

lines <- readLines(metrics_file)
skip <- which(str_detect(lines, "^## METRICS CLASS"))
n_max <- Inf
empty_lines <- which(str_detect(lines[(skip + 1):length(lines)], "^$"))
if (!is_empty(empty_lines)) n_max <- empty_lines[1] - 2

metrics <- read_tsv(metrics_file, skip = skip, n_max = n_max, col_types = cols(.default = col_character()))

metrics %>%
  mutate(ID = id) %>%
  select(ID, everything()) %>%
  write_tsv(output_file, na = "")

