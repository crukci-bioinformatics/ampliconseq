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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
{
  message("Usage: extract_picard_metrics.R id metrics_file output_file")
  quit(status = 1)
}

id <- args[1]
metrics_file <- args[2]
output_file <- args[3]

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

