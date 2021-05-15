#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
#
# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Computes background noise levels for substitutions at each amplicon
# position fitting allele fractions for the substitution across all
# samples/libraries to a beta distribution and applying a probability-based
# upper threshold.

# Also computes the background noise level for each substitution type within
# a sample/library by fitting the allele fractions for all substitutions of that
# type across all corresponding amplicon positions within the sample/library.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--pileup-counts"), dest = "pileup_counts_file",
              help = "Pileup counts file with read counts for each allele at every amplicon position"),

  make_option(c("--variants"), dest = "variants_file",
              help = "File containing variants used to restrict positions at which to compute background noise (Amplicon, Chromosome, Position, Ref and Alt columns required)"),

  make_option(c("--position-thresholds"), dest = "position_thresholds_file",
              help = "Output thresholds file for background noise for each substitution at each amplicon position"),

  make_option(c("--library-thresholds"), dest = "library_thresholds_file",
              help = "Output thresholds file for background noise for each substitution type within each sample/library"),

  make_option(c("--minimum-depth"), dest = "minimum_depth", type = "integer", default = 100,
              help = "Minimum depth for amplicon loci to be included in the fitting of background noise distributions (default: %default)"),

  make_option(c("--exclude-highest-fraction"), dest = "exclude_highest_fraction", type = "double", default = 0.1,
              help = "Fraction of allele fraction observations with highest values to be excluded from fitting  (default: %default)"),

  make_option(c("--maximum-allele-fraction"), dest = "maximum_allele_fraction", type = "double", default = 0.03,
              help = "Maximum value above which to exclude allele fractions from fitting (default: %default)"),

  make_option(c("--minimum-number-for-fitting"), dest = "minimum_number_for_fitting", type = "integer", default = 10,
              help = "Minimum number of allele fractions required for fitting (default: %default)"),

  make_option(c("--chunk-size"), dest = "chunk_size", type = "integer", default = 500000,
              help = "Maximum number of pileup count rows to process in a chunk (default: %default)"),

  make_option(c("--read-chunk-size"), dest = "read_chunk_size", type = "integer", default = 100000,
              help = "Chunk size for reading pileup count records prior to chunking for processing (default: %default)")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

pileup_counts_file <- opt$pileup_counts_file
variants_file <- opt$variants_file
position_thresholds_file <- opt$position_thresholds_file
library_thresholds_file <- opt$library_thresholds_file
minimum_depth <- opt$minimum_depth
exclude_highest_fraction <- opt$exclude_highest_fraction
maximum_allele_fraction <- opt$maximum_allele_fraction
minimum_number_for_fitting <- opt$minimum_number_for_fitting
chunk_size <- opt$chunk_size
read_chunk_size <- opt$read_chunk_size

if (is.null(pileup_counts_file)) stop("Pileup counts file must be specified")
if (is.null(position_thresholds_file)) stop("Position noise thresholds output file must be specified")
if (is.null(library_thresholds_file)) stop("Library noise thresholds output file must be specified")

if (!is.integer(minimum_depth) || minimum_depth <= 0) {
  stop("Invalid minimum depth of coverage for positions to be included in fitting noise distribution")
}

if (!is.double(exclude_highest_fraction) || exclude_highest_fraction < 0 || exclude_highest_fraction >= 1) {
  stop("Invalid value given for highest fraction of allele fractions to exclude from fitting noise distribution")
}

if (!is.double(maximum_allele_fraction) || maximum_allele_fraction <= 0 || maximum_allele_fraction > 1) {
  stop("Invalid value given for maximum allele fraction to include in fitting noise distributions")
}

if (!is.integer(minimum_number_for_fitting) || minimum_number_for_fitting < 5) {
  stop("Invalid value given for minimum number of allele fractions required for fitting noise distributions")
}

if (!is.integer(chunk_size) || chunk_size < 100000) {
  stop("Invalid chunk size for fitting noise distributions")
}

if (!is.integer(read_chunk_size) || read_chunk_size < 100000) {
  stop("Invalid chunk size for reading pileup count data")
}

suppressPackageStartupMessages({
  library(fitdistrplus)
  library(tidyverse)
})


# function to compute a threshold allele fraction from a set of allele
# fractions and a threshold probability
compute_threshold <- function(
  allele_fractions,
  p = 0.9999,
  exclude_zero_values = TRUE,
  exclude_highest_fraction = 0.0,
  maximum_allele_fraction = 1.0,
  minimum_number_for_fitting = 10)
{
  if (is.null(allele_fractions) ||
      length(allele_fractions) < minimum_number_for_fitting ||
      exclude_highest_fraction >= 1) {
    return(NA)
  }

  filtered_allele_fractions <- sort(allele_fractions)

  if (exclude_highest_fraction > 0)
  {
    n <- length(allele_fractions)
    filtered_allele_fractions <- filtered_allele_fractions[1:(n - as.integer(n * exclude_highest_fraction))]
  }

  if (exclude_zero_values) {
    filtered_allele_fractions <- filtered_allele_fractions[filtered_allele_fractions != 0]
  }

  filtered_allele_fractions <- filtered_allele_fractions[filtered_allele_fractions <= maximum_allele_fraction]

  if (length(filtered_allele_fractions) < minimum_number_for_fitting) return(NA)

  fit <- suppressWarnings(fitdist(filtered_allele_fractions, "beta", method = "mge"))

  qbeta(p, shape1 = fit$estimate[1], shape2 = fit$estimate[2])
}


message(Sys.time(), "  Reading pileup counts file")

total <- 0

amplicon_position_row_counts <- tibble(
  Amplicon = character(),
  Chromosome = character(),
  Position = integer(),
  n = integer()
)

library_reference_base_row_counts <- tibble(
  ID = character(),
  Ref = character(),
  n = integer()
)

# function to count pileup records for each amplicon position and for each
# sample/library and reference base
collect_row_counts <- function(data, pos) {
  total <<- total + nrow(data)

  data <- rename(data, Ref = `Reference base`)

  amplicon_position_row_counts <<- data %>%
    count(Amplicon, Chromosome, Position) %>%
    bind_rows(amplicon_position_row_counts) %>%
    count(Amplicon, Chromosome, Position, wt = n)

  library_reference_base_row_counts <<- data %>%
    count(ID, Ref) %>%
    bind_rows(library_reference_base_row_counts) %>%
    count(ID, Ref, wt = n)
}

result <- read_tsv_chunked(pileup_counts_file, SideEffectChunkCallback$new(collect_row_counts), chunk_size = read_chunk_size, col_types = "cccdcddddddd")

library_row_counts <- count(library_reference_base_row_counts, ID, wt = n)

number_of_libraries <- nrow(library_row_counts)
number_of_positions <- nrow(amplicon_position_row_counts)

message("Total number of pileup count records: ", total)
message("Number of libraries: ", number_of_libraries)
message("Distinct target positions: ", number_of_positions)

compute_position_thresholds <- nrow(amplicon_position_row_counts) > 0 && max(amplicon_position_row_counts$n) >= minimum_number_for_fitting
compute_library_thresholds <- nrow(library_reference_base_row_counts) > 0 && max(library_reference_base_row_counts$n) >= minimum_number_for_fitting

chunk_file_prefix <- tempfile("pileup_counts.", getwd())


# Compute position substitution thresholds
# ----------------------------------------

create_position_chunk_files <- function(data, pos) {
  data <- rename(data, Ref = `Reference base`)
  for (chunk in 1:number_of_chunks) {
    chunk_file = str_c(chunk_file_prefix, ".", chunk)
    chunk_data <- filter(amplicon_position_row_counts, Chunk == chunk)
    data %>%
      semi_join(chunk_data, by = c("Amplicon", "Chromosome", "Position")) %>%
      write_tsv(chunk_file, append = pos > 1)
  }
}

if (compute_position_thresholds) {

  message(as.character(Sys.time()), "  Creating position chunk files")

  # restrict positions to those for the variants if provided
  variants <- NULL

  if (!is.null(variants_file)) {
    variants <- read_tsv(variants_file, col_types = cols(.default = "c"))

    # check for expected columns
    expected_columns <- c("Amplicon", "Chromosome", "Position", "Ref", "Alt")
    missing_columns <- setdiff(expected_columns, colnames(variants))
    if (length(missing_columns) > 0) {
      stop("missing columns found in ", variants_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
    }

    variants <- variants %>%
      filter(nchar(Ref) == 1, nchar(Alt) == 1) %>%
      select(all_of(expected_columns)) %>%
      distinct() %>%
      mutate(Position = parse_integer(Position))

    positions <- distinct(variants, Amplicon, Chromosome, Position)

    amplicon_position_row_counts <- semi_join(amplicon_position_row_counts, positions, by = c("Amplicon", "Chromosome", "Position"))

    number_of_positions <- nrow(amplicon_position_row_counts)
    message("Distinct target positions (variants only): ", number_of_positions)
  }

  number_of_chunks <- ceiling(sum(amplicon_position_row_counts$n) / chunk_size)
  number_of_rows_per_chunk <- ceiling(nrow(amplicon_position_row_counts) / number_of_chunks)
  amplicon_position_row_counts <- amplicon_position_row_counts %>%
    arrange(Amplicon, Chromosome, Position) %>%
    mutate(Chunk = rep(1:number_of_chunks, each = number_of_rows_per_chunk)[1:nrow(amplicon_position_row_counts)])
  number_of_chunks <- max(amplicon_position_row_counts$Chunk)
  message("Number of chunks: ", number_of_chunks)

  result <- read_tsv_chunked(pileup_counts_file, SideEffectChunkCallback$new(create_position_chunk_files), chunk_size = read_chunk_size, col_types = "cccdcddddddd")

  message(as.character(Sys.time()), "  Computing position substitution thresholds")

  for (chunk in 1:number_of_chunks) {

    message(as.character(Sys.time()), "  Chunk ", chunk, " of ", number_of_chunks)

    chunk_file <- str_c(chunk_file_prefix, ".", chunk)

    pileup_counts <- read_tsv(chunk_file, col_types = "cccdcddddddd")

    pileup_counts <- pileup_counts %>%
      filter(Depth >= minimum_depth) %>%
      select(ID, Amplicon, Chromosome, Position, Ref, Depth, A = `A count`, C = `C count`, G = `G count`, T = `T count`) %>%
      pivot_longer(c(A, C, G, T), names_to = "Alt", values_to = "Count") %>%
      filter(Ref != Alt)

    if (!is.null(variants)) {
      pileup_counts <- semi_join(pileup_counts, variants, by = c("Amplicon", "Chromosome", "Position", "Ref", "Alt"))
    }

    allele_fractions <- pileup_counts %>%
      mutate(`Allele fraction` = Count / Depth) %>%
      select(ID, Amplicon, Chromosome, Position, Ref, Alt, `Allele fraction`) %>%
      arrange(ID, Amplicon, Chromosome, Position, Ref, Alt)

    time_summary <- system.time(
      thresholds <- allele_fractions %>%
        group_by(Amplicon, Chromosome, Position, Ref, Alt) %>%
        summarize(
          `Allele fraction threshold` = compute_threshold(
            `Allele fraction`,
            exclude_highest_fraction = exclude_highest_fraction,
            maximum_allele_fraction = maximum_allele_fraction,
            minimum_number_for_fitting = minimum_number_for_fitting),
          .groups = "drop") %>%
        filter(!is.na(`Allele fraction threshold`)) %>%
        mutate(`Allele fraction threshold` = sprintf("%.6f", `Allele fraction threshold`))
    )

    message("User time:    ", time_summary[["user.self"]], "s")
    message("Elapsed time: ", time_summary[["elapsed"]], "s")

    write_tsv(thresholds, position_thresholds_file, append = chunk > 1)

    file.remove(chunk_file)
  }

} else {

  message("Insufficient pileup count data to model background noise for position substitutions")

  write_tsv(
    tibble(
      Amplicon = character(),
      Chromosome = character(),
      Position = integer(),
      Ref = character(),
      Alt = character(),
      `Allele fraction threshold` = double()
    ),
    position_thresholds_file
  )
}


# Compute library substitution thresholds
# ---------------------------------------

create_library_chunk_files <- function(data, pos) {
  data <- rename(data, Ref = `Reference base`)
  for (chunk in 1:number_of_chunks) {
    chunk_file = str_c(chunk_file_prefix, ".", chunk)
    chunk_data <- filter(library_row_counts, Chunk == chunk)
    data %>%
      semi_join(chunk_data, by = "ID") %>%
      write_tsv(chunk_file, append = pos > 1)
  }
}

if (compute_library_thresholds) {

  message(as.character(Sys.time()), "  Creating library chunk files")

  number_of_chunks <- ceiling(sum(library_row_counts$n) / chunk_size)
  number_of_rows_per_chunk <- ceiling(nrow(library_row_counts) / number_of_chunks)
  library_row_counts <- library_row_counts %>%
    arrange(ID) %>%
    mutate(Chunk = rep(1:number_of_chunks, each = number_of_rows_per_chunk)[1:nrow(library_row_counts)])
  number_of_chunks <- max(library_row_counts$Chunk)
  message("Number of chunks: ", number_of_chunks)

  result <- read_tsv_chunked(pileup_counts_file, SideEffectChunkCallback$new(create_library_chunk_files), chunk_size = read_chunk_size, col_types = "cccdcddddddd")

  message(as.character(Sys.time()), "  Computing library substitution thresholds")

  for (chunk in 1:number_of_chunks) {

    message(as.character(Sys.time()), "  Chunk ", chunk, " of ", number_of_chunks)

    chunk_file <- str_c(chunk_file_prefix, ".", chunk)

    pileup_counts <- read_tsv(chunk_file, col_types = "cccdcddddddd")

    allele_fractions <- pileup_counts %>%
      filter(Depth >= minimum_depth) %>%
      select(ID, Amplicon, Chromosome, Position, Ref, Depth, A = `A count`, C = `C count`, G = `G count`, T = `T count`) %>%
      pivot_longer(c(A, C, G, T), names_to = "Alt", values_to = "Count") %>%
      filter(Ref != Alt) %>%
      mutate(`Allele fraction` = Count / Depth) %>%
      select(ID, Amplicon, Chromosome, Position, Ref, Alt, `Allele fraction`) %>%
      arrange(ID, Amplicon, Chromosome, Position, Ref, Alt)

    time_summary <- system.time(
      thresholds <- allele_fractions %>%
        group_by(ID, Ref, Alt) %>%
        summarize(
          `Allele fraction threshold` = compute_threshold(
            `Allele fraction`,
            exclude_highest_fraction = exclude_highest_fraction,
            maximum_allele_fraction = maximum_allele_fraction,
            minimum_number_for_fitting = minimum_number_for_fitting),
          .groups = "drop") %>%
        filter(!is.na(`Allele fraction threshold`)) %>%
        mutate(`Allele fraction threshold` = sprintf("%.6f", `Allele fraction threshold`))
    )

    message("User time:    ", time_summary[["user.self"]], "s")
    message("Elapsed time: ", time_summary[["elapsed"]], "s")

    write_tsv(thresholds, library_thresholds_file, append = chunk > 1)

    file.remove(chunk_file)
  }

} else {

  message("Insufficient pileup count data to model background noise for library substitutions")

  write_tsv(
    tibble(
      ID = character(),
      Ref = character(),
      Alt = character(),
      `Allele fraction threshold` = double()
    ),
    library_thresholds_file
  )
}


message(as.character(Sys.time()), "  Finished")

