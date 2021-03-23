#!/usr/bin/env Rscript

# Creates groups of non-overlapping groups of amplicons in which no amplicon
# within a group overlaps with another amplicon in the same group.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
{
  stop("Usage: Rscript create_non_overlapping_amplicon_groups.R amplicons_file reference_sequence_index")
}

amplicons_file <- args[1]
reference_sequence_index_file <- args[2]

suppressPackageStartupMessages(library(tidyverse))


# read and check amplicons file
if (str_ends(str_to_lower(amplicons_file), "\\.csv")) {
  amplicons <- read_csv(amplicons_file, col_types = cols(.default = col_character()))
} else {
  amplicons <- read_tsv(amplicons_file, col_types = cols(.default = col_character()))
}

expected_columns <- c("ID", "Chromosome", "AmpliconStart", "AmpliconEnd", "TargetStart", "TargetEnd", "Gene")
missing_columns <- setdiff(expected_columns, colnames(amplicons))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", amplicons_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}
amplicons <- select(amplicons, all_of(expected_columns))

if (nrow(amplicons) == 0) {
  stop("no amplicon intervals found in: ", amplicons_file)
}

if (nrow(filter(amplicons, is.na(amplicons$ID))) > 0) {
  stop("missing IDs in ", amplicons_file)
}

duplicates <- amplicons %>%
  count(ID) %>%
  filter(n > 1)
if (nrow(duplicates) > 0) {
  stop("duplicate IDs found in ", amplicons_file, ": '", str_c(duplicates$ID, collapse = "', '"), "'")
}

missing_values <- filter(amplicons, if_any(Chromosome:TargetEnd, ~ is.na(.x)))
if (nrow(missing_values) > 0) {
  stop("amplicons with missing values found in ", amplicons_file, ": '", str_c(missing_values$ID, collapse = "', '"), "'")
}

duplicates <- amplicons %>%
  count(Chromosome, AmpliconStart, AmpliconEnd) %>%
  filter(n > 1) %>%
  left_join(amplicons, by = c("Chromosome", "AmpliconStart", "AmpliconEnd"))
if (nrow(duplicates) > 0) {
  stop("amplicons with same genomic coordinates in ", amplicons_file, ": '", str_c(duplicates$ID, collapse = "', '"), "'")
}


# convert coordinates to integer values
amplicons <- mutate(amplicons, across(AmpliconStart:TargetEnd, parse_integer))
missing_values <- filter(amplicons, if_any(AmpliconStart:TargetEnd, ~ is.na(.x)))
if (nrow(missing_values) > 0) {
  stop("amplicons with non-integer coordinates found in ", amplicons_file, ": '", str_c(missing_values$ID, collapse = "', '"), "'")
}

incorrect_amplicon_coordinates <- filter(amplicons, AmpliconStart < 1 | AmpliconStart > AmpliconEnd)
if (nrow(incorrect_amplicon_coordinates)) {
  stop("invalid amplicon start and/or end coordinates in ", amplicons_file, " for following amplicons: '", str_c(incorrect_amplicon_coordinates$ID, collapse = "', '"), "'")
}

incorrect_target_coordinates <- filter(amplicons, TargetStart < AmpliconStart | TargetEnd > AmpliconEnd | TargetStart > TargetEnd)
if (nrow(incorrect_target_coordinates)) {
  stop("invalid target start and/or end coordinates in ", amplicons_file, " for following amplicons: '", str_c(incorrect_target_coordinates$ID, collapse = "', '"), "'")
}

# read reference genome index file
chromosomes <- read_tsv(reference_sequence_index_file, col_types = "cnnnn", col_names = c("Chromosome", "Length", "Offset", "Linebases", "Linewidth"))

# additional checks based on chromosomes from reference genome
missing_chromosomes <- anti_join(amplicons, chromosomes, by = "Chromosome")
if (nrow(missing_chromosomes) > 0) {
  stop("Amplicons found with chromosomes that don't match the reference genome: '", str_c(missing_chromosomes$ID, collapse = "', '"), "'")
}

invalid_coordinates <- amplicons %>%
  left_join(chromosomes, by = "Chromosome") %>%
  filter(AmpliconEnd > Length)
if (nrow(invalid_coordinates) > 0) {
  stop("Amplicons found with out-of-bounds coordinates: '", str_c(invalid_coordinates$ID, collapse = "', '"), "'")
}


# sort into chromosome/position order
amplicons <- amplicons %>%
  mutate(Chromosome = factor(Chromosome, levels = chromosomes$Chromosome)) %>%
  arrange(Chromosome, AmpliconStart, AmpliconEnd)

message("Number of amplicons: ", nrow(amplicons))


# group amplicons into non-overlapping groups

find_overlaps <- function(query_id, query_chromosome, query_start, query_end) {
  amplicons %>%
    filter(Chromosome == query_chromosome, AmpliconStart <= query_end, AmpliconEnd >= query_start) %>%
    transmute(QueryID = query_id, OverlapID = ID) %>%
    select(ID = QueryID, OverlapID)
}

overlaps <- amplicons %>%
  select(query_id = ID, query_chromosome = Chromosome, query_start = AmpliconStart, query_end = AmpliconEnd) %>%
  pmap_dfr(find_overlaps)

overlap_counts <- count(overlaps, ID)
max_overlap_count <- max(overlap_counts$n)
message("Maximum number of overlapping amplicons: ", max_overlap_count)

number_of_groups <- max_overlap_count

while (TRUE) {
  groups <- amplicons %>%
    select(ID) %>%
    mutate(Group = rep_len(1:number_of_groups, nrow(.)))

  # check whether any overlaps between amplicons assigned to different groups
  overlapping_groups <- overlaps %>%
    filter(ID != OverlapID) %>%
    left_join(groups, by = "ID") %>%
    left_join(select(groups, OverlapID = ID, OverlapGroup = Group), by = "OverlapID") %>%
    filter(Group == OverlapGroup)

  if (nrow(overlapping_groups) == 0) break

  number_of_groups <- number_of_groups + 1
}

message("Number of non-overlapping interval groups: ", number_of_groups)

amplicons <- left_join(amplicons, groups, by = "ID")

write_tsv(amplicons, "amplicon_groups.txt")


# write separate amplicon details files for each group
for (group in 1:number_of_groups) {
  amplicons %>%
    filter(Group == group) %>%
    write_tsv(str_c("amplicon_groups.", group, ".txt"))
}

