args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
{
  stop("Usage: Rscript check_samples_file.R amplicon_intervals_file target_intervals_file reference_genome_index")
}

amplicon_intervals_file <- args[1]
target_intervals_file <- args[2]
reference_genome_index_file <- args[3]

suppressPackageStartupMessages(library(tidyverse))


# read and check amplicon intervals file
expected_columns <- c("ID", "Chromosome", "Start", "End")

if (str_ends(str_to_lower(amplicon_intervals_file), "\\.csv")) {
  amplicons <- read_csv(amplicon_intervals_file)

  missing_columns <- setdiff(expected_columns, colnames(amplicons))
  if (length(missing_columns) > 0) {
    stop("missing columns found in ", amplicon_intervals_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
  }
} else {
  amplicons <- read_tsv(amplicon_intervals_file, comment = "@", col_types = "cnncc", col_names = c("Chromosome", "Start", "End", "Strand", "ID")) %>%
    select(-Strand)
}

amplicons <- select(amplicons, all_of(expected_columns))

if (nrow(amplicons) == 0) {
  stop("no amplicon intervals found in: ", amplicon_intervals_file)
}

if (nrow(filter(amplicons, is.na(amplicons$ID))) > 0) {
  stop("missing IDs in ", amplicon_intervals_file)
}

duplicates <- amplicons %>%
  count(ID) %>%
  filter(n > 1)
if (nrow(duplicates) > 0) {
  stop("duplicate IDs found in ", amplicon_intervals_file, ": '", str_c(duplicates$ID, collapse = "', '"), "'")
}

incomplete <- filter(amplicons, if_any(everything(), ~ is.na(.x)))
if (nrow(incomplete) > 0) {
  stop("missing values found in ", amplicon_intervals_file)
}

incorrect_start_and_end <- filter(amplicons, Start < 1 | Start > End)
if (nrow(incorrect_start_and_end)) {
  stop("invalid start and/or end coordinates in ", amplicon_intervals_file, " for following amplicons: '", str_c(incorrect_start_and_end$ID, collapse = "', '"), "'")
}

duplicates <- amplicons %>%
  count(Chromosome, Start, End) %>%
  filter(n > 1) %>%
  left_join(amplicons, by = c("Chromosome", "Start", "End"))
if (nrow(duplicates) > 0) {
  stop("amplicons with same genomic coordinates in ", amplicon_intervals_file, ": '", str_c(duplicates$ID, collapse = "', '"), "'")
}


# read and check target intervals file
if (str_ends(str_to_lower(target_intervals_file), "\\.csv")) {
  targets <- read_csv(target_intervals_file)
  
  missing_columns <- setdiff(expected_columns, colnames(targets))
  if (length(missing_columns) > 0) {
    stop("missing columns found in ", target_intervals_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
  }
} else {
  targets <- read_tsv(target_intervals_file, comment = "@", col_types = "cnncc", col_names = c("Chromosome", "Start", "End", "Strand", "ID")) %>%
    select(-Strand)
}

targets <- select(targets, all_of(expected_columns))

if (nrow(targets) == 0) {
  stop("no target intervals found in: ", target_intervals_file)
}

if (nrow(filter(targets, is.na(targets$ID))) > 0) {
  stop("missing IDs in ", target_intervals_file)
}

incomplete <- filter(targets, if_any(everything(), ~ is.na(.x)))
if (nrow(incomplete) > 0) {
  stop("missing values found in ", target_intervals_file)
}

incorrect_start_and_end <- filter(targets, Start < 1 | Start > End)
if (nrow(incorrect_start_and_end)) {
  stop("invalid start and end coordinates in ", target_intervals_file, " for following targets: '", str_c(incorrect_start_and_end$ID, collapse = "', '"), "'")
}


# check for inconsistencies between target and amplicon intervals
missing <- anti_join(targets, amplicons, by = "ID")
if (nrow(missing) > 0) {
  stop("IDs found in ", target_intervals_file, " without corresponding entry in ", amplicon_intervals_file, ": '", str_c(missing$ID, collapse = "', '"), "'")
}

missing <- anti_join(amplicons, targets, by = "ID")
if (nrow(missing) > 0) {
  stop("IDs found in ", amplicon_intervals_file, " without corresponding entry in ", target_intervals_file, ": '", str_c(missing$ID, collapse = "', '"), "'")
}

targets <- rename(targets, TargetChromosome = Chromosome, TargetStart = Start, TargetEnd = End)
amplicons <- inner_join(amplicons, targets, by = "ID")

mismatched_chromosomes <- filter(amplicons, Chromosome != TargetChromosome)
if (nrow(mismatched_chromosomes) > 0) {
  stop("mismatched chromosomes in ", amplicon_intervals_file, " and ", target_intervals_file, " for the following amplicons: '", str_c(mismatched_chromosomes$ID, collapse = "', '"), "'")
}

targets_outside_amplicons <- filter(amplicons, TargetStart < Start | TargetEnd > End)
if (nrow(targets_outside_amplicons) > 0) {
  stop("Target intervals not within amplicon intervals: '", str_c(targets_outside_amplicons$ID, collapse = "', '"), "'")
}

amplicons <- select(amplicons, ID, Chromosome, Start, End, TargetStart, TargetEnd)


# read reference genome index file
chromosomes <- read_tsv(reference_genome_index_file, col_types = "cnnnn", col_names = c("Chromosome", "Length", "Offset", "Linebases", "Linewidth"))

# additional checks based on chromosomes from reference genome
missing_chromosomes <- anti_join(amplicons, chromosomes, by = "Chromosome")
if (nrow(missing_chromosomes) > 0) {
  stop("Amplicons found with chromosomes that don't match the reference genome: '", str_c(missing_chromosomes$ID, collapse = "', '"), "'")
}

invalid_coordinates <- amplicons %>%
  left_join(chromosomes, by = "Chromosome") %>%
  filter(End > Length)
if (nrow(invalid_coordinates) > 0) {
  stop("Amplicons found with out-of-bounds coordinates: '", str_c(invalid_coordinates$ID, collapse = "', '"), "'")
}


# sort into chromosome/position order
amplicons <- amplicons %>%
  mutate(Chromosome = factor(Chromosome, levels = chromosomes$Chromosome)) %>%
  arrange(Chromosome, Start, End)

message("Number of amplicons: ", nrow(amplicons))


# group amplicons into non-overlapping groups

find_overlaps <- function(query_id, query_chromosome, query_start, query_end) {
  amplicons %>%
    filter(Chromosome == query_chromosome, Start <= query_end, End >= query_start) %>%
    transmute(QueryID = query_id, OverlapID = ID) %>%
    select(ID = QueryID, OverlapID)
}

overlaps <- amplicons %>%
  select(query_id = ID, query_chromosome = Chromosome, query_start = Start, query_end = End) %>%
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

write_csv(amplicons, "amplicon_details.csv")


# write BED files for each group
for (i in 1:number_of_groups) {
  amplicon_group <- filter(amplicons, Group == i)

  amplicon_group %>%
    transmute(Chromosome, Start = Start - 1, End) %>%
    write_tsv(str_c("amplicons.", i, ".bed"), col_names = FALSE)

  amplicon_group %>%
    transmute(Chromosome, Start = TargetStart - 1, End = TargetEnd) %>%
    write_tsv(str_c("targets.", i, ".bed"), col_names = FALSE)
}

