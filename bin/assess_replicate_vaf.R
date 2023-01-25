#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
#
# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Assess correlation of SNV allele fractions for sample replicates by creating
# heatmaps and attempting to identify mismatched replicates.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--samples"), dest = "samples_file",
              help = "Sample sheet file (ID and Sample columns required)"),

  make_option(c("--pileup-counts"), dest = "pileup_counts_file",
              help = "Pileup counts file"),

  make_option(c("--output-prefix"), dest = "output_prefix", default = "",
              help = "Prefix for output files including the allele fraction table, heatmap plots and mismatched replicate table")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

samples_file <- opt$samples_file
pileup_counts_file <- opt$pileup_counts_file
output_prefix <- opt$output_prefix

if (is.null(samples_file)) stop("Samples file must be specified")
if (is.null(pileup_counts_file)) stop("Pileup counts file must be specified")

suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(scales)
  library(svglite)
  library(rsvg)
})

allele_fraction_file <- str_c(output_prefix, "allele_fractions.txt")
vaf_heatmap_prefix <- str_c(output_prefix, "vaf_heatmap")
vaf_correlation_heatmap_prefix <- str_c(output_prefix, "vaf_correlation_heatmap")
replicate_mismatch_file <- str_c(output_prefix, "vaf_mismatched_replicates.txt")

minimum_depth <- 100
minimum_variant_allele_fraction <- 0.1
minimum_allele_fraction_range <- 0.25
minimum_proportion_of_variants <- 0.75
minimum_proportion_of_libraries <- 0.9
number_of_variants_to_sample <- 1000
maximum_number_of_variants <- 100
minimum_correlation_difference <- 0.05
minimum_mismatch_correlation <- 0.9

samples <- read_tsv(samples_file, col_types = cols(.default = "c"))
samples <- select(samples, ID, Sample)

ids <- NULL
amplicon_loci <- NULL
variants <- NULL
allele_fractions <- NULL

pileup_col_types <- cols(Position = "i", Depth = "i", `A count` = "i", `C count` = "i", `G count` = "i", `T count` = "i", .default = "c")

# function used in chunked reading of pileup file to find amplicon
# loci where there is a variant in at least one library
collect_variants <- function(data, pos) {
  chunk_ids <- distinct(data, ID)

  ids <<- ids %>%
    bind_rows(chunk_ids) %>%
    distinct()

  chunk_amplicon_loci <- distinct(data, Amplicon, Chromosome, Position)
  amplicon_loci <<- amplicon_loci %>%
    bind_rows(chunk_amplicon_loci) %>%
    distinct()

  chunk_variants <- data %>%
    filter(Depth >= minimum_depth) %>%
    select(ID, Amplicon, Chromosome, Position, Ref = `Reference base`, `A count`:`T count`, Depth) %>%
    pivot_longer(`A count`:`T count`, names_to = "Alt", values_to = "Count") %>%
    mutate(Alt = str_remove(Alt, " count$")) %>%
    filter(Ref != Alt) %>%
    mutate(`Allele fraction` = Count / Depth) %>%
    filter(`Allele fraction` >= minimum_variant_allele_fraction) %>%
    distinct(Amplicon, Chromosome, Position, Ref, Alt)

  variants <<- variants %>%
    bind_rows(chunk_variants) %>%
    distinct()
}

message("Reading pileup file and identifying variant loci")
# time_summary <- system.time(
result <- read_tsv_chunked(pileup_counts_file, SideEffectChunkCallback$new(collect_variants), chunk_size = 100000, col_types = pileup_col_types, progress = TRUE)
# )
# message("User time:    ", round(time_summary[["user.self"]]), "s")
# message("Elapsed time: ", round(time_summary[["elapsed"]]), "s")

message("Libraries: ", nrow(ids))
message("Amplicon loci: ", nrow(amplicon_loci))
message("Variants: ", nrow(variants))

# function used in second pass of chunked reading of the pileup file to extract
# the allele fractions for all libraries for those variants identified in the
# first pass
collect_allele_fractions_for_variants <- function(data, pos) {

  chunk_allele_fractions <- data %>%
    filter(Depth >= minimum_depth) %>%
    select(ID, Amplicon, Chromosome, Position, Ref = `Reference base`, `A count`:`T count`, Depth) %>%
    pivot_longer(`A count`:`T count`, names_to = "Alt", values_to = "Count") %>%
    mutate(Alt = str_remove(Alt, " count$")) %>%
    semi_join(variants, by = c("Amplicon", "Chromosome", "Position", "Ref", "Alt")) %>%
    mutate(`Allele fraction` = Count / Depth)

  allele_fractions <<- bind_rows(allele_fractions, chunk_allele_fractions)
}

message("Reading pileup file and obtaining variant allele fractions for all libraries")
# time_summary <- system.time(
result <- read_tsv_chunked(pileup_counts_file, SideEffectChunkCallback$new(collect_allele_fractions_for_variants), chunk_size = 100000, col_types = pileup_col_types, progress = TRUE)
# )
# message("User time:    ", round(time_summary[["user.self"]]), "s")
# message("Elapsed time: ", round(time_summary[["elapsed"]]), "s")

allele_fractions <- allele_fractions %>%
  left_join(samples, by = "ID") %>%
  select(ID, Sample, everything()) %>%
  arrange(ID, Chromosome, Position, Ref, Alt)

missing_sample_names <- allele_fractions %>%
  filter(is.na(Sample)) %>%
  distinct(ID)
if (nrow(missing_sample_names) > 0) {
  stop("missing sample names for ", str_c(missing_sample_names$ID, collapse = ", "))
}

allele_fractions %>%
  mutate(`Allele fraction` = sprintf("%.5f", `Allele fraction`)) %>%
  write_tsv(allele_fraction_file)

allele_fractions <- allele_fractions %>%
  mutate(Variant = str_c(Amplicon, " ", Chromosome, ":", Position, " ", Ref, ">", Alt), `Allele fraction`) %>%
  select(ID, Sample, Variant, `Allele fraction`)

message("Libraries: ", nrow(distinct(allele_fractions, ID)))
message("Variants: ", nrow(distinct(allele_fractions, Variant)))
message("Allele fractions: ", nrow(allele_fractions))

# plot distribution of allele fractions
# allele_fractions %>%
#   filter(`Allele fraction` >= minimum_variant_allele_fraction) %>%
#   ggplot(aes(x = `Allele fraction`)) +
#   geom_density()

# remove variants with narrow range of allele fractions across all samples
message("Excluding variants with narrow range of allele fractions across all samples")

allele_fractions <- allele_fractions %>%
  group_by(Variant) %>%
  filter((max(`Allele fraction`) - min(`Allele fraction`)) >= minimum_allele_fraction_range) %>%
  ungroup()

message("Libraries: ", nrow(distinct(allele_fractions, ID)))
message("Variants: ", nrow(distinct(allele_fractions, Variant)))

# exclude variants with too many missing allele fractions
message("Excluding variants with too many missing samples/allele fractions")

allele_fractions <- allele_fractions %>%
  add_count(Variant) %>%
  filter(n >= minimum_proportion_of_libraries * nrow(distinct(allele_fractions, ID))) %>%
  select(!n)

message("Libraries: ", nrow(distinct(allele_fractions, ID)))
message("Variants: ", nrow(distinct(allele_fractions, Variant)))

# sample to reduce the number of variants
variants <- distinct(allele_fractions, Variant)
if (nrow(variants) > number_of_variants_to_sample) {
  message("Sampling random subset of variants")
  variants <- sample_n(variants, number_of_variants_to_sample)
  allele_fractions <- semi_join(allele_fractions, variants, by = "Variant")
  message("Libraries: ", nrow(distinct(allele_fractions, ID)))
  message("Variants: ", nrow(distinct(allele_fractions, Variant)))
}

# remove highly correlated variants
variants <- distinct(allele_fractions, Variant)
variant_count <- nrow(variants)

if (variant_count > maximum_number_of_variants) {
  message("Removing highly correlated variants")

  correlations <- allele_fractions %>%
    select(ID, Variant, `Allele fraction`) %>%
    pivot_wider(id_cols = ID, names_from = Variant, values_from = `Allele fraction`) %>%
    column_to_rownames(var = "ID") %>%
    as.matrix() %>%
    cor(use = "pairwise.complete.obs") %>%
    as.data.frame() %>%
    rownames_to_column(var = "Variant1") %>%
    as_tibble() %>%
    pivot_longer(!Variant1, names_to = "Variant2", values_to = "Correlation") %>%
    filter(Variant1 < Variant2) %>%
    arrange(desc(Correlation))

  while (variant_count > maximum_number_of_variants) {
    variant_to_remove <- correlations$Variant2[1]
    correlations <- filter(correlations, Variant1 != variant_to_remove, Variant2 != variant_to_remove)
    variant_count <- variant_count - 1
  }

  variants <- unique(c(correlations$Variant1, correlations$Variant2))
  allele_fractions <- filter(allele_fractions, Variant %in% variants)

  message("Libraries: ", nrow(distinct(allele_fractions, ID)))
  message("Variants: ", nrow(distinct(allele_fractions, Variant)))
}

# exclude libraries with too many missing allele fractions
message("Excluding libraries with too many missing variant allele fractions")

allele_fractions <- allele_fractions %>%
  add_count(ID) %>%
  filter(n >= minimum_proportion_of_variants * nrow(distinct(allele_fractions, Variant))) %>%
  select(!n)

message("Libraries: ", nrow(distinct(allele_fractions, ID)))
message("Variants: ", nrow(distinct(allele_fractions, Variant)))

# allele fraction heatmap
message("Creating variant allele fraction heatmap")

heatmap_width <- 12
heatmap_height <- 6

if (nrow(allele_fractions) > 0) {

  matrix <- allele_fractions %>%
    transmute(ID, Sample, DisplayID = str_c(ID, Sample, sep = "  "), Variant, `Allele fraction`) %>%
    pivot_wider(id_cols = c(ID, Sample, DisplayID), names_from = Variant, values_from = `Allele fraction`)

  display_ids <- select(matrix, ID, Sample, DisplayID)

  groups <- display_ids$Sample
  group_colours <- hue_pal()(length(groups))
  names(group_colours) <- groups

  heatmapAnnotation <- HeatmapAnnotation(
    df = data.frame(Group = groups),
    col = list(Group = group_colours),
    simple_anno_size = unit(2, "mm"),
    show_annotation_name = FALSE
  )

  matrix <- matrix %>%
    select(!c(ID, Sample)) %>%
    column_to_rownames(var = "DisplayID") %>%
    as.matrix() %>%
    t()

  row_label_size <- max(8 - floor(nrow(matrix) / 10), 1.5)
  column_label_size <- max(8 - floor(ncol(matrix) / 25), 1.5)

  heatmap <- Heatmap(
    matrix,
    name = "Allele fraction",
    col = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
    row_names_side = "left",
    row_names_gp = gpar(fontsize = row_label_size),
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = column_label_size),
    column_dend_height = unit(30, "mm"),
    show_row_dend = FALSE,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 9),
      labels_gp = gpar(fontsize = 8),
      grid_width = unit(4, "mm"),
      grid_height = unit(4, "mm")
    ),
    top_annotation = heatmapAnnotation
  )

  pdf(str_c(vaf_heatmap_prefix, ".pdf"), width = heatmap_width, height = heatmap_height)
  draw(heatmap, show_annotation_legend = FALSE)
  dev.off()

  svglite(str_c(vaf_heatmap_prefix, ".svg"), width = heatmap_width, height = heatmap_height)
  draw(heatmap, show_annotation_legend = FALSE)
  dev.off()

  rsvg_png(str_c(vaf_heatmap_prefix, ".svg"), str_c(vaf_heatmap_prefix, ".png"), width = heatmap_width * 1000, height = heatmap_height * 1000)

} else {

  empty_plot <- ggplot() + labs(title = "Insufficient data to create VAF heatmap - too few variant loci to cluster samples/libraries")

  pdf(str_c(vaf_heatmap_prefix, ".pdf"), width = heatmap_width, height = heatmap_height)
  print(empty_plot)
  dev.off()

  svglite(str_c(vaf_heatmap_prefix, ".svg"), width = heatmap_width, height = heatmap_height)
  print(empty_plot)
  dev.off()

  rsvg_png(str_c(vaf_heatmap_prefix, ".svg"), str_c(vaf_heatmap_prefix, ".png"), width = heatmap_width * 1000, height = heatmap_height * 1000)
}

# correlation heatmap
message("Creating correlation heatmap")

correlation_heatmap_width <- 8
correlation_heatmap_height <- 10

if (nrow(allele_fractions) > 0 && nrow(matrix) > 1 && ncol(matrix) > 1) {

  correlation_matrix <- cor(matrix, use = "pairwise.complete.obs")

  label_size <- max(8 - floor(ncol(correlation_matrix) / 25), 1.25)

  heatmap <- Heatmap(
    correlation_matrix,
    name = "r",
    col = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
    row_names_side = "left",
    row_names_gp = gpar(fontsize = label_size),
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = label_size),
    column_dend_height = unit(30, "mm"),
    show_row_dend = FALSE,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 9),
      labels_gp = gpar(fontsize = 8),
      grid_width = unit(4, "mm"),
      grid_height = unit(4, "mm")
    ),
    top_annotation = heatmapAnnotation
  )

  pdf(str_c(vaf_correlation_heatmap_prefix, ".pdf"), width = correlation_heatmap_width, height = correlation_heatmap_height)
  draw(heatmap, show_annotation_legend = FALSE)
  dev.off()

  svglite(str_c(vaf_correlation_heatmap_prefix, ".svg"), width = correlation_heatmap_width, height = correlation_heatmap_height)
  draw(heatmap, show_annotation_legend = FALSE)
  dev.off()

  rsvg_png(str_c(vaf_correlation_heatmap_prefix, ".svg"), str_c(vaf_correlation_heatmap_prefix, ".png"), width = correlation_heatmap_width * 1000, height = correlation_heatmap_height * 1000)

  # find pairs of sample replicates with low correlation but which have high
  # correlation to a library from another sample
  correlations <- correlation_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = "DisplayID1") %>%
    pivot_longer(!DisplayID1, names_to = "DisplayID2", values_to = "Correlation") %>%
    left_join(display_ids, by = c("DisplayID1" = "DisplayID")) %>%
    rename(ID1 = ID, Sample1 = Sample) %>%
    left_join(display_ids, by = c("DisplayID2" = "DisplayID")) %>%
    rename(ID2 = ID, Sample2 = Sample) %>%
    select(ID1, Sample1, ID2, Sample2, Correlation) %>%
    filter(ID1 != ID2)

  # find replicates more highly correlated with libraries from another sample
  mismatched_replicates <- correlations %>%
    filter(Sample1 == Sample2) %>%
    filter(Correlation < (1 - minimum_correlation_difference)) %>%
    select(Sample = Sample1, ID = ID1, `Replicate ID` = ID2, Correlation) %>%
    left_join(select(correlations, ID = ID1, `Sample 2` = Sample2, `ID 2` = ID2, `Correlation 2` = Correlation), by = "ID") %>%
    filter(Sample != `Sample 2`) %>%
    filter(`Correlation 2` >= minimum_mismatch_correlation) %>%
    filter((`Correlation 2` - Correlation) >= minimum_correlation_difference) %>%
    group_by(ID, `Replicate ID`) %>%
    slice_max(order_by = `Correlation 2`, n = 5) %>%
    ungroup()

  mismatched_replicates %>%
    arrange(Sample, ID, `Replicate ID`, `Sample 2`, `ID 2`) %>%
    mutate(across(c(Correlation, `Correlation 2`), round, digits = 3)) %>%
    write_tsv(replicate_mismatch_file)

} else {

  empty_plot <- ggplot() + labs(title = "Insufficient data to create correlation heatmap")

  pdf(str_c(vaf_correlation_heatmap_prefix, ".pdf"), width = correlation_heatmap_width, height = correlation_heatmap_height)
  print(empty_plot)
  dev.off()

  svglite(str_c(vaf_correlation_heatmap_prefix, ".svg"), width = correlation_heatmap_width, height = correlation_heatmap_height)
  print(empty_plot)
  dev.off()

  rsvg_png(str_c(vaf_correlation_heatmap_prefix, ".svg"), str_c(vaf_correlation_heatmap_prefix, ".png"), width = correlation_heatmap_width * 1000, height = correlation_heatmap_height * 1000)

  tibble(
    Sample = character(0),
    ID = character(0),
    `Replicate ID` = character(0),
    Correlation = numeric(0),
    `Sample 2` = character(0),
    `ID 2` = character(0),
    `Correlation 2` = numeric(0)
  ) %>%
    write_tsv(replicate_mismatch_file)
}
