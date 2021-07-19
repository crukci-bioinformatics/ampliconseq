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
minimum_proportion_of_variants <- 0.5
minimum_proportion_of_libraries <- 0.75
minimum_replicate_correlation <- 0.95

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

ids <- distinct(allele_fractions, ID)
variants <- distinct(allele_fractions, Amplicon, Chromosome, Position, Ref, Alt)
message("Libraries: ", nrow(ids))
message("Variants: ", nrow(variants))
message("Allele fractions: ", nrow(allele_fractions))

# plot distribution of allele fractions
# allele_fractions %>%
#   filter(`Allele fraction` >= minimum_variant_allele_fraction) %>%
#   ggplot(aes(x = `Allele fraction`)) +
#   geom_density()

# remove variants with narrow range of allele fractions across all samples
message("Excluding variants with narrow range of allele fractions across all samples")

allele_fractions <- allele_fractions %>%
  group_by(Amplicon, Chromosome, Position, Ref, Alt) %>%
  filter((max(`Allele fraction`) - min(`Allele fraction`)) >= minimum_allele_fraction_range) %>%
  ungroup()

ids <- distinct(allele_fractions, ID)
variants <- distinct(allele_fractions, Amplicon, Chromosome, Position, Ref, Alt)
message("Libraries: ", nrow(ids))
message("Variants: ", nrow(variants))

# exclude libraries with too many missing allele fractions
message("Excluding libraries with too many missing variant allele fractions")

allele_fractions <- allele_fractions %>%
  add_count(ID) %>%
  filter(n >= minimum_proportion_of_variants * nrow(variants)) %>%
  select(!n)

ids <- distinct(allele_fractions, ID)
variants <- distinct(allele_fractions, Amplicon, Chromosome, Position, Ref, Alt)
message("Libraries: ", nrow(ids))
message("Variants: ", nrow(variants))

# exclude variants with too many missing allele fractions
message("Excluding variants with too many missing samples/allele fractions")

allele_fractions <- allele_fractions %>%
  add_count(Amplicon, Chromosome, Position, Ref, Alt) %>%
  filter(n >= minimum_proportion_of_libraries * nrow(ids)) %>%
  select(!n)

ids <- distinct(allele_fractions, ID)
variants <- distinct(allele_fractions, Amplicon, Chromosome, Position, Ref, Alt)
message("Libraries: ", nrow(ids))
message("Variants: ", nrow(variants))

# allele fraction heatmap
message("Creating variant allele fraction heatmap")

matrix <- allele_fractions %>%
  transmute(ID, Sample, DisplayID = str_c(ID, Sample, sep = "  "), SNV = str_c(Amplicon, " ", Chromosome, ":", Position, " ", Ref, ">", Alt), `Allele fraction`) %>%
  pivot_wider(id_cols = c(ID, Sample, DisplayID), names_from = SNV, values_from = `Allele fraction`)

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

label_size <- max(8 - floor(ncol(matrix) / 25), 2)

heatmap <- Heatmap(
  matrix,
  name = "Allele fraction",
  col = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
  column_names_gp = gpar(fontsize = label_size),
  column_dend_height = unit(30, "mm"),
  show_row_names = FALSE,
  show_row_dend = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 8),
    grid_width = unit(4, "mm"),
    grid_height = unit(4, "mm")
  ),
  top_annotation = heatmapAnnotation
)

heatmap_width <- 12
heatmap_height <- 6

pdf(str_c(vaf_heatmap_prefix, ".pdf"), width = heatmap_width, height = heatmap_height)
draw(heatmap, show_annotation_legend = FALSE)
dev.off()

svglite(str_c(vaf_heatmap_prefix, ".svg"), width = heatmap_width, height = heatmap_height)
draw(heatmap, show_annotation_legend = FALSE)
dev.off()

rsvg_png(str_c(vaf_heatmap_prefix, ".svg"), str_c(vaf_heatmap_prefix, ".png"), width = heatmap_width * 1000, height = heatmap_height * 1000)

# correlation heatmap
message("Creating correlation heatmap")

hc <- hclust(dist(scale(t(matrix))))
# plot(hc)

correlation_matrix <- cor(matrix, use = "pairwise.complete.obs")

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

correlation_matrix <- correlation_matrix[hc$order, hc$order]

groups <- display_ids[hc$order,]$Sample
group_colours <- hue_pal()(length(groups))
names(group_colours) <- groups

heatmapAnnotation <- HeatmapAnnotation(
  df = data.frame(Group = groups),
  col = list(Group = group_colours),
  simple_anno_size = unit(2, "mm"),
  show_annotation_name = FALSE
)

label_size <- max(8 - floor(ncol(matrix) / 25), 2)

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

correlation_heatmap_width <- 8
correlation_heatmap_height <- 10 

pdf(str_c(vaf_correlation_heatmap_prefix, ".pdf"), width = correlation_heatmap_width, height = correlation_heatmap_height)
draw(heatmap, show_annotation_legend = FALSE)
dev.off()

svglite(str_c(vaf_correlation_heatmap_prefix, ".svg"), width = correlation_heatmap_width, height = correlation_heatmap_height)
draw(heatmap, show_annotation_legend = FALSE)
dev.off()

rsvg_png(str_c(vaf_correlation_heatmap_prefix, ".svg"), str_c(vaf_correlation_heatmap_prefix, ".png"), width = correlation_heatmap_width * 1000, height = correlation_heatmap_height * 1000)

# find pairs of sample replicates with low correlation but which have high
# correlation to a library from another sample
mismatched_replicates <- correlations %>%
  filter(Sample1 == Sample2) %>%
  filter(Correlation < minimum_replicate_correlation)

correlations %>%
  semi_join(mismatched_replicates, by = "Sample1") %>%
  filter(Sample1 == Sample2 | Correlation >= minimum_replicate_correlation) %>%
  select(Sample1, ID1, Sample2, ID2, Correlation) %>%
  arrange(Sample1, ID1, desc(Correlation)) %>%
  mutate(Correlation = round(Correlation, digits = 3)) %>%
  write_tsv(replicate_mismatch_file)

