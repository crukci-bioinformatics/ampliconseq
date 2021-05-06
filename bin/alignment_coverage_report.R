#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
#
# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Create alignment summary report and merged alignment and coverage metrics
# table from the collated Picard alignment and targeted PCR metrics and the
# amplicon coverage metrics.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--samples"), dest = "samples_file",
              help = "CSV file containing details of sample datasets (ID and Sample columns required)"),

  make_option(c("--alignment-metrics"), dest = "alignment_metrics_file",
              help = "Alignment metrics from Picard CollectAlignmentSummaryMetrics in collated tabular format"),

  make_option(c("--targeted-pcr-metrics"), dest = "targeted_pcr_metrics_file",
              help = "Targeted PCR metrics from Picard CollectTargetedPcrMetrics in collated tabular format"),

  make_option(c("--amplicon-coverage"), dest = "amplicon_coverage_file",
              help = "Amplicon coverage metrics file"),

  make_option(c("--output-metrics"), dest = "output_metrics_file",
              help = "Amplicon coverage metrics file"),

  make_option(c("--output-report"), dest = "output_report_file",
              help = "Output HTML report")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

samples_file <- opt$samples_file
alignment_metrics_file <- opt$alignment_metrics_file
targeted_pcr_metrics_file <- opt$targeted_pcr_metrics_file
amplicon_coverage_file <- opt$amplicon_coverage_file
output_metrics_file <- opt$output_metrics_file
output_report_file <- opt$output_report_file

if (is.null(samples_file)) stop("Samples file must be specified")
if (is.null(alignment_metrics_file)) stop("Alignment metrics file must be specified")
if (is.null(targeted_pcr_metrics_file)) stop("Targeted PCR metrics file must be specified")
if (is.null(amplicon_coverage_file)) stop("Alignment coverage file must be specified")
if (is.null(output_metrics_file)) stop("Output merged metrics file must be specified")
if (is.null(output_report_file)) stop("Output report file must be specified")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Nozzle.R1)
  library(base64)
})

# function for creating base64-encoded image that can be embedded within the report
to_embedded_image <- function(file, mimetype)
{
  tempFile <- tempfile(pattern = "base64_", fileext = ".txt")
  encode(file, tempFile)
  encoded <- readChar(tempFile, file.info(tempFile)$size)
  file.remove(tempFile)
  encoded <- str_c("data:", mimetype, ";base64,", encoded)
  encoded
}

# sample metadata
samples <- read_csv(samples_file)

expected_columns <- c("ID", "Sample")
missing_columns <- setdiff(expected_columns, colnames(samples))
if (length(missing_columns) > 0) {
  stop("missing columns found in ", samples_file, ": '", str_c(missing_columns, collapse = "', '"), "'")
}

samples <- samples %>%
    select(all_of(expected_columns)) %>%
    arrange(ID)

ids <- samples$ID

# alignment metrics
alignment_metrics <- read_tsv(alignment_metrics_file) %>%
  filter(CATEGORY %in% c("PAIR", "UNPAIRED")) %>%
  select(-any_of(c("SAMPLE", "LIBRARY", "READ_GROUP", "CATEGORY"))) %>%
  arrange(ID)

if (nrow(alignment_metrics) != length(ids)) {
  stop("unexpected number of entries in alignment metrics file")
}

if (any(alignment_metrics$ID != ids)) {
  stop("inconsistent IDs in samples file and alignment metrics file")
}

# targeted PCR metrics
targeted_pcr_metrics <- read_tsv("targeted_pcr_metrics.txt") %>%
  select(-any_of(c("SAMPLE", "LIBRARY", "READ_GROUP"))) %>%
  arrange(ID)

if (nrow(targeted_pcr_metrics) != length(ids)) {
  stop("unexpected number of entries in targeted PCR metrics file")
}

if (any(targeted_pcr_metrics$ID != ids)) {
  stop("inconsistent IDs in samples file and targeted PCR metrics file")
}

if (any(alignment_metrics$TOTAL_READS != targeted_pcr_metrics$TOTAL_READS)) {
  stop("inconsistent TOTAL_READS in alignment metrics and targeted PCR metrics files")
}

if (any(alignment_metrics$PF_READS != targeted_pcr_metrics$PF_READS)) {
  stop("inconsistent PF_READS in alignment metrics and targeted PCR metrics files")
}

# remove duplicated columns from targeted PCR metrics
targeted_pcr_metrics <- targeted_pcr_metrics %>%
  select(-one_of("TOTAL_READS", "PF_READS", "PCT_PF_READS"))

# amplicon coverage
amplicon_coverage <- read_tsv(amplicon_coverage_file)

amplicon_coverage_ids <- amplicon_coverage %>%
  distinct(ID) %>%
  arrange(ID)

if (nrow(amplicon_coverage_ids) != length(ids)) {
  stop("unexpected number of entries in amplicon coverage file")
}

if (any(amplicon_coverage_ids != ids)) {
  stop("inconsistent IDs in samples file and amplicon coverage file")
}

# compute total number of assigned read pairs per dataset
assigned_read_pairs <- amplicon_coverage %>%
  group_by(ID) %>%
  summarize(ASSIGNED_READ_PAIRS = sum(`Read pairs`))

# merge metrics into single data frame
merged_metrics <- samples %>%
  left_join(alignment_metrics, by = "ID") %>%
  left_join(targeted_pcr_metrics, by = "ID") %>%
  left_join(assigned_read_pairs, by = "ID")

# compute derived metrics
merged_metrics <- merged_metrics %>%
  mutate(PCT_ASSIGNED_READ_PAIRS = ifelse(TOTAL_READS == 0, 0, 2 * ASSIGNED_READ_PAIRS / TOTAL_READS)) %>%
  mutate(PCT_ALIGNED_BASES = ifelse(PF_BASES == 0, 0, PF_ALIGNED_BASES / PF_BASES)) %>%
  mutate(PCT_USABLE_READS = ifelse(PF_READS == 0, 0, PF_UQ_READS_ALIGNED / PF_READS)) %>%
  mutate(PCT_USABLE_BASES = ifelse(PF_BASES == 0, 0, PF_UQ_BASES_ALIGNED / PF_BASES)) %>%
  mutate(PCT_USABLE_BASES_ON_AMPLICON = ifelse(PF_BASES == 0, 0, ON_AMPLICON_BASES / PF_BASES)) %>%
  mutate(YIELD = PF_BASES / 1e9) %>%
  mutate(YIELD_ALIGNED = PCT_ALIGNED_BASES * YIELD) %>%
  mutate(YIELD_UNALIGNED = pmax(0, YIELD - YIELD_ALIGNED)) %>%
  mutate(YIELD_USABLE_BASES = PCT_USABLE_BASES * YIELD) %>%
  mutate(YIELD_UNUSABLE = pmax(0, YIELD_ALIGNED - YIELD_USABLE_BASES)) %>%
  mutate(YIELD_USABLE_BASES_ON_AMPLICON = PCT_USABLE_BASES_ON_AMPLICON * YIELD) %>%
  mutate(YIELD_OFF_TARGET = pmax(0, YIELD_USABLE_BASES - YIELD_USABLE_BASES_ON_AMPLICON))

# write merged metrics
if (str_ends(str_to_lower(output_metrics_file), "\\.csv")) {
  write_csv(merged_metrics, output_metrics_file)
} else {
  write_tsv(merged_metrics, output_metrics_file)
}

# correct percentages
merged_metrics <- merged_metrics %>%
  mutate_at(vars(matches("^PCT_")), ~ . * 100) %>%
  mutate_at(vars(matches("_RATE$")), ~ . * 100)

# create labels from ID and SAMPLE fields
merged_metrics <- merged_metrics %>%
  mutate(row = row_number()) %>%
  mutate(label = str_c(ID, Sample, sep = " ")) %>%
  mutate(label = as_factor(label))

# alignment metrics for table
alignment_metrics <- merged_metrics %>%
  select(
    ID,
    Sample,
    `Mean read length` = MEAN_READ_LENGTH,
    `No. reads` = PF_READS,
    `% reads mapped` = PCT_PF_READS_ALIGNED,
    `% bases mapped` = PCT_ALIGNED_BASES,
    `% mismatch rate` = PF_MISMATCH_RATE,
    `% usable reads` = PCT_USABLE_READS,
    `% usable bases` = PCT_USABLE_BASES,
    `% bases off amplicon` = PCT_OFF_AMPLICON,
    `% assigned read pairs` = PCT_ASSIGNED_READ_PAIRS
  )

# yield plot
yield_metrics <- merged_metrics %>%
  select(
    id = label,
    `on-amplicon` = YIELD_USABLE_BASES_ON_AMPLICON,
    `off-amplicon` = YIELD_OFF_TARGET,
    unusable = YIELD_UNUSABLE,
    unaligned = YIELD_UNALIGNED
  ) %>%
  pivot_longer(-id, names_to = "variable", values_to = "value") %>%
  mutate(variable = as_factor(variable)) %>%
  mutate(id = fct_rev(id))

yield_plot <-
  ggplot(data = yield_metrics, aes(x = id, y = value, fill = fct_rev(variable))) +
  geom_bar(stat = "identity", width = 0.66, color = "dodgerblue4") +
  scale_fill_manual(values = rev(alpha(c("dodgerblue", "dodgerblue3", "darkgrey", "white"), 0.75))) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, max(merged_metrics$YIELD) * 1.025)) +
  coord_flip() +
  ylab("Yield (Gbases)") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.position = "bottom",
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_text(size = 6, hjust = 0, vjust = 0.5),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  guides(fill = guide_legend(keywidth = 1, keyheight = 0.25, reverse = TRUE))

yield_plot_file <- tempfile(pattern = "aligned_yield_", fileext = ".svg")
ggsave(yield_plot_file, plot = yield_plot, width = 7, height = max(nrow(merged_metrics) / 10.0, 5.0))

# target coverage plot
mean_target_coverage_plot <-
  ggplot(data = merged_metrics, aes(x = row, y = MEAN_TARGET_COVERAGE, group = 1)) +
  geom_line(color = "chartreuse3") +
  geom_point(color = "chartreuse3") +
  xlab("Samples/Libraries") +
  ylab("Mean target coverage") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 10)
  ) +
  scale_x_continuous(expand = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, max(merged_metrics$MEAN_TARGET_COVERAGE * 1.1)))

mean_target_coverage_plot_file <- tempfile(pattern = "mean_target_coverage_", fileext = ".svg")
ggsave(mean_target_coverage_plot_file, plot = mean_target_coverage_plot, width = 7, height = 6)

# usable bases off-amplicon plot
usable_bases_off_amplicon_plot <-
  ggplot(data = merged_metrics, aes(x = row, y = PCT_OFF_AMPLICON, group = 1)) +
  geom_line(color = "firebrick") +
  geom_point(color = "firebrick") +
  xlab("Samples/Libraries") +
  ylab("% usable bases off amplicon") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 10)
  ) +
  scale_x_continuous(expand = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, min(c(100, max(merged_metrics$PCT_OFF_AMPLICON) * 1.1))))

usable_bases_off_amplicon_plot_file <- tempfile(pattern = "usable_bases_off_amplicon_", fileext = ".svg")
ggsave(usable_bases_off_amplicon_plot_file, plot = usable_bases_off_amplicon_plot, width = 7, height = 6)

# amplicon coverage box plots
amplicon_coverage_plot_data <- amplicon_coverage %>%
  select(Amplicon, `Mean coverage`) %>%
  arrange(Amplicon) %>%
  mutate(Amplicon = factor(Amplicon)) %>%
  mutate(Amplicon = fct_rev(Amplicon)) %>%
  mutate(`Mean coverage` = pmax(`Mean coverage`, 1.0))

amplicon_coverage_plot <-
  ggplot(amplicon_coverage_plot_data, aes(x = Amplicon, y = `Mean coverage`)) +
  geom_boxplot(fill = "darkolivegreen3", outlier.colour = "grey") +
  coord_flip() +
  xlab("Amplicon") +
  ylab("Mean coverage (log scale)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_text(size = 6, hjust = 0, vjust = 0.5),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_log10(expand = c(0.01, 0), breaks = c(10, 100, 1000, 10000), limits = c(1, max(amplicon_coverage_plot_data$`Mean coverage`)))

amplicon_coverage_plot_file <- tempfile(pattern = "amplicon_coverage_", fileext = ".svg")
number_of_amplicons <- amplicon_coverage_plot_data %>% distinct(Amplicon) %>% nrow()
ggsave(amplicon_coverage_plot_file, plot = amplicon_coverage_plot, width = 7, height = max(number_of_amplicons / 10.0, 5.0))

# Nozzle report
report <- newCustomReport("Amplicon Sequencing Alignment and Target Coverage Report")
# report <- setReportSubTitle(report, title)

alignment_metrics_section <- newSection("Alignment Metrics")

alignment_metrics_section <- addTo(alignment_metrics_section, newParagraph("
Alignment metrics were computed using Picard CollectAlignmentSummaryMetrics. The percentages of usable reads and bases were obtained by running Picard CollectTargetedPcrMetrics.
"))

alignment_metrics_section <- addTo(alignment_metrics_section, newParagraph("
Usable reads exclude reads that are unmapped or have a mapping quality of zero indicating that they align equally well to more than one location.
"))

alignment_metrics_table <- newTable(as.data.frame(alignment_metrics), "Alignment Metrics", significantDigits = 3)
alignment_metrics_section <- addTo(alignment_metrics_section, alignment_metrics_table)

yield_section <- newSection("Sample Yields")

yield_section <- addTo(yield_section, newParagraph("
The yield is the total number of bases sequenced per sample, i.e. the number of reads multiplied by the read length. The proportion of bases that align, are usuable and are on target (within a target amplicon) are shown as coloured segments of the bar for each sample in the figure below.
"))

yield_figure <- newFigure(to_embedded_image(yield_plot_file, "image/svg+xml"), "Total yield for each sample/library and the aligned portion of reads.")
yield_section <- addTo(yield_section, yield_figure)

target_coverage_section <- newSection("Target Coverage Metrics")

target_coverage_section <- addTo(target_coverage_section, newParagraph("
Target coverage metrics were computed using Picard CollectTargetedPcrMetrics.
"))

target_coverage_section <- addTo(target_coverage_section, newParagraph("
In computing target coverage metrics, Picard CollectTargetedPcrMetrics does not include reads that are unmapped or have a mapping quality of zero indicating that they align equally well to more than one location. It also excludes supplemental records, e.g. for split read alignments from BWA-MEM.
"))

mean_target_coverage_figure <- newFigure(to_embedded_image(mean_target_coverage_plot_file, "image/svg+xml"), "Mean target coverage for each sample.")
target_coverage_section <- addTo(target_coverage_section, mean_target_coverage_figure)

usable_bases_off_amplicon_figure <- newFigure(to_embedded_image(usable_bases_off_amplicon_plot_file, "image/svg+xml"), "Percentage of usable bases off amplicon.")
target_coverage_section <- addTo(target_coverage_section, usable_bases_off_amplicon_figure)

amplicon_coverage_figure <- newFigure(to_embedded_image(amplicon_coverage_plot_file, "image/svg+xml"), "Mean coverage for each target amplicon across all samples based on reads assigned to amplicons.")
target_coverage_section <- addTo(target_coverage_section, amplicon_coverage_figure)

report <- addTo(report, alignment_metrics_section)
report <- addTo(report, yield_section)
report <- addTo(report, target_coverage_section)

temp_report_file <- tempfile(pattern = "report_", fileext = "")
writeReport(report, filename = temp_report_file, output = HTML.REPORT)
file.copy(str_c(temp_report_file, ".html"), output_report_file)

file.remove(yield_plot_file)
file.remove(mean_target_coverage_plot_file)
file.remove(usable_bases_off_amplicon_plot_file)
file.remove(amplicon_coverage_plot_file)
file.remove(temp_report_file)

