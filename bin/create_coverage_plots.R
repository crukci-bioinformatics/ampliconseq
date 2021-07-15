#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
#
# Licensed under the MIT license (http://opensource.org/licenses/MIT).
# This file may not be copied, modified, or distributed except according
# to those terms.
# ------------------------------------------------------------------------------

# Creates coverage and yield plots.

suppressPackageStartupMessages(library(optparse))

option_list <- list(

  make_option(c("--alignment-metrics"), dest = "alignment_metrics_file",
              help = "Coverage metrics file"),

  make_option(c("--amplicon-coverage"), dest = "amplicon_coverage_file",
              help = "Coverage metrics file"),

  make_option(c("--output-prefix"), dest = "output_prefix", default = "",
              help = "Prefix for output plot files")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)
opt <- parse_args(option_parser)

alignment_metrics_file <- opt$alignment_metrics_file
amplicon_coverage_file <- opt$amplicon_coverage_file
output_prefix <- opt$output_prefix

if (is.null(alignment_metrics_file)) stop("Alignment coverage metrics file must be specified")
if (is.null(amplicon_coverage_file)) stop("Amplicon coverage file must be specified")

suppressPackageStartupMessages(library(tidyverse))

yield_plot_prefix <- str_c(output_prefix, "yield")
amplicon_coverage_plot_prefix <- str_c(output_prefix, "amplicon_coverage")

alignment_metrics <- read_csv(alignment_metrics_file, col_types = cols(ID = "c", Sample = "c", .default = "d"))
amplicon_coverage <- read_tsv(amplicon_coverage_file, col_types = cols(ID = "c", Sample = "c", Amplicon = "c", Chromosome = "c", .default = "d"))

# yield stacked bar plot

yield_metrics <- alignment_metrics %>%
  mutate(row = row_number()) %>%
  mutate(label = str_c(ID, Sample, sep = " ")) %>%
  mutate(label = as_factor(label)) %>%
  transmute(
    id = label,
    `on target (usable)` = `Bases on target (usable)`,
    `on target (filtered)` = `Bases on target` - `Bases on target (usable)`,
    `on amplicon (primer)` = `Bases on amplicon` - `Bases on target`,
    `near amplicon` = `Bases near amplicon`,
    `off amplicon` = `Bases off amplicon`,
    unaligned = `Bases` - `Bases aligned`,
    total = `Bases`
  ) %>%
  mutate(across(-id, ~ . / 1e9))

max_yield <- max(yield_metrics$total)

yield_metrics <- yield_metrics %>%
  select(-total) %>%
  pivot_longer(-id, names_to = "category", values_to = "yield") %>%
  mutate(category = as_factor(category)) %>%
  mutate(id = fct_rev(id))

yield_plot <- yield_metrics %>%
  ggplot(aes(x = id, y = yield, fill = fct_rev(category))) +
  geom_bar(stat = "identity", width = 0.66, color = "dodgerblue4") +
  scale_fill_manual(values = rev(alpha(c("dodgerblue", "dodgerblue3", "darkgrey", "orange", "red", "white"), 0.75))) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, max_yield * 1.025)) +
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

yield_plot_height <- nlevels(yield_metrics$id) / 7.5 + 1.0
ggsave(str_c(yield_plot_prefix, ".svg"), plot = yield_plot, width = 8, height = yield_plot_height, limitsize = FALSE)
ggsave(str_c(yield_plot_prefix, ".pdf"), plot = yield_plot, width = 8, height = yield_plot_height, limitsize = FALSE)

# amplicon coverage box plots

amplicon_coverage <- amplicon_coverage %>%
  select(Amplicon, `Mean coverage`) %>%
  arrange(Amplicon) %>%
  mutate(Amplicon = factor(Amplicon)) %>%
  mutate(Amplicon = fct_rev(Amplicon)) %>%
  mutate(`Mean coverage` = pmax(`Mean coverage`, 1.0))

amplicon_coverage_plot <- amplicon_coverage %>%
  ggplot(aes(x = Amplicon, y = `Mean coverage`)) +
  geom_boxplot(fill = "dodgerblue", outlier.colour = "grey") +
  scale_y_log10(
    expand = c(0.01, 0),
    breaks = c(10, 100, 1000, 10000, 100000),
    labels = scales::format_format(big.mark = ",", scientific = FALSE),
    limits = c(1, 1.025 * max(amplicon_coverage$`Mean coverage`))
  ) +
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
  )

number_of_amplicons <- amplicon_coverage %>% distinct(Amplicon) %>% nrow()

amplicon_coverage_plot_height <- number_of_amplicons / 7.5 + 0.5

ggsave(str_c(amplicon_coverage_plot_prefix, ".svg"), plot = amplicon_coverage_plot, width = 8, height = amplicon_coverage_plot_height, limitsize = FALSE)
ggsave(str_c(amplicon_coverage_plot_prefix, ".pdf"), plot = amplicon_coverage_plot, width = 8, height = amplicon_coverage_plot_height, limitsize = FALSE)


