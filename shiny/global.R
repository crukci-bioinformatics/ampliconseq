library(fitdistrplus)
library(tidyr)
library(dplyr)
library(shiny)
library(DT)
library(highcharter)

source("plots.R")
source("fit_distribution.R")

options(shiny.maxRequestSize = 1024 * 1024 * 1024)

thresholdProbabilities <- c(0.99, 0.999, 0.9999)
initialThresholdProbability <- last(thresholdProbabilities)

bases <- c("A", "C", "G", "T")
substitutions <- c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G")

readCountColumns <- c("ID", "Sample", "Amplicon", "Chromosome", "Position", "Reference base", "Depth unfiltered", "Depth", "A count", "C count", "G count", "T count")

emptyReadCounts <- tibble(
  ID = character(),
  Sample = character(),
  Amplicon = character(),
  Chromosome = character(),
  Position = integer(),
  `Reference base` = character(),
  `Depth unfiltered` = integer(),
  `Depth` = integer(),
  `A count` = integer(),
  `C count` = integer(),
  `G count` = integer(),
  `T count` = integer()
)

snvColumns <- c("ID", "Amplicon", "Chromosome", "Position", "Reference base", "Alternate allele", "Filters", "Confidence")

emptySnvs <- tibble(
  ID = character(),
  Amplicon = character(),
  Chromosome = character(),
  Position = integer(),
  `Reference base` = character(),
  `Alternate allele` = character(),
  Filters = character(),
  Confidence = character()
)

