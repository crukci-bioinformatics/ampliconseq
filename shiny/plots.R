library(tidyr)
library(dplyr)
library(stringr)
library(highcharter)

# allele fraction density plot

# input data should be a data frame with the following columns:
#   Allele fraction
#   Density
#   Series

alleleFractionDensityPlot <- function(
  data,
  series = NULL,
  type = NULL,
  colours = c("#7cb5ec", "#434348", "#90ed7d", "#f7a35c", "#8085e9",
              "#f15c80", "#e4d354", "#2b908f", "#f45b5b", "#91e8e1"),
  maximumAlleleFraction = NULL,
  thresholdAlleleFraction = -999,
  thresholdLabel = "",
  decimalPlaces = 3
)
{
  hc <- highchart() %>%
    hc_chart(animation = TRUE, zoomType = "x")

  if (!is.data.frame(data) || nrow(data) == 0) return(hc)
  if (length(setdiff(c("Allele fraction", "Density", "Series"), colnames(data))) > 0) return(hc)

  if (is.null(series))
    series <- data %>% select(Series) %>% distinct %>% arrange(Series) %>% unlist() %>% as.character()
  if (!is.vector(series) || length(series) == 0) return(hc)

  hc <- hc %>%
    hc_xAxis(
      title = list(text = "allele fraction"),
      min = 0.0, max = maximumAlleleFraction,
      plotLines = list(list(
        color = "#FF0000",
        width = 2,
        value = thresholdAlleleFraction,
        label = list(text = str_c("<br>", thresholdLabel), useHTML = TRUE, rotation = 0, textAlign = "center")
      ))
    )

  hc <- hc %>% hc_yAxis(
    title = list(text = "density"),
    min = 0, max = NULL, minRange = 0.9
  )

  for (i in 1:length(series))
  {
    hc <- hc %>% hc_add_series(
      data = data %>%
        filter(Series == series[i]) %>%
        select(-Series) %>%
        rename(x = `Allele fraction`) %>%
        rename(y = Density) %>%
        list_parse,
      name = series[i],
      type = ifelse(is.null(type[i]) || is.na(type[i]), "spline", type[i]),
      color = colours[i],
      marker = list(enabled = FALSE)
    )
  }

  hc <- hc %>%
    hc_legend(layout = "vertical", align = "right", verticalAlign = "middle", floating = TRUE, borderWidth = 1)

  hc <- hc %>%
    hc_tooltip(followPointer = TRUE, crosshairs = TRUE, formatter = JS(str_c("function() { return (this.x.toFixed(", decimalPlaces, ")) }")))

  hc
}


# allele fraction scatter/box plot

alleleFractionScatterBoxPlot <- function(
  data,
  series = NULL,
  xlabel = NULL,
  colours = c("#7cb5ec", "#90ed7d", "#f7a35c", "#8085e9",
              "#f15c80", "#e4d354", "#2b908f", "#f45b5b", "#91e8e1"),
  maximumAlleleFraction = NULL,
  thresholdAlleleFraction = -999,
  thresholdLabel = "",
  clicked = NULL
)
{
  hc <- highchart() %>%
    hc_chart(animation = TRUE, zoomType = "y")

  if (!is.data.frame(data) || nrow(data) == 0) return(hc)
  if (length(setdiff(c("ID", "Allele fraction", "Series", "tooltip"), colnames(data))) > 0) return(hc)

  if (is.null(series))
    series <- data %>% select(Series) %>% distinct %>% arrange(Series) %>% unlist() %>% as.character()
  if (!is.vector(series) || length(series) == 0) return(hc)

  hc <- hc %>%
    hc_xAxis(
      title = list(text = xlabel),
      labels = list(enabled = FALSE),
      min = -1.02, max = 1.02,
      tickLength = 0
    )

  hc <- hc %>%
    hc_yAxis(
      title = list(text = "Allele fraction"),
      min = 0,
      max = maximumAlleleFraction,
      plotLines = list(list(
        color = "#FF0000",
        width = 2,
        value = thresholdAlleleFraction,
        label = list(text = thresholdLabel)
      ))
    )

  stats <- data %>%
    select(`Allele fraction`) %>%
    unlist %>%
    summary %>%
    as.numeric
  stats <- stats[-4]

  hc <- hc %>%
    hc_add_series(
      data = list(stats),
      type = "boxplot",
      name = "box plot",
      color = hex_to_rgba("#434348", alpha = 1.0),
      fillColor = hex_to_rgba("#FFFFFF", alpha = 0.0),
      showInLegend = FALSE,
      tooltip = list(headerFormat = "", valueDecimals = 4)
    )

  data <- data %>%
    mutate(x = seq(-1.0, 1.0, length.out = nrow(data)))

  for (i in 1:length(series))
  {
    seriesData <- data %>%
      filter(Series == series[i]) %>%
      select(-Series)
    if (nrow(seriesData) == 0) next

    colour <- colours[i]
    if (is.na(colour)) colour <- NULL

    hc <- hc %>% hc_add_series(
      data = seriesData %>%
        rename(y = `Allele fraction`) %>%
        list_parse,
      type = "scatter",
      name = series[i],
      color = colour,
      marker = list(symbol = "circle"),
      showInLegend = TRUE,
      tooltip = list(pointFormatter = JS("function() { return (this.tooltip) }"))
    )
  }

  if (!is.null(clicked))
  {
    fn <- str_c("function() { Shiny.onInputChange('", clicked, "', this.ID) }")
    hc <- hc %>%
      hc_plotOptions(
        series = list(
          cursor = "pointer",
          point = list(
            events = list(
              click = JS(fn)
            )
          )
        )
      )
  }

  hc
}
