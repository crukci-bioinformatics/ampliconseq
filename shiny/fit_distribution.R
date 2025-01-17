library(fitdistrplus)
library(dplyr)

fitDistribution <- function(
  alleleFractions,
  excludeHighestProportion = 0.025,
  maximumAlleleFraction = 1.0,
	includeZeroAlleleFractionValues = FALSE,
  distribution = "Beta",
  thresholdProbabilities = c(0.99, 0.999, 0.999),
  numberOfPoints = 1024,
  minimumNumberOfAlleleFractions = 5)
{
  result <- list(
    alleleFractions = alleleFractions,
    density = NULL,
    filteredAlleleFractions = NULL,
    filteredDensity = NULL,
    maximumAlleleFraction = NULL,
    fitted = NULL,
    thresholds = NULL
  )

  if (is.null(alleleFractions) || length(alleleFractions) == 0) return(result)

  result$maximumAlleleFraction <- min(max(alleleFractions) * 1.1, 1.0)

  density <- density(alleleFractions, n = numberOfPoints)
  result$density <- tibble(`Allele fraction` = density$x, Density = density$y)

  filteredAlleleFractions <- alleleFractions

  if (excludeHighestProportion > 0)
  {
    filteredAlleleFractions <- sort(filteredAlleleFractions)
    n <- length(alleleFractions)
    m <- as.integer(excludeHighestProportion * n)
    m <- n - min(m, n)
    filteredAlleleFractions <- filteredAlleleFractions[1:m]
  }

	if (!includeZeroAlleleFractionValues) filteredAlleleFractions <- filteredAlleleFractions[filteredAlleleFractions != 0]

  if (!is.na(maximumAlleleFraction) && !is.null(maximumAlleleFraction))
    filteredAlleleFractions <- filteredAlleleFractions[filteredAlleleFractions <= maximumAlleleFraction]

  result$filteredAlleleFractions <- filteredAlleleFractions

  density <- density(filteredAlleleFractions, n = numberOfPoints)
  result$filteredDensity <- tibble(`Allele fraction` = density$x, Density = density$y)

  if (length(filteredAlleleFractions) < minimumNumberOfAlleleFractions) return(result)

  if (distribution == "Normal")
  {
    fit <- fitdist(filteredAlleleFractions, "norm", method = "mme")
    result$thresholds <- qnorm(thresholdProbabilities, mean = fit$estimate[1], sd = fit$estimate[2])
    result$maximumAlleleFraction <- min(max(c(result$maximumAlleleFraction, result$thresholds)) * 1.1, 1.0)
    x <- seq(0.0, result$maximumAlleleFraction, length.out = numberOfPoints)
    fitted <- dnorm(x, mean = fit$estimate[1], sd = fit$estimate[2])
    result$fitted <- tibble(`Allele fraction` = x, Density = fitted)
  }
  else if (distribution == "Log-normal")
  {
    fit <- fitdist(filteredAlleleFractions, "lnorm", method = "mge")
    result$thresholds <- qlnorm(thresholdProbabilities, meanlog = fit$estimate[1], sdlog = fit$estimate[2])
    result$maximumAlleleFraction <- min(max(c(result$maximumAlleleFraction, result$thresholds)) * 1.1, 1.0)
    x <- seq(0.0, result$maximumAlleleFraction, length.out = numberOfPoints)
    fitted <- dlnorm(x, meanlog = fit$estimate[1], sdlog = fit$estimate[2])
    result$fitted <- tibble(`Allele fraction` = x, Density = fitted)
  }
  else if (distribution == "Beta")
  {
    fit <- fitdist(filteredAlleleFractions, "beta", method = "mge")
    result$thresholds <- qbeta(thresholdProbabilities, shape1 = fit$estimate[1], shape2 = fit$estimate[2])
    result$maximumAlleleFraction <- min(max(c(result$maximumAlleleFraction, result$thresholds)) * 1.1, 1.0)
    x <- seq(0.0, result$maximumAlleleFraction, length.out = numberOfPoints)
    fitted <- dbeta(x, shape1 = fit$estimate[1], shape2 = fit$estimate[2])
    result$fitted <- tibble(`Allele fraction` = x, Density = fitted)
  }
  else if (distribution == "Gamma")
  {
    fit <- fitdist(filteredAlleleFractions, "gamma", method = "mge")
    result$thresholds <- qgamma(thresholdProbabilities, shape = fit$estimate[1], rate = fit$estimate[2])
    result$maximumAlleleFraction <- min(max(c(result$maximumAlleleFraction, result$thresholds)) * 1.1, 1.0)
    x <- seq(0.0, result$maximumAlleleFraction, length.out = numberOfPoints)
    fitted <- dgamma(x, shape = fit$estimate[1], rate = fit$estimate[2])
    result$fitted <- tibble(`Allele fraction` = x, Density = fitted)
  }
  else if (distribution == "Weibull")
  {
    fit <- fitdist(filteredAlleleFractions, "weibull", method = "mge")
    result$thresholds <- qweibull(thresholdProbabilities, fit$estimate[1], fit$estimate[2])
    result$maximumAlleleFraction <- min(max(c(result$maximumAlleleFraction, result$thresholds)) * 1.1, 1.0)
    x <- seq(0.0, result$maximumAlleleFraction, length.out = numberOfPoints)
    fitted <- dweibull(x, fit$estimate[1], fit$estimate[2])
    result$fitted <- tibble(`Allele fraction` = x, Density = fitted)
  }

  return(result)
}
