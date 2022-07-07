# Copyright 2022 Observational Health Data Sciences and Informatics
#
# This file is part of DiscoverySystemSimulator
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# library(dplyr)
# library(purrr)
# source("C:/Users/mschuemi/git/DiscoverySystemSimulator/R/SimulationSettings.R")

evaluateSignals <- function(signals, simulationSettings) {
  groups <- signals %>%
    group_by(.data$alpha) %>%
    group_split()

  computePerAlphaMetrics <- function(group, simulationSettings) {
    metrics <- computeConfusionMatrix(group, simulationSettings) %>%
      inner_join(computeUtility(group, simulationSettings), by = "label") %>%
      mutate(alpha = group$alpha[1])
  }
  perAlphaMetrics <- map_dfr(groups, computePerAlphaMetrics, simulationSettings = simulationSettings)

  return(perAlphaMetrics)
}



computeConfusionMatrix <- function(signals, simulationSettings) {
  negativeControlIds <- getNegativeControlIds(simulationSettings)

  signalsPerExposureOutcome <- signals %>%
    group_by(.data$exposureOutcomeId) %>%
    summarise(signalMaxSprt = any(.data$signalMaxSprt),
              signalCalibratedMaxSprt = any(.data$signalCalibratedMaxSprt),
              signalP = any(.data$signalP),
              signalCalibratedP = any(.data$signalCalibratedP),
              .groups = "drop") %>%
    mutate(groundTruth = !.data$exposureOutcomeId %in% negativeControlIds)

  confusionMatrix <- bind_rows(computeMatrix(signalsPerExposureOutcome$signalCalibratedMaxSprt,
                                             signalsPerExposureOutcome$groundTruth,
                                             "Calibrated MaxSPRT"),
                               computeMatrix(signalsPerExposureOutcome$signalMaxSprt,
                                             signalsPerExposureOutcome$groundTruth,
                                             "MaxSPRT"),
                               computeMatrix(signalsPerExposureOutcome$signalP,
                                             signalsPerExposureOutcome$groundTruth,
                                             "P"),
                               computeMatrix(signalsPerExposureOutcome$signalCalibratedP,
                                             signalsPerExposureOutcome$groundTruth,
                                             "calibratedP"))
  return(confusionMatrix)
}

computeMatrix <- function(signal, groundTruth, label) {
  tp <- sum(signal & groundTruth)
  fp <- sum(signal & !groundTruth)
  tn <- sum(!signal & !groundTruth)
  fn <- sum(!signal & groundTruth)
  tibble(tp = tp,
         fp = fp,
         tn = tn,
         fn = fn,
         type1 = fp / (fp + tn),
         type2 = fn / (tp + fn),
         fdr = fp / (tp + fp),
         label = label) %>%
    return()
}

computeUtility <- function(signals, simulationSettings) {
  utilityPerExposureOutcome <- map_dfr(simulationSettings$exposureOutcomeSettings, computeUtilityPerExposureOutcome, simulationSettings = simulationSettings) %>%
    mutate(exposureOutcomeId = row_number())

  signalsPerExposureOutcome <- signals %>%
    group_by(.data$exposureOutcomeId) %>%
    summarise(
      signalMaxSprt = any(.data$signalMaxSprt),
      signalCalibratedMaxSprt = any(.data$signalCalibratedMaxSprt),
      signalP = any(.data$signalP),
      signalCalibratedP = any(.data$signalCalibratedP),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = starts_with("signal"),
    )


  utility <- signalsPerExposureOutcome %>%
    inner_join(utilityPerExposureOutcome, by = "exposureOutcomeId") %>%
    mutate(utility = if_else(.data$value, .data$uPositive, .data$uNegative)) %>%
    group_by(.data$name) %>%
    summarise(
      utility = sum(utility),
      .groups = "drop"
    ) %>%
    mutate(label = case_when(
      .data$name == "signalCalibratedMaxSprt" ~ "Calibrated MaxSPRT",
      .data$name == "signalCalibratedP" ~ "calibratedP",
      .data$name == "signalMaxSprt" ~ "MaxSPRT",
      .data$name == "signalP" ~ "P"
    )) %>%
    select(-.data$name)
  return(utility)
}

computeUtilityPerExposureOutcome <- function(exposureOutcomeSetting, simulationSettings) {
  sumDbMultipliers <- sum(unlist(ParallelLogger::selectFromList(simulationSettings$databaseSettings, "sampleSizeMultiplier")))
  negative <- exposureOutcomeSetting$logRrMean == 0 & exposureOutcomeSetting$logRrSd == 0
  exposedCases <- exposureOutcomeSetting$nTarget *
    exposureOutcomeSetting$backgroundRate *
    (exposureOutcomeSetting$riskEnd - exposureOutcomeSetting$riskStart + 1) *
    sumDbMultipliers *
    exp(exposureOutcomeSetting$logRrMean)
  if (negative) {
    fictitiousRr <- 2
    utility <- tibble(
      uPositive = -(exposedCases - (exposedCases/fictitiousRr)),
      uNegative = 0
    )
  } else {
    utility <- tibble(
      uPositive = 0,
      uNegative = -(exposedCases - (exposedCases/exp(exposureOutcomeSetting$logRrMean)))
    )
  }
  return(utility)
}


#' Plot distribution of false positives and negatives
#'
#' @param evaluation        An object generated by [evaluateIterations()].
#' @param labels            Which labels to plot? (Values of the `label` column). If `NULL` then all values are plotted.
#' @param cumulative        Should cumulative probabilities be shown?
#' @param alphas            Which alphas to plot? (Values of the `alpha` column). If `NULL` then all alphas are plotted.
#' @param fileName          Optional: the name of the file to save the plot to.
#'
#' @return
#' A GGPlot object.
#'
#' @export
plotFalsePositiveNegatives <- function(evaluation,
                                       labels = NULL,
                                       cumulative = TRUE,
                                       alphas = NULL,
                                       fileName = NULL) {
  if (is.null(labels)) {
    labels <- unique(evaluation$label)
  }

  if (is.null(alphas)) {
    alphas <- unique(evaluation$alpha)
  }

  nIterations <- length(unique(evaluation$iteration))

  plotData <- evaluation %>%
    filter(.data$label %in% labels & .data$alpha %in% alphas) %>%
    select(.data$fp, .data$fn, .data$label, .data$alpha, .data$iteration) %>%
    tidyr::pivot_longer(cols = c("fp", "fn"), values_to = "count", names_to = "type") %>%
    group_by(.data$label, .data$alpha, .data$type, .data$count) %>%
    summarize(percentage = 100* n() / nIterations, .groups = "drop") %>%
    mutate(type = recode(.data$type, fp = "False positive", fn = "False negative"))

  if (cumulative) {
    plotData <- plotData %>%
      filter(.data$count != 0)

    toAdd <- plotData %>%
      group_by(.data$label, .data$alpha, .data$type) %>%
      summarize(minCount = min(.data$count) - 1, .groups = "drop") %>%
      filter(.data$minCount > 0)
    if (nrow(toAdd) > 0) {
      for (i in 1:nrow(toAdd)) {
        plotData <- bind_rows(tibble(label = toAdd$label[i],
                                     alpha = toAdd$alpha[i],
                                     type = toAdd$type[i],
                                     count = 1:toAdd$minCount[i],
                                     percentage = 0),
                              plotData)
      }
    }

    plotData <- plotData %>%
      group_by(.data$label, .data$alpha, .data$type) %>%
      arrange(desc(.data$count)) %>%
      mutate(percentage = cumsum(.data$percentage)) %>%
      ungroup()
    yLabel <- "Cumulative probability (%)"
  } else {
    yLabel <- "Probability (%)"
  }
  if (length(unique(plotData$label)) > 1) {
    if (length(unique(plotData$alpha)) > 1) {
      plotData$text <- sprintf("%s, alpha = %0.2f", plotData$label, plotData$alpha)
    } else {
      plotData$text <- sprintf("%s", plotData$label)
    }
  } else {
    plotData$text <- sprintf("alpha = %0.2f", plotData$alpha)
  }

  plotData$text <- factor(plotData$text, levels = sort(unique(plotData$text)))
  plot <- ggplot2::ggplot(plotData, ggplot2::aes(x = as.factor(.data$count), y = .data$percentage)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_x_discrete("Number of misclassified estimates") +
    ggplot2::scale_y_continuous(yLabel) +
    ggplot2::facet_grid(text ~ type, scales = "free_x")

  if (!is.null(fileName)) {
    ggplot2::ggsave(plot = plot,
                    filename = fileName,
                    width = 6,
                    height = 1.25 + length(unique(plotData$text)) * 1.25,
                    dpi =  200)
  }

  return(plot)
}

#' Plot false discovery rate
#'
#' @param evaluation        An object generated by [evaluateIterations()].
#' @param labels            Which labels to plot? (Values of the `label` column). If `NULL` then all values are plotted.
#' @param alphas            Which alphas to plot? (Values of the `alpha` column). If `NULL` then all alphas are plotted.
#' @param fileName          Optional: the name of the file to save the plot to.
#'
#' @return
#' A GGPlot object.
#'
#' @export
plotFalseDiscoveryRate <- function(evaluation,
                                   labels = NULL,
                                   alphas = NULL,
                                   fileName = NULL) {
  if (is.null(labels)) {
    labels <- unique(evaluation$label)
  }

  if (is.null(alphas)) {
    alphas <- unique(evaluation$alpha)
  }

  plotData <- evaluation %>%
    filter(.data$label %in% labels & .data$alpha %in% alphas) %>%
    mutate(fdr = .data$fp / (.data$fp + .data$tp)) %>%
    select(.data$fdr, .data$label, .data$alpha, .data$iteration) %>%
    mutate(alpha = as.factor(sprintf("%0.2f", .data$alpha)))

  plot <- ggplot2::ggplot(plotData, ggplot2::aes(x = as.factor(.data$alpha), y = .data$fdr)) +
    ggplot2::geom_violin(scale = "width", color = "#e85847", fill = "#e85847", alpha = 0.5) +
    ggplot2::scale_x_discrete("Alpha") +
    ggplot2::scale_y_continuous("False discovery rate") +
    ggplot2::facet_grid(.data$label~.)

  if (!is.null(fileName)) {
    ggplot2::ggsave(plot = plot,
                    filename = fileName,
                    width = 2 + length(unique(plotData$alpha)) * 0.5,
                    height = 1.25 + length(unique(plotData$label)) * 1.25,
                    dpi =  200)
  }

  return(plot)
}

#' Plot receiver operator curve (ROC)
#'
#' @param evaluation        An object generated by [evaluateIterations()].
#' @param labels            Which labels to plot? (Values of the `label` column). If `NULL` then all values are plotted.
#' @param alphas            Which alphas to plot? (Values of the `alpha` column). If `NULL` then all alphas are plotted.
#' @param fileName          Optional: the name of the file to save the plot to.
#'
#' @return
#' A GGPlot object.
#'
#' @export
plotRoc <- function(evaluation,
                    labels = NULL,
                    alphas = NULL,
                    fileName = NULL) {
  if (is.null(labels)) {
    labels <- unique(evaluation$label)
  }

  if (is.null(alphas)) {
    alphas <- unique(evaluation$alpha)
  }

  nIterations <- length(unique(evaluation$iteration))

  plotData <- evaluation %>%
    filter(.data$label %in% labels & .data$alpha %in% alphas) %>%
    mutate(fpr = .data$fp / (.data$fp + .data$tn),
           tpr = .data$tp / (.data$fn + .data$tp)) %>%
    select(.data$fpr, .data$tpr, .data$label, .data$alpha, .data$iteration) %>%
    mutate(alpha = as.factor(sprintf("%0.2f", .data$alpha)))


  plot <- ggplot2::ggplot(plotData, ggplot2::aes(x = .data$fpr, y = .data$tpr, group = .data$iteration)) +
    ggplot2::geom_abline(slope = 1, linetype = "dashed") +
    ggplot2::geom_line(color = "#444444", alpha = 0.3) +
    ggplot2::geom_point(ggplot2::aes(color = .data$alpha), alpha = 0.3) +
    ggplot2::scale_x_continuous("False positive rate", limits = c(0,1)) +
    ggplot2::scale_y_continuous("True positive rate", limits = c(0,1)) +
    ggplot2::facet_grid(label~.)

  if (!is.null(fileName)) {
    ggplot2::ggsave(plot = plot,
                    filename = fileName,
                    width = 4.5,
                    height = 2 + length(unique(plotData$label)) * 1,
                    dpi =  200)
  }

  return(plot)
}


#' Plot utility
#'
#' @param evaluation        An object generated by [evaluateIterations()].
#' @param labels            Which labels to plot? (Values of the `label` column). If `NULL` then all values are plotted.
#' @param alphas            Which alphas to plot? (Values of the `alpha` column). If `NULL` then all alphas are plotted.
#' @param fileName          Optional: the name of the file to save the plot to.
#'
#' @return
#' A GGPlot object.
#'
#' @export
plotUtility <- function(evaluation,
                        labels = NULL,
                        alphas = NULL,
                        fileName = NULL) {
  if (is.null(labels)) {
    labels <- unique(evaluation$label)
  }

  if (is.null(alphas)) {
    alphas <- unique(evaluation$alpha)
  }

  plotData <- evaluation %>%
    filter(.data$label %in% labels & .data$alpha %in% alphas) %>%
    select(.data$utility  , .data$label, .data$alpha, .data$iteration) %>%
    mutate(alpha = as.factor(sprintf("%0.2f", .data$alpha)))

  plot <- ggplot2::ggplot(plotData, ggplot2::aes(x = as.factor(.data$alpha), y = .data$utility  )) +
    ggplot2::geom_violin(scale = "width", color = "#e85847", fill = "#e85847", alpha = 0.5) +
    ggplot2::scale_x_discrete("Alpha") +
    ggplot2::scale_y_continuous("Utility") +
    ggplot2::facet_grid(.data$label~.)
plot
  if (!is.null(fileName)) {
    ggplot2::ggsave(plot = plot,
                    filename = fileName,
                    width = 2 + length(unique(plotData$alpha)) * 0.5,
                    height = 1.5 + length(unique(plotData$label)) * 1.4,
                    dpi =  200)
  }

  return(plot)
}
