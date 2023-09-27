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

doTransform <- function(signals) {
  # Should become obsolete when fixing performCalibratedMaxSprt()
  pivot <- signals %>%
    tidyr::pivot_longer(cols = starts_with("signal"),
                        names_to = "detectionMethod",
                        values_to = "signal") %>%
    mutate(label = case_when(
      .data$detectionMethod == "signalCalibratedMaxSprt" ~ "Calibrated MaxSPRT",
      .data$detectionMethod == "signalCalibratedP" ~ "Calibrated P",
      .data$detectionMethod == "signalMaxSprt" ~ "MaxSPRT",
      .data$detectionMethod == "signalP" ~ "P"
    )) %>%
    select(-"detectionMethod")
  signals <- pivot %>%
    filter(.data$signal) %>%
    group_by(.data$exposureOutcomeId, .data$label, .data$alpha) %>%
    summarize(lookId = min(.data$lookId))
  signals <- pivot %>%
    distinct(.data$exposureOutcomeId, .data$label, .data$alpha) %>%
    left_join(signals, by = join_by("exposureOutcomeId", "label", "alpha")) %>%
    mutate(lookId = if_else(is.na(.data$lookId), Inf, .data$lookId))
  return(signals)
}

evaluateSignals <- function(signals, simulationSettings) {
  signals <- doTransform(signals)
  signalAll <- tibble(
    exposureOutcomeId = seq_along(simulationSettings$exposureOutcomeSettings),
    label = "Signal all",
    lookId = 0)
  signalNone <- tibble(
    exposureOutcomeId = seq_along(simulationSettings$exposureOutcomeSettings),
    label = "Signal none",
    lookId = Inf)
  groups <- signals %>%
    group_by(.data$alpha) %>%
    group_split()
  computePerAlphaMetrics <- function(group, simulationSettings) {
    group <- group %>%
      bind_rows(signalAll, signalNone)
    metrics <- computeConfusionMatrix(group, simulationSettings) %>%
      left_join(computeAttributableRisk(group, simulationSettings), by = "label") %>%
      mutate(alpha = group$alpha[1])
    # metrics <- computeAttributableRisk(group, simulationSettings) %>%
    #   mutate(alpha = group$alpha[1])
  }
  perAlphaMetrics <- map_dfr(groups, computePerAlphaMetrics, simulationSettings = simulationSettings)

  return(perAlphaMetrics)
}

computeConfusionMatrix <- function(signals, simulationSettings) {
  negativeControlIds <- getNegativeControlIds(simulationSettings)
  labels <- signals %>%
    distinct(.data$label) %>%
    pull()
  fullGrid <- expand.grid(label = labels,
                          exposureOutcomeId = seq_along(simulationSettings$exposureOutcomeSettings),
                          lookId = simulationSettings$looks) %>%
    # lookId = seq_len(simulationSettings$looks)) %>%
    left_join(signals %>%
                rename(signalLookId = "lookId"),
              by = join_by("label", "exposureOutcomeId")) %>%
    mutate(signal = lookId >= signalLookId,
           groundTruth = !.data$exposureOutcomeId %in% negativeControlIds)
  # Probably could replace group_split and lapply with just some more dplyr here:
  groups <- fullGrid %>%
    group_by(.data$label, .data$lookId) %>%
    group_split
  confusionMatrices <- lapply(groups, computeMatrix) %>%
    bind_rows()
  return(confusionMatrices)
}

computeMatrix <- function(group) {
  tp <- sum(group$signal & group$groundTruth)
  fp <- sum(group$signal & !group$groundTruth)
  tn <- sum(!group$signal & !group$groundTruth)
  fn <- sum(!group$signal & group$groundTruth)
  tibble(tp = tp,
         fp = fp,
         tn = tn,
         fn = fn,
         type1 = fp / (fp + tn),
         type2 = fn / (tp + fn),
         fdr = fp / (tp + fp),
         label = group$label[1],
         lookId = group$lookId[1]) %>%
    return()
}

computeBaselineAttributableRisk <- function(simulationSettings) {
  attributableRiskPerExposureOutcome <- map_dbl(simulationSettings$exposureOutcomeSettings,
                                                computeAttributableRiskPerExposureOutcome,
                                                simulationSettings = simulationSettings)
  # Note: current simulation assumes equal number of exposures at each look, so attributable
  # risk after n looks is just n * (attributable risk in one look):
  attributableRiskPerExposureOutcome <- tibble(
    attributableRisk = attributableRiskPerExposureOutcome,
    exposureOutcomeId = seq_along(attributableRiskPerExposureOutcome)) %>%
    cross_join(tibble(lookId = c(0, seq_len(simulationSettings$looks)))) %>%
    mutate(attributableRisk = attributableRisk * (simulationSettings$looks - .data$lookId))
  return(attributableRiskPerExposureOutcome)
}

computeAttributableRisk <- function(signals, simulationSettings) {
  attributableRiskPerExposureOutcome <- computeBaselineAttributableRisk(
    simulationSettings = simulationSettings
  )
  attributableRisk <- signals %>%
    filter(!is.infinite(.data$lookId)) %>%
    inner_join(attributableRiskPerExposureOutcome, by = join_by("exposureOutcomeId", "lookId")) %>%
    group_by(.data$label) %>%
    summarise(
      attributableRisk = sum(.data$attributableRisk),
      .groups = "drop"
    )
  return(attributableRisk)
}

computeAttributableRiskPerExposureOutcome <- function(exposureOutcomeSetting, simulationSettings) {
  sumDbMultipliers <- sum(unlist(ParallelLogger::selectFromList(simulationSettings$databaseSettings, "sampleSizeMultiplier")))
  negative <- exposureOutcomeSetting$logRrMean == 0 & exposureOutcomeSetting$logRrSd == 0
  exposedCases <- exposureOutcomeSetting$nTarget *
    exposureOutcomeSetting$backgroundRate *
    (exposureOutcomeSetting$riskEnd - exposureOutcomeSetting$riskStart + 1) *
    sumDbMultipliers *
    exp(exposureOutcomeSetting$logRrMean)
  if (negative) {
    attributableRisk <- 0
  } else {
    attributableRisk <- exposedCases - (exposedCases/exp(exposureOutcomeSetting$logRrMean))
  }
  return(attributableRisk)
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
plotDecisionCurves <- function(evaluation,
                               labels = NULL,
                               alphas = NULL,
                               fileName = NULL) {
  if (is.null(labels)) {
    labels <- unique(evaluation$label)
  }

  if (is.null(alphas)) {
    alphas <- unique(evaluation$alpha)
  }
  # evaluation %>%
  #   filter(label == "Calibrated MaxSPRT") %>%
  #   mutate(p = 1 - (.data$alpha / nExposureOutcomes)) %>%
  #   mutate(odds = p / (1-p))

  nExposureOutcomes <- evaluation %>%
    mutate(n = .data$tp + .data$fp + .data$tn + .data$fn) %>%
    summarize(n = median(n)) %>%
    pull()
  # plotData <- evaluation %>%
  #   filter(.data$iteration == 1) %>%
  #   mutate(p = 1 - (.data$alpha / nExposureOutcomes)) %>%
  #   mutate(netBenefit = (.data$tp / nExposureOutcomes) - (.data$fp / nExposureOutcomes) * (p / (1-p)))
  plotData <- evaluation %>%
    filter(.data$iteration == 1, abs(.data$alpha - 5) < 0.01) %>%
    cross_join(tibble(p = seq(0, 1, by = 0.01))) %>%
    mutate(netBenefit = (.data$tp / nExposureOutcomes) - (.data$fp / nExposureOutcomes) * (p / (1-p)))
  maxBenefit <- max(plotData$netBenefit, na.rm = TRUE)
  plot <- ggplot2::ggplot(plotData, ggplot2::aes(x = p, y = .data$netBenefit, group = .data$label, color = .data$label)) +
    ggplot2::geom_line(size = 2) +
    ggplot2::scale_x_continuous("P") +
    ggplot2::scale_y_continuous("Net benefit", limits = c(-maxBenefit/5, maxBenefit))
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
