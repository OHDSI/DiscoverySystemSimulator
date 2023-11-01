# Copyright 2023 Observational Health Data Sciences and Informatics
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
    # metrics <- computeConfusionMatrix(group, simulationSettings) %>%
    #   left_join(computeAttributableRisk(group, simulationSettings), by = "label") %>%
    #   mutate(alpha = group$alpha[1])
    metrics <- computeAttributableRisk(group, simulationSettings) %>%
      mutate(alpha = group$alpha[1])
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

computeAttributableRiskPerExposureOutcome <- function(exposureOutcomeSetting, simulationSettings) {
  sumDbMultipliers <- sum(unlist(ParallelLogger::selectFromList(simulationSettings$databaseSettings, "sampleSizeMultiplier")))
  meanLogRr <- mean(unlist(ParallelLogger::selectFromList(simulationSettings$exposureOutcomeSettings, "logRrMean")))
  negative <- exposureOutcomeSetting$logRrMean == 0 & exposureOutcomeSetting$logRrSd == 0
  exposedCases <- exposureOutcomeSetting$nTarget *
    exposureOutcomeSetting$backgroundRate *
    (exposureOutcomeSetting$riskEnd - exposureOutcomeSetting$riskStart + 1) *
    sumDbMultipliers *
    exp(exposureOutcomeSetting$logRrMean)
  if (negative) {
    tibble(attributableRiskPositive = 0,
           exposedCasesPositive = 0,
           exposedPositive = 0,
           timePositive = 0,
           countPositive = 0,
           attributableRiskNegative = (exp(meanLogRr) - 1) * exposedCases,
           exposedCasesNegative = exposedCases,
           exposedNegative = exposureOutcomeSetting$nTarget * sumDbMultipliers,
           timeNegative = 1,
           countNegative = 1) %>%
      return()

  } else {
    tibble(attributableRiskPositive = exposedCases - (exposedCases/exp(exposureOutcomeSetting$logRrMean)),
           exposedCasesPositive = exposedCases,
           exposedPositive = exposureOutcomeSetting$nTarget * sumDbMultipliers,
           timePositive = 1,
           countPositive = 1,
           attributableRiskNegative = 0,
           exposedCasesNegative = 0,
           exposedNegative = 0,
           timeNegative = 0,
           countNegative = 0) %>%
      return()
  }
}

computeBaselineAttributableRisk <- function(simulationSettings) {
  attributableRiskPerExposureOutcome <- map_dfr(simulationSettings$exposureOutcomeSettings,
                                                computeAttributableRiskPerExposureOutcome,
                                                simulationSettings = simulationSettings) %>%
    mutate(exposureOutcomeId = row_number())
  # Note: current simulation assumes equal number of exposures at each look, so attributable
  # risk after n looks is just n * (attributable risk in one look):
  attributableRiskPerExposureOutcome <- attributableRiskPerExposureOutcome %>%
    cross_join(tibble(lookId = c(0, seq_len(simulationSettings$looks)))) %>%
    mutate(attributableRiskPositive = attributableRiskPositive * (simulationSettings$looks - .data$lookId),
           exposedCasesPositive = exposedCasesPositive * (simulationSettings$looks - .data$lookId),
           exposedPositive = exposedPositive * (simulationSettings$looks - .data$lookId),
           timePositive = timePositive * (simulationSettings$looks - .data$lookId),
           attributableRiskNegative = attributableRiskNegative * (simulationSettings$looks - .data$lookId),
           exposedCasesNegative = exposedCasesNegative * (simulationSettings$looks - .data$lookId),
           exposedNegative = exposedNegative * (simulationSettings$looks - .data$lookId),
           timeNegative = timeNegative * (simulationSettings$looks - .data$lookId))
  return(attributableRiskPerExposureOutcome)
}

computeAttributableRisk <- function(signals, simulationSettings) {
  attributableRiskPerExposureOutcome <- computeBaselineAttributableRisk(
    simulationSettings = simulationSettings
  )
  attributableRisk <- signals %>%
    # filter(!is.infinite(.data$lookId)) %>%
    left_join(attributableRiskPerExposureOutcome, by = join_by("exposureOutcomeId", "lookId")) %>%
    group_by(.data$label) %>%
    summarise(
      attributableRiskTp = sum(.data$attributableRiskPositive, na.rm = TRUE),
      exposedCasesTp = sum(.data$exposedCasesPositive, na.rm = TRUE),
      exposedTp = sum(.data$exposedPositive, na.rm = TRUE),
      timeTp = sum(.data$timePositive, na.rm = TRUE),
      countTp = sum(.data$countPositive, na.rm = TRUE),
      attributableRiskFp = sum(.data$attributableRiskNegative, na.rm = TRUE),
      exposedCasesFp = sum(.data$exposedCasesNegative, na.rm = TRUE),
      exposedFp = sum(.data$exposedNegative, na.rm = TRUE),
      timeFp = sum(.data$timeNegative, na.rm = TRUE),
      countFp = sum(.data$countNegative, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    cross_join(
      attributableRiskPerExposureOutcome %>%
        filter(.data$lookId == 0) %>%
        summarise(
          attributableRiskPositive = sum(.data$attributableRiskPositive),
          exposedCasesPositive = sum(.data$exposedCasesPositive),
          exposedPositive = sum(.data$exposedPositive),
          timePositive = sum(.data$timePositive),
          countPositive = sum(.data$countPositive),
          attributableRiskNegative = sum(.data$attributableRiskNegative),
          exposedCasesNegative = sum(.data$exposedCasesNegative),
          exposedNegative = sum(.data$exposedNegative),
          timeNegative = sum(.data$timeNegative),
          countNegative = sum(.data$countNegative),
          .groups = "drop"
        )
    ) %>%
    mutate(
      attributableRiskFn = .data$attributableRiskPositive - .data$attributableRiskTp,
      exposedCasesFn = .data$exposedCasesPositive - .data$exposedCasesTp,
      exposedFn = .data$exposedPositive - .data$exposedTp,
      timeFn = .data$timePositive - .data$timeTp,
      countFn = .data$countPositive - .data$countTp,
      attributableRiskTn = .data$attributableRiskNegative - .data$attributableRiskFp,
      exposedCasesTn = .data$exposedCasesNegative - .data$exposedCasesFp,
      exposedTn = .data$exposedNegative - .data$exposedFp,
      timeTn = .data$timeNegative - .data$timeFp,
      countTn = .data$countNegative - .data$countFp
    ) %>%
    select(-"attributableRiskPositive",
           -"exposedCasesPositive",
           -"exposedPositive",
           -"timePositive",
           -"countPositive",
           -"attributableRiskNegative",
           -"exposedCasesNegative",
           -"exposedNegative",
           -"timeNegative",
           -"countNegative")
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
#' @param impactWeighting   Type of impact weighting to use. Can be "none", "time",
#'                          "exposed", "exposed cases", or "attributable cases".
#' @param labels            Which labels to plot? (Values of the `label` column). If `NULL` then all values are plotted.
#' @param alphas            Which alphas to plot? (Values of the `alpha` column). If `NULL` then all alphas are plotted.
#' @param fileName          Optional: the name of the file to save the plot to.
#'
#' @return
#' A GGPlot object.
#'
#' @export
plotDecisionCurves <- function(evaluation,
                               impactWeighting = "time",
                               pickOptimalAlpha = FALSE,
                               showQuartiles = TRUE,
                               labels = NULL,
                               alphas = NULL,
                               fileName = NULL) {
  # temp = evaluation
  # alphas <- sort(unique(evaluation$alpha))[c(1, 3, 7, 10)]
  if (!showQuartiles) {
    evaluation <- evaluation %>%
      filter(.data$iteration == min(.data$iteration))
  }
  if (!is.null(labels)) {
    evaluation <- evaluation %>%
      filter(.data$label %in% labels)
  }
  if (!is.null(alphas)) {
    evaluation <- evaluation %>%
      filter(.data$alpha %in% alphas)
  }
  if (impactWeighting == "none") {
    evaluation <- evaluation %>%
      mutate(
        tp = .data$countTp,
        fp = .data$countFp,
        tn = .data$countTn,
        fn = .data$countFn
      )
  } else if (impactWeighting == "time") {
    evaluation <- evaluation %>%
      mutate(
        tp = .data$timeTp,
        fp = .data$timeFp,
        tn = .data$timeTn,
        fn = .data$timeFn
      )
  } else if (impactWeighting == "exposed") {
    evaluation <- evaluation %>%
      mutate(
        tp = .data$exposedTp,
        fp = .data$exposedFp,
        tn = .data$exposedTn,
        fn = .data$exposedFn
      )
  } else if (impactWeighting == "exposed cases") {
    evaluation <- evaluation %>%
      mutate(
        tp = .data$exposedCasesTp,
        fp = .data$exposedCasesFp,
        tn = .data$exposedCasesTn,
        fn = .data$exposedCasesFn
      )
  } else if (impactWeighting == "attributable cases") {
    evaluation <- evaluation %>%
      mutate(
        tp = .data$attributableRiskTp,
        fp = .data$attributableRiskFp,
        tn = .data$attributableRiskTn,
        fn = .data$attributableRiskFn
      )
  }
  nExposureOutcomes <- evaluation %>%
    mutate(n = .data$tp + .data$fp + .data$tn + .data$fn) %>%
    summarize(n = median(n)) %>%
    pull()
  plotData <- evaluation %>%
    select("label", "tp", "fp", "tn", "fn", "alpha") %>%
    cross_join(tibble(p = seq(0, 1, by = 0.01))) %>%
    mutate(netBenefit = (.data$tp / nExposureOutcomes) - (.data$fp / nExposureOutcomes) * (p / (1-p))) %>%
    mutate(alpha = sprintf("alpha = %0.2f", .data$alpha)) %>%
    group_by(.data$label, .data$p, .data$alpha) %>%
    summarise(
      lb = quantile(.data$netBenefit, 0.25, na.rm = TRUE),
      ub = quantile(.data$netBenefit, 0.75, na.rm = TRUE),
      netBenefit = median(.data$netBenefit),
      .groups = "drop"
    )
  if (pickOptimalAlpha) {
    plotData <- plotData %>%
      group_by(.data$p, .data$label) %>%
      filter(.data$netBenefit == max(.data$netBenefit)) %>%
      ungroup() %>%
      mutate(alpha = 0)
  }

  maxBenefit <- max(plotData$ub, na.rm = TRUE)
  plot <- ggplot2::ggplot(plotData, ggplot2::aes(x = .data$p,
                                                 y = .data$netBenefit,
                                                 ymin = .data$lb,
                                                 ymax = .data$ub,
                                                 group = .data$label,
                                                 color = .data$label,
                                                 fill = .data$label))
  if (showQuartiles) {
    plot <- plot +
      ggplot2::geom_ribbon(color = rgb(0, 0, 0, alpha = 0), alpha = 0.4)
  }
  plot <- plot +
    ggplot2::geom_line(alpha = 0.8, size = 0.5) +
    ggplot2::scale_color_manual(values = c(wesanderson::wes_palette("Darjeeling1", 5), gray(0.5))) +
    ggplot2::scale_fill_manual(values = c(wesanderson::wes_palette("Darjeeling1", 5), gray(0.5))) +
    ggplot2::coord_cartesian(ylim = c(-maxBenefit/5, maxBenefit)) +
    ggplot2::scale_x_continuous("P") +
    ggplot2::scale_y_continuous("Net benefit") +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  if (length(unique(plotData$alpha)) > 1) {
    plot <- plot +
      ggplot2::facet_wrap("alpha")
      # ggplot2::facet_grid(.data$alpha ~ .)

  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(filename = fileName,
                    plot = plot,
                    width = 5,
                    # height = 1 + 2 * length(unique(plotData$alpha)),
                    height = 3,
                    dpi = 200)
  }
  return(plot)
}
