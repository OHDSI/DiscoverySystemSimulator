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

#' Run multiple iterations of the same simulation.
#'
#' @param simulationsFolder  The folder where the simulation objects should be written.
#' @param simulationSettings An object created using [createSimulationSettings()].
#' @param iterations         The number of iterations.
#' @param threads            The number of threads to use. Ideally, this number should
#'                           be lower than the number of CPU cores.
#'
#' @return
#' This function does not return anything, but writes all output to the specified
#' folder. This folder will contain one file per iteration.
#'
#' @export
runSimulationIterations <- function(simulationsFolder,
                                    simulationSettings = createSimulationSettings(),
                                    iterations = 100,
                                    threads = 4) {
  startTime <- Sys.time()
  if (!dir.exists(simulationsFolder)) {
    message(sprintf("Folder '%s' does not exist, so creating it.", simulationsFolder))
    dir.create(simulationsFolder, recursive = TRUE)
  }

  ParallelLogger::addDefaultFileLogger(file.path(simulationsFolder, "log.txt"))
  ParallelLogger::addDefaultErrorReportLogger(file.path(simulationsFolder, "errorReport.txt"))
  on.exit(ParallelLogger::unregisterLogger("DEFAULT_FILE_LOGGER"))
  on.exit(ParallelLogger::unregisterLogger("DEFAULT_ERRORREPORT_LOGGER"), add = TRUE)

  cluster <- ParallelLogger::makeCluster(threads)
  on.exit(ParallelLogger::stopCluster(cluster), add = TRUE)
  dummy <- ParallelLogger::clusterApply(cluster, 1:iterations, runIteration, simulationSettings = simulationSettings, simulationsFolder = simulationsFolder)

  delta <- Sys.time() - startTime
  message(paste("Running simulations took", signif(delta, 3), attr(delta, "units")))
}

runIteration <- function(iteration, simulationSettings, simulationsFolder) {
  fileName <- file.path(simulationsFolder, sprintf("Simulation_i%d.rds", iteration))
  if (!file.exists(fileName)) {
    set.seed(iteration)
    simulation <- simulateDiscoverySystem(simulationSettings)
    saveRDS(simulation, fileName)
  }
}

#' Run signal discovery on multiple iterations of the same simulation.
#'
#' @details
#' Requires [runSimulationIterations()] to have already been executed.
#'
#' @param simulationsFolder  The folder where the simulation objects were written.
#' @param signalsFolder      The folder where the signal objects will be written.
#' @param discoverySystemSettings An object generated by [createDiscoverySystemSettings()].
#' @param threads            The number of threads to use. Ideally, this number should
#'                           be lower than the number of CPU cores.
#' @param useCaching         Cache intermediate artifacts in the `simulationsFolder`? This will
#'                           greatly speed up subsequent runs.
#'
#' @return
#' This function does not return anything, but writes all output to the specified
#' folder. This folder will contain one file per iteration.
#'
#' @export
runDiscoverySystemIterations <- function(simulationsFolder,
                                         signalsFolder,
                                         discoverySystemSettings = createDiscoverySystemSettings(),
                                         threads = 4,
                                         useCaching = TRUE) {
   startTime <- Sys.time()
  if (!dir.exists(signalsFolder)) {
    message(sprintf("Folder '%s' does not exist, so creating it.", signalsFolder))
    dir.create(signalsFolder, recursive = TRUE)
  }

  ParallelLogger::addDefaultFileLogger(file.path(signalsFolder, "log.txt"))
  ParallelLogger::addDefaultErrorReportLogger(file.path(signalsFolder, "errorReport.txt"))
  on.exit(ParallelLogger::unregisterLogger("DEFAULT_FILE_LOGGER"))
  on.exit(ParallelLogger::unregisterLogger("DEFAULT_ERRORREPORT_LOGGER"), add = TRUE)

  simulationFiles <- list.files(path = simulationsFolder, pattern = "Simulation_.*.rds")

  cluster <- ParallelLogger::makeCluster(threads)
  on.exit(ParallelLogger::stopCluster(cluster), add = TRUE)
  dummy <- ParallelLogger::clusterApply(cluster = cluster,
                                        x = simulationFiles,
                                        fun = runDiscoverySystemIteration,
                                        discoverySystemSettings = discoverySystemSettings,
                                        simulationsFolder = simulationsFolder,
                                        signalsFolder = signalsFolder,
                                        useCaching = useCaching)

  delta <- Sys.time() - startTime
  message(paste("Running discovery system took", signif(delta, 3), attr(delta, "units")))
}

# simulationFile = simulationFiles[1]
runDiscoverySystemIteration <- function(simulationFile,
                                        discoverySystemSettings,
                                        simulationsFolder,
                                        signalsFolder,
                                        useCaching) {
  signalsFile <- gsub("Simulation_", "Signals_", simulationFile)

  if (!file.exists(file.path(signalsFolder, signalsFile))) {
    simulation <- readRDS(file.path(simulationsFolder, simulationFile))

    if (useCaching) {
      cacheFolder <- file.path(simulationsFolder, gsub(".rds", "", gsub("Simulation_", "Cache_", simulationFile)))
    } else {
      cacheFolder <- NULL
    }
    signals <- runDiscoverySystem(simulation = simulation,
                                  discoverySystemSettings = discoverySystemSettings,
                                  cacheFolder = cacheFolder)
    saveRDS(signals, file.path(signalsFolder, signalsFile))
    simulationSettingsFile <- file.path(signalsFolder, "simulationSettings.rds")
    if (!file.exists(simulationSettingsFile)) {
      saveRDS(attr(simulation, "simulationSettings"), simulationSettingsFile)
    }
  }
}

#' Evaluate iterations
#'
#' @details
#' Requires [runSimulationIterations()] to have already been executed.
#'
#' @param signalsFolder      The folder containing the signal objects.
#'
#' @return
#' A tibble combining the confusion matrices of all iterations.
#'
#' @export
evaluateIterations <- function(signalsFolder) {
  simulationSettings <- readRDS(file.path(signalsFolder, "simulationSettings.rds"))
  signalFiles <- list.files(path = signalsFolder, pattern = "Signals_.*.rds", full.names = TRUE)

  doEvaluation <- function(signalFile, simulationSettings) {
    iteration <- as.numeric(gsub("^.*_i", "", gsub(".rds", "", signalFile)))
    signals <- readRDS(signalFile)
    evaluateSignals(
      signals = signals,
      simulationSettings = simulationSettings
    ) %>%
      mutate(iteration = !!iteration) %>%
      return()
  }
  map_dfr(signalFiles, doEvaluation, simulationSettings = simulationSettings) %>%
    return()
}
