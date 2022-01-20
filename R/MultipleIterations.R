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

#' Run multiple iterations of the same simulation.
#'
#' @param outputFolder       The folder where the simulation objects should be written.
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
runSimulationIterations <- function(outputFolder,
                                    simulationSettings = createSimulationSettings(),
                                    iterations = 100,
                                    threads = 4) {
  startTime <- Sys.time()
  if (!dir.exists(outputFolder)) {
    message(sprintf("Folder '%s' does not exist, so creating it.", outputFolder))
    dir.create(outputFolder, recursive = TRUE)
  }

  cluster <- ParallelLogger::makeCluster(threads)
  on.exit(ParallelLogger::stopCluster(threads))
  dummy <- ParallelLogger::clusterApply(cluster, 1:iterations, runIteration, simulationSettings = simulationSettings, outputFolder = outputFolder)

  delta <- Sys.time() - startTime
  message(paste("Running simulations took", signif(delta, 3), attr(delta, "units")))
}

runIteration <- function(iteration, simulationSettings, outputFolder) {
  fileName <- file.path(outputFolder, sprintf("Simulation_i%d.rds", iteration))
  if (!file.exists(fileName)) {
    set.seed(iteration)
    simulation <- simulateDiscoverySystem(simulationSettings)
    saveRDS(simulation, fileName)
  }
}
