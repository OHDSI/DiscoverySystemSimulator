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


#' Create exposure-outcome settings
#'
#' @details
#' For each database, an effect size is sampled from the relative risk distribution. The relative risk
#' indicates the increase in the risk during the risk window.
#'
#' @param nTarget        Number of subjects in the target population.
#' @param nComparator    Number of subjects in the comparator (counterfactual) population.
#' @param backgroundRate Poisson background rate of the outcome.
#' @param logRrMean      The mean of the log relative distribution across databases.
#' @param logRrSd        The standard deviation (SD) of the log relative distribution across databases.
#' @param riskStart      Start of the true risk window (relative to exposure start).
#' @param riskEnd        End of the true risk window (relative to exposure start).
#'
#' @return
#' A settings object
#'
#' @export
createExposureOutcomeSettings <- function(nTarget = 100,
                                          nComparator = nTarget * 2,
                                          backgroundRate = 0.0001,
                                          logRrMean = 0,
                                          logRrSd = 0,
                                          riskStart = 0,
                                          riskEnd = 21) {
  settings <- list()
  for (name in names(formals(createExposureOutcomeSettings))) {
    settings[[name]] <- get(name)
  }
  return(settings)
}

#' Create time-at-risk settings
#'
#' @param start    The start of the time-at-risk (relative to exposure start).
#' @param end      The end of the time-at-risk (relative to exposure start).
#'
#' @return
#' A settings object
#'
#' @export
createTimeAtRiskSettings <- function(start = 0,
                                     end = 21) {
  settings <- list()
  for (name in names(formals(createTimeAtRiskSettings))) {
    settings[[name]] <- get(name)
  }
  return(settings)
}

#' Create method settings
#'
#' @param systematicErrorMean The mean of the systematic error distribution.
#' @param systematicErrorSd  The standard deviation (SD) of the systematic error distribution.
#'
#' @return
#' A settings object
#'
#' @export
createMethodSettings <- function(systematicErrorMean = 0,
                                 systematicErrorSd = 0) {
  settings <- list()
  for (name in names(formals(createMethodSettings))) {
    settings[[name]] <- get(name)
  }
  return(settings)
}

#' Create database settings
#'
#' @param sampleSizeMultiplier The relative sample size of the database.
#'
#' @return
#' A settings object
#'
#' @export
createDatabaseSettings <- function(sampleSizeMultiplier = 1) {
  settings <- list()
  for (name in names(formals(createDatabaseSettings))) {
    settings[[name]] <- get(name)
  }
  return(settings)
}

#' Create simulation settings
#'
#' @param exposureOutcomeSettings A list of objects created using [createExposureOutcomeSettings()].
#' @param timeAtRiskSettings      A list of objects created using [createTimeAtRiskSettings()].
#' @param methodSettings          A list of objects created using [createMethodSettings()].
#' @param databaseSettings        A list of objects created using [createDatabaseSettings()].
#' @param looks                   The number of looks over time.
#'
#' @return
#' A settings object
#'
#' @export
createSimulationSettings <- function(exposureOutcomeSettings = c(
                                       lapply(rep(1000, 90), createExposureOutcomeSettings, logRrMean = 0, logRrSd = 0),
                                       lapply(rep(1000, 10), createExposureOutcomeSettings, logRrMean = log(2), logRrSd = 0.25)
                                     ),
                                     timeAtRiskSettings = list(
                                       createTimeAtRiskSettings(0, 7),
                                       createTimeAtRiskSettings(0, 21),
                                       createTimeAtRiskSettings(0, 42),
                                       createTimeAtRiskSettings(0, 90)
                                     ),
                                     methodSettings = list(
                                       createMethodSettings(0.05, 0.05),
                                       createMethodSettings(0.10, 0.10),
                                       createMethodSettings(0.20, 0.20)
                                     ),
                                     databaseSettings = list(
                                       createDatabaseSettings(1.0),
                                       createDatabaseSettings(0.5),
                                       createDatabaseSettings(2.0)
                                     ),
                                     looks = 10) {
  settings <- list()
  for (name in names(formals(createSimulationSettings))) {
    settings[[name]] <- get(name)
  }
  return(settings)
}
