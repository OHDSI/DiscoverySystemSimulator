# Some code to evaluate whether we can formulate the critical value as a p-value
library(dplyr)
library(ggplot2)
simulationsFolder <- "s:/DiscoverSytemSimulations/Simulations"
simulation <- readRDS(file.path(simulationsFolder, "Simulation_i1.rds"))

simulationSettings <- attr(simulation, "simulationSettings")
negativeControlIds <- c()
for (i in 1:length(simulationSettings$exposureOutcomeSettings)) {
  exposureOutcomeSetting <- simulationSettings$exposureOutcomeSettings[[i]]
  if (exposureOutcomeSetting$logRrMean == 0 &&
      exposureOutcomeSetting$logRrSd == 0) {
    negativeControlIds <- c(negativeControlIds, i)
  }
}

profiles <- attr(simulation, "profiles")
attr(simulation, "profiles") <- NULL
subsets <- simulation %>%
  filter(methodId == 2 & databaseId == 1 & timeAtRiskId == 1)
subsets <- split(subsets, subsets$lookId)

for (i in 1:length(subsets)) {
  subset <- subsets[[i]]
  ncs <- subset %>%
    filter(exposureOutcomeId %in% negativeControlIds)

  ncProfiles <- profiles[ncs$profileIdx]
  ncProfiles <- ncProfiles[!sapply(ncProfiles, is.null)]
  null <- EmpiricalCalibration::fitNullNonNormalLl(ncProfiles)

  subset$calibratedP <- EmpiricalCalibration::calibrateP(null, subset$logRr, subset$seLogRr, twoSided = FALSE)
  subset <- subset %>%
    mutate(calibratedP = tidyr::replace_na(calibratedP, 1))

  subset$calibratedPfromLlr <- 1
  subsetProfiles <- profiles[subset$profileIdx]
  idx <- !sapply(subsetProfiles, is.null)
  subset$calibratedPfromLlr[idx] <- EmpiricalCalibration:::computePFromLlr(EmpiricalCalibration::calibrateLlr(null, subsetProfiles[idx]), subset$logRr[idx])


  subset <- subset %>%
    mutate(fractionP = 1 / n()) %>%
    arrange(calibratedP) %>%
    mutate(fractionP = cumsum(fractionP))
  ggplot(subset, aes(x = calibratedP, y =  fractionP)) +
    geom_abline(slope = 1) +
    geom_line() +
    xlim(0, 0.5) +
    ylim(0, 0.5)

  subset <- subset %>%
    mutate(fractionPFromLlr = 1 / n()) %>%
    arrange(calibratedPfromLlr) %>%
    mutate(fractionPFromLlr = cumsum(fractionPFromLlr))
  ggplot(subset, aes(x = calibratedPfromLlr, y =  fractionPFromLlr)) +
    geom_abline(slope = 1) +
    geom_line() +
    xlim(0, 0.5) +
    ylim(0, 0.5)


  # plot(subsetProfiles[[100]]$point, subsetProfiles[[100]]$value)
  # subset[100, ]

}
