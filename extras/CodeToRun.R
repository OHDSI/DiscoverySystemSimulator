library(DiscoverySystemSimulator)
library(dplyr)
maxCores <- 25
mainFolder <- "d:/DiscoverSytemSimulations"
simulationsFolder <- file.path(mainFolder, "Simulations")
signalsFolder <- file.path(mainFolder, "SignalsBaseline")

# Create settings --------------------------------------------------------------
set.seed(1)
createEos <- function(i) {
  if (i %% 2 == 0) {
    # Established product
    nTargetFirstLook <- round(exp(runif(1, log(100), log(1000))))
    nTargetNextLooks <- round(exp(runif(1, log(100), log(1000))))
    nTargetNextLooksDelta <- max(round(nTargetNextLooks * rnorm(1, sd = 0.5)),
                                 -nTargetNextLooks / 9)
  } else {
    # New to market
    nTargetFirstLook <- round(exp(runif(1, log(1), log(100))))
    nTargetNextLooks <- round(exp(runif(1, log(100), log(1000))))
    nTargetNextLooksDelta <- abs(round(nTargetNextLooks * rnorm(1, sd = 0.5)))
  }
  backgroundRate <- exp(runif(1, log(0.0001), log(0.001)))
  if (i <= 90) {
    logRrMean = 0
    logRrSd = 0
  } else {
    logRrMean = runif(1, log(1.25), log(4))
    logRrSd = 0.2
  }
  riskEnd <- round(runif(1, 7, 90))
  settings <- createExposureOutcomeSettings(nTargetFirstLook = nTargetFirstLook,
                                            nTargetNextLooks = nTargetNextLooks,
                                            nTargetNextLooksDelta = nTargetNextLooksDelta,
                                            nComparatorMultiplier = 2,
                                            backgroundRate = backgroundRate,
                                            logRrMean = logRrMean,
                                            logRrSd = logRrSd,
                                            riskStart = 0,
                                            riskEnd = riskEnd)
  # print(format(unlist(settings), scientific = FALSE))
  return(settings)
}
eosList <- lapply(1:100, createEos)
simulationSettings <- createSimulationSettings(exposureOutcomeSettings = eosList)

# Run simulation ---------------------------------------------------------------
runSimulationIterations(
  simulationSettings = simulationSettings,
  simulationsFolder = simulationsFolder,
  threads = maxCores)

# Run baseline detection system ------------------------------------------------
discoverySystemSettings <- createDiscoverySystemSettings()
runDiscoverySystemIterations(simulationsFolder = simulationsFolder,
                             signalsFolder = signalsFolder,
                             discoverySystemSettings = discoverySystemSettings,
                             threads = maxCores)
