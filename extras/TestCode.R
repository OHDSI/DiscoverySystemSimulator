library(DiscoverySystemSimulator)
maxCores <- 25
cvCacheFile <- "s:/DiscoverSytemSimulations/cvCache.rds"

# Large test ------------------------------------------------------------------
simulationSettings <- createSimulationSettings()
simulationsFolder <- "s:/DiscoverSytemSimulations/Simulations"
runSimulationIterations(simulationsFolder = simulationsFolder,
                        threads = maxCores)

discoverySystemSettings <- createDiscoverySystemSettings()
signalsFolder <- "s:/DiscoverSytemSimulations/SignalsAlpha0.5"
runDiscoverySystemIterations(simulationsFolder = simulationsFolder,
                             signalsFolder = signalsFolder,
                             discoverySystemSettings = discoverySystemSettings,
                             threads = maxCores)

confusionMatrices <- evaluateIterations(signalsFolder = signalsFolder)

confusionMatrices <- evaluateIterations(signalsFolder = signalsFolder, level = "across looks")
mean(confusionMatrices$type1[confusionMatrices$label == "Calibrated MaxSPRT"] ) * 10000
mean(confusionMatrices$type1[confusionMatrices$label == "MaxSPRT"]) * 10000
mean(confusionMatrices$type1[confusionMatrices$label == "P"]) * 10000
alpha <- 0.5
maxSprtAlpha <- alpha /
  length(simulationSettings$exposureOutcomeSettings) /
  length(simulationSettings$databaseSettings) /
  length(simulationSettings$methodSettings) /
  length(simulationSettings$timeAtRiskSettings)
maxSprtAlpha * 10000
# confusionMatrices <- readRDS(file.path(signalsFolder, "ConfusionMatrices.rds"))
plotFalsePositiveNegatives(confusionMatrices,
                           cumulative = TRUE,
                           alpha = 0.5,
                           fileName = file.path(signalsFolder, "fnfpPlot.png"))

# simulation <- simulateDiscoverySystem(simulationSettings)
# simulation <- readRDS(file.path(simulationsFolder, "Simulation_i1.rds"))
signals <- runDiscoverySystem(simulation, discoverySystemSettings)
# signals <- readRDS(file.path(signalsFolder, "Signals_i1.rds"))
# computeConfusionMatrix(signals = signals,
#                        simulationSettings = simulationSettings)

level = "across looks"
simulationSettings <- readRDS(file.path(signalsFolder, "simulationSettings.rds"))
signalFiles <- list.files(path = signalsFolder, pattern = "Signals_.*.rds", full.names = TRUE)

doEvaluation <- function(signalFile) {
  iteration <- as.numeric(gsub("^.*_i", "", gsub(".rds", "", signalFile)))
  signals <- readRDS(signalFile) %>%
    filter(.data$methodId == 1)
  matrix <- computeConfusionMatrix(signals = signals,
                                   simulationSettings = simulationSettings,
                                   level = level) %>%
    mutate(iteration = !!iteration) %>%
    return()
}
confusionMatrices <-map_dfr(signalFiles, doEvaluation)

# Small Simulation test --------------------------------------------------------
simulationsFolder <- "s:/DiscoverSytemSimulations/smallSimulations"
simulationSettings <- createSimulationSettings(
  exposureOutcomeSettings = c(
    lapply(rep(1000, 90), createExposureOutcomeSettings, logRrMean = 0, logRrSd = 0),
    lapply(rep(1000, 10), createExposureOutcomeSettings, logRrMean = log(2), logRrSd = 0)
  ),
  timeAtRiskSettings = list(
    createTimeAtRiskSettings(0, 21)
  ),
  methodSettings = list(
    createMethodSettings(0.10, 0.10)
  ),
  databaseSettings = list(
    createDatabaseSettings(1.0)
  ),
  looks = 10
)

signalsFolder <- "s:/DiscoverSytemSimulations/smallSignals"
discoverySystemSettings <- createDiscoverySystemSettings()

runSimulationIterations(
  simulationsFolder = simulationsFolder,
  threads = maxCores
)
runDiscoverySystemIterations(
  simulationsFolder = simulationsFolder,
  signalsFolder = signalsFolder,
  discoverySystemSettings = discoverySystemSettings,
  threads = maxCores,
  cvCacheFile = cvCacheFile
)
