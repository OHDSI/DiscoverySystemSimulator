library(DiscoverySystemSimulator)
maxCores <- 25
mainFolder <- "d:/DiscoverSytemSimulations"
cvCacheFile <- file.path(mainFolder, "cvCache.rds")

# Large test ------------------------------------------------------------------
simulationSettings <- createSimulationSettings()
simulationsFolder <- file.path(mainFolder, "Simulations")
runSimulationIterations(simulationsFolder = simulationsFolder,
                        threads = maxCores)

discoverySystemSettings <- createDiscoverySystemSettings()
signalsFolder <- file.path(mainFolder, "SignalsBaseline")
runDiscoverySystemIterations(simulationsFolder = simulationsFolder,
                             signalsFolder = signalsFolder,
                             discoverySystemSettings = discoverySystemSettings,
                             threads = maxCores,
                             cvCacheFile = cvCacheFile)

# Small Simulation test --------------------------------------------------------
simulationsFolder <- file.path(mainFolder, "smallSimulations")
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

signalsFolder <- file.path(mainFolder, "smallSignals")
discoverySystemSettings <- createDiscoverySystemSettings()

runSimulationIterations(
  simulationsFolder = simulationsFolder,
  simulationSettings = simulationSettings,
  threads = maxCores
)
runDiscoverySystemIterations(
  simulationsFolder = simulationsFolder,
  signalsFolder = signalsFolder,
  discoverySystemSettings = discoverySystemSettings,
  threads = maxCores,
  cvCacheFile = cvCacheFile
)

# Do evaluations ---------------------------------------------------------------
evaluation <- evaluateIterations(signalsFolder)
saveRDS(evaluation, file.path(signalsFolder, "evaluation.rds"))
plotFalsePositiveNegatives(
  evaluation = evaluation,
  labels = "Calibrated MaxSPRT",
  alphas = NULL,
  fileName = file.path(signalsFolder, "fpfn_CalibratedMaxSprt.png")
)

plotFalsePositiveNegatives(
  evaluation = evaluation,
  labels = "Calibrated MaxSPRT",
  alphas = unique(evaluation$alpha)[c(1, 6, 10)],
  fileName = file.path(signalsFolder, "fpfn_CalibratedMaxSprt3Alphas.png")
)

plotFalseDiscoveryRate(
  evaluation = evaluation,
  labels = "Calibrated MaxSPRT",
  alphas = NULL,
  fileName = file.path(signalsFolder, "fdr_CalibratedMaxSprt.png")
)
plotFalseDiscoveryRate(
  evaluation = evaluation,
  labels = "Calibrated MaxSPRT",
  alphas = unique(evaluation$alpha)[c(1, 6, 10)],
  fileName = file.path(signalsFolder, "fdr_CalibratedMaxSprt3Alphas.png")
)

plotRoc(
  evaluation = evaluation,
  labels = "Calibrated MaxSPRT",
  alphas = NULL,
  fileName = file.path(signalsFolder, "roc_CalibratedMaxSprt.png")
)

plotDecisionCurves(
  evaluation = evaluation,
  fileName = file.path(signalsFolder, "decisionCurves_CalibratedMaxSprt.png")
)


