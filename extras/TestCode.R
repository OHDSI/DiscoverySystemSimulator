library(DiscoverySystemSimulator)
maxCores <- 25
mainFolder <- "d:/DiscoverSytemSimulations_old"
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
evaluation <- readRDS(file.path(signalsFolder, "evaluation.rds"))
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
  impactWeighting = "none",
  fileName = file.path(signalsFolder, "decisionCurves.png")
)
plotDecisionCurves(
  evaluation = evaluation,
  impactWeighting = "time",
  fileName = file.path(signalsFolder, "decisionCurvesTime.png")
)
plotDecisionCurves(
  evaluation = evaluation,
  impactWeighting = "exposed",
  fileName = file.path(signalsFolder, "decisionCurvesExposed.png")
)
plotDecisionCurves(
  evaluation = evaluation,
  impactWeighting = "exposed cases",
  fileName = file.path(signalsFolder, "decisionCurveExposedCases.png")
)
plotDecisionCurves(
  evaluation = evaluation,
  impactWeighting = "attributable cases",
  fileName = file.path(signalsFolder, "decisionCurvesAttributable cases.png")
)
plotDecisionCurves(
  evaluation = evaluation,
  impactWeighting = "none",
  pickOptimalAlpha = TRUE,
  fileName = file.path(signalsFolder, "decisionCurvesOpt.png")
)
plotDecisionCurves(
  evaluation = evaluation,
  impactWeighting = "time",
  pickOptimalAlpha = TRUE,
  fileName = file.path(signalsFolder, "decisionCurvesTimeOpt.png")
)

alpha <- unique(evaluation$alpha)[6]
plotDecisionCurves(evaluation = evaluation,
                   impactWeighting = "none",
                   pickOptimalAlpha = FALSE,
                   showQuartiles = FALSE,
                   alphas = min(evaluation$alpha),
                   fileName = file.path(signalsFolder, "decisionCurve1.png"))
plotDecisionCurves(evaluation = evaluation,
                   impactWeighting = "none",
                   pickOptimalAlpha = FALSE,
                   showQuartiles = TRUE,
                   alphas = min(evaluation$alpha),
                   fileName = file.path(signalsFolder, "decisionCurve2.png"))
plotDecisionCurves(evaluation = evaluation,
                   impactWeighting = "none",
                   pickOptimalAlpha = FALSE,
                   showQuartiles = TRUE,
                   alphas = unique(evaluation$alpha)[c(1, 3, 7, 10)],
                   fileName = file.path(signalsFolder, "decisionCurve3.png"))
plotDecisionCurves(evaluation = evaluation,
                   impactWeighting = "none",
                   pickOptimalAlpha = TRUE,
                   showQuartiles = TRUE,
                   fileName = file.path(signalsFolder, "decisionCurve4.png"))
plotDecisionCurves(evaluation = evaluation,
                   impactWeighting = "time",
                   pickOptimalAlpha = TRUE,
                   showQuartiles = TRUE,
                   fileName = file.path(signalsFolder, "decisionCurve5.png"))
plotDecisionCurves(evaluation = evaluation,
                   impactWeighting = "exposed",
                   pickOptimalAlpha = TRUE,
                   showQuartiles = TRUE,
                   fileName = file.path(signalsFolder, "decisionCurve6.png"))
plotDecisionCurves(evaluation = evaluation,
                   impactWeighting = "exposed cases",
                   pickOptimalAlpha = TRUE,
                   showQuartiles = TRUE,
                   fileName = file.path(signalsFolder, "decisionCurve7.png"))

