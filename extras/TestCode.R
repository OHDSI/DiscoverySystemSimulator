library(DiscoverySystemSimulator)


maxCores <- 20


simulationSettings <- createSimulationSettings()
simulationsFolder <- "s:/temp/Simulations"
runSimulationIterations(simulationsFolder = simulationsFolder,
                        threads = maxCores)

discoverySystemSettings <- createDiscoverySystemSettings()
discoverySystemSettings$alpha <- 0.50
signalsFolder <- "s:/temp/SignalsAlpha05"
runDiscoverySystemIterations(simulationsFolder = simulationsFolder,
                             signalsFolder = signalsFolder,
                             discoverySystemSettings = discoverySystemSettings,
                             threads = maxCores)

confusionMatrices <- evaluateIterations(signalsFolder = signalsFolder)

confusionMatrices <- evaluateIterations(signalsFolder = signalsFolder, level = "across looks")
mean(confusionMatrices$type1[confusionMatrices$label == "Calibrated MaxSPRT"])
mean(confusionMatrices$type1[confusionMatrices$label == "MaxSPRT"]) * 10000
mean(confusionMatrices$type1[confusionMatrices$label == "P"]) * 10000
alpha <- 0.05
maxSprtAlpha <- alpha /
  length(simulationSettings$exposureOutcomeSettings) /
  length(simulationSettings$databaseSettings) /
  length(simulationSettings$methodSettings) /
  length(simulationSettings$timeAtRiskSettings)
maxSprtAlpha
# confusionMatrices <- readRDS(file.path(signalsFolder, "ConfusionMatrices.rds"))
plotFalsePositiveNegatives(confusionMatrices,
                           strategies = c("Calibrated MaxSPRT"),
                           cumulative = TRUE,
                           alpha = 0.5,
                           fileName = file.path(signalsFolder, "fnfpPlot.png"))

# simulation <- simulateDiscoverySystem(simulationSettings)
# simulation <- readRDS(file.path(simulationsFolder, "Simulation_i1.rds"))
# signals <- runDiscoverySystem(simulation, discoverySystemSettings)
# signals <- readRDS(file.path(signalsFolder, "Signals_i1.rds"))
# computeConfusionMatrix(signals = signals,
#                        simulationSettings = simulationSettings)

discoverySystemSettings <- createDiscoverySystemSettings()
discoverySystemSettings$alpha <- 180
signalsFolder <- "s:/temp/SignalsAlpha180"
runDiscoverySystemIterations(simulationsFolder = simulationsFolder,
                             signalsFolder = signalsFolder,
                             discoverySystemSettings = discoverySystemSettings,
                             threads = maxCores)
