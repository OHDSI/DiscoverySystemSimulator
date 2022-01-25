library(DiscoverySystemSimulator)


maxCores <- 20


simulationSettings <- createSimulationSettings()
simulationsFolder <- "s:/temp/Simulations"
runSimulationIterations(simulationsFolder = simulationsFolder,
                        threads = maxCores)

discoverySystemSettings <- createDiscoverySystemSettings()
discoverySystemSettings$alpha <- 0.50
signalsFolder <- "s:/temp/SignalsAlpha50"
runDiscoverySystemIterations(simulationsFolder = simulationsFolder,
                             signalsFolder = signalsFolder,
                             discoverySystemSettings = discoverySystemSettings,
                             threads = maxCores)

evaluateIterations(signalsFolder = signalsFolder)
confusionMatrices <- readRDS(file.path(signalsFolder, "ConfusionMatrices.rds"))
plotFalsePositiveNegatives(confusionMatrices,
                           strategies = c("Calibrated MaxSPRT"),
                           cumulative = TRUE,
                           alpha = 0.5,
                           fileName = file.path(signalsFolder, "fnfpPlot.png"))

# simulation <- simulateDiscoverySystem(simulationSettings)
# simulation <- readRDS(file.path(outputFolder, "Simulation_i1.rds"))
# signals <- runDiscoverySystem(simulation, discoverySystemSettings)
# signals <- readRDS(file.path(signalsFolder, "Signals_i1.rds"))
# computeConfusionMatrix(signals = signals,
#                        simulationSettings = simulationSettings)
