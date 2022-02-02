# A simple simulation study exploring the effect of the alpha parameter
library(DiscoverySystemSimulator)

maxCores <- 25
folder <- "s:/DiscoverSytemSimulations"
simulationsFolder <- file.path(folder, "Simulations")


# Simulate estimates ----------------------------------------------------
simulationSettings <- createSimulationSettings()
runSimulationIterations(simulationsFolder = simulationsFolder,
                        threads = maxCores)

# Create discovery system settings --------------------------------------
discoverySystemSettingsList <- list()
for (alpha in c(0.05, 0.5, 1, 2, 5, 10)) {
  discoverySystemSettings <- createDiscoverySystemSettings()
  discoverySystemSettings$alpha <- alpha
  discoverySystemSettings$signalsFolder <- sprintf("SignalsAlpha%s", alpha)
  discoverySystemSettings$label <- sprintf("alpha = %s", alpha)
  discoverySystemSettingsList[[length(discoverySystemSettingsList) + 1]] <- discoverySystemSettings
}

# Run discovery system simulations --------------------------------------
for (discoverySystemSettings in discoverySystemSettingsList) {
  message(sprintf("Running %s", discoverySystemSettings$signalsFolder))
  runDiscoverySystemIterations(simulationsFolder = simulationsFolder,
                               signalsFolder = file.path(folder, discoverySystemSettings$signalsFolder),
                               discoverySystemSettings = discoverySystemSettings,
                               threads = maxCores)
}

# Compute confusion matrices --------------------------------------------
confusionMatricesList <- list()
for (discoverySystemSettings in discoverySystemSettingsList) {
  confusionMatrices <- evaluateIterations(signalsFolder = file.path(folder, discoverySystemSettings$signalsFolder))
  confusionMatrices <- confusionMatrices[confusionMatrices$label == "Calibrated MaxSPRT", ]
  confusionMatrices$label <- discoverySystemSettings$label
  confusionMatricesList[[length(confusionMatricesList) + 1]] <- confusionMatrices
}
confusionMatrices <- do.call(rbind, confusionMatricesList)
saveRDS(confusionMatrices, file.path(folder, "confusionMatrices.rds"))

# Plot false negatives and false positives ------------------------------
confusionMatrices <- readRDS(file.path(folder, "confusionMatrices.rds"))
plotFalsePositiveNegatives(confusionMatrices,
                           cumulative = TRUE,
                           alpha = NA,
                           fileName = file.path(folder, "fnfpPlot.png"))



# confusionMatrices <- evaluateIterations(signalsFolder = signalsFolder, level = "across looks")
# mean(confusionMatrices$type1[confusionMatrices$label == "Calibrated MaxSPRT"]) * 10000
# mean(confusionMatrices$type1[confusionMatrices$label == "MaxSPRT"]) * 10000
# mean(confusionMatrices$type1[confusionMatrices$label == "P"]) * 10000
# alpha <- 0.5
# maxSprtAlpha <- alpha /
#   length(simulationSettings$exposureOutcomeSettings) /
#   length(simulationSettings$databaseSettings) /
#   length(simulationSettings$methodSettings) /
#   length(simulationSettings$timeAtRiskSettings)
# maxSprtAlpha * 10000
# confusionMatrices <- readRDS(file.path(signalsFolder, "ConfusionMatrices.rds"))
# plotFalsePositiveNegatives(confusionMatrices,
#                            strategies = c("Calibrated MaxSPRT"),
#                            cumulative = TRUE,
#                            alpha = NA,
#                            fileName = file.path(signalsFolder, "fnfpPlot.png"))
# simulation <- simulateDiscoverySystem(simulationSettings)
# simulation <- readRDS(file.path(simulationsFolder, "Simulation_i1.rds"))
# signals <- runDiscoverySystem(simulation, discoverySystemSettings)
# signals <- readRDS(file.path(signalsFolder, "Signals_i1.rds"))
# computeConfusionMatrix(signals = signals,
#                        simulationSettings = simulationSettings)
