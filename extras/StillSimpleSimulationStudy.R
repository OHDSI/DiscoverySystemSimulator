# A simple simulation study exploring the effect of eliminating a method or a database
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

discoverySystemSettings <- createDiscoverySystemSettings(alpha = 1)
discoverySystemSettings$signalsFolder <- "SignalsAlpha1"
discoverySystemSettings$label <- "Baseline"
discoverySystemSettingsList[[length(discoverySystemSettingsList) + 1]] <- discoverySystemSettings

discoverySystemSettings <- createDiscoverySystemSettings(alpha = 1,
                                                         databaseIdsToIgnore = 2)
discoverySystemSettings$signalsFolder <- "DropDb2"
discoverySystemSettings$label <- "Drop smallest DB"
discoverySystemSettingsList[[length(discoverySystemSettingsList) + 1]] <- discoverySystemSettings

discoverySystemSettings <- createDiscoverySystemSettings(alpha = 1,
                                                         databaseIdsToIgnore = 3)
discoverySystemSettings$signalsFolder <- "DropDb3"
discoverySystemSettings$label <- "Drop largest DB"
discoverySystemSettingsList[[length(discoverySystemSettingsList) + 1]] <- discoverySystemSettings

discoverySystemSettings <- createDiscoverySystemSettings(alpha = 1,
                                                         methodIdsToIgnore = 1)
discoverySystemSettings$signalsFolder <- "DropMethod1"
discoverySystemSettings$label <- "Drop best method"
discoverySystemSettingsList[[length(discoverySystemSettingsList) + 1]] <- discoverySystemSettings

discoverySystemSettings <- createDiscoverySystemSettings(alpha = 1,
                                                         methodIdsToIgnore = 3)
discoverySystemSettings$signalsFolder <- "DropMethod3"
discoverySystemSettings$label <- "Drop worst method"
discoverySystemSettingsList[[length(discoverySystemSettingsList) + 1]] <- discoverySystemSettings

# Run discovery system simulations --------------------------------------
for (discoverySystemSettings in discoverySystemSettingsList) {
  message(sprintf("Running '%s', writing to '%s'.", discoverySystemSettings$label, discoverySystemSettings$signalsFolder))
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
saveRDS(confusionMatrices, file.path(folder, "confusionMatrices2.rds"))

# Plot false negatives and false positives ------------------------------
confusionMatrices <- readRDS(file.path(folder, "confusionMatrices2.rds"))
plotFalsePositiveNegatives(confusionMatrices,
                           cumulative = TRUE,
                           alpha = NA,
                           fileName = file.path(folder, "fnfpPlot2.png"))



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
