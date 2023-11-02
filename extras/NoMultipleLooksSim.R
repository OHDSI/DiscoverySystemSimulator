# A simple simulation with 1 look and 1 TAR only

library(DiscoverySystemSimulator)

maxCores <- 4
folder <- "e:/DiscoverSytemSimulations"
simulationsFolder <- file.path(folder, "Simulations")
signalsFolder <-  file.path(folder, "Signals")

# Simulate estimates ----------------------------------------------------
simulationSettings <- createSimulationSettings(
  exposureOutcomeSettings = c(
    lapply(rep(100, 900), createExposureOutcomeSettings, logRrMean = 0, logRrSd = 0),
    lapply(rep(100, 100), createExposureOutcomeSettings, logRrMean = log(2), logRrSd = 0.25)
  ),
  timeAtRiskSettings = list(
    createTimeAtRiskSettings(0, 21)
  ),
  methodSettings = list(
    createMethodSettings(0.05, 0.05),
    createMethodSettings(0.10, 0.10)
  ),
  databaseSettings = list(
    createDatabaseSettings(1.0),
    createDatabaseSettings(0.5),
    createDatabaseSettings(2.0),
    createDatabaseSettings(1.0),
    createDatabaseSettings(0.5),
    createDatabaseSettings(2.0),
    createDatabaseSettings(1.0),
    createDatabaseSettings(0.5),
    createDatabaseSettings(2.0)
  ),
  looks = 1)
runSimulationIterations(simulationsFolder = simulationsFolder,
                        simulationSettings = simulationSettings,
                        threads = maxCores,
                        iterations = 1)

# Create discovery system settings --------------------------------------
discoverySystemSettings <- createDiscoverySystemSettings(alpha = 0.5)

# Run discovery system simulations --------------------------------------
runDiscoverySystemIterations(simulationsFolder = simulationsFolder,
                             signalsFolder = signalsFolder,
                             discoverySystemSettings = discoverySystemSettings,
                             threads = maxCores,
                             useCaching = TRUE)

# Compute confusion matrices --------------------------------------------
confusionMatrices <- evaluateIterations(signalsFolder = signalsFolder)

saveRDS(confusionMatrices, file.path(folder, "confusionMatrices.rds"))
