# A simple simulation with 1 look and 1 TAR only

library(DiscoverySystemSimulator)

maxCores <- 16
folder <- "d:/DiscoverSytemSimulations"
simulationsFolder <- file.path(folder, "Simulations")
signalsFolder <-  file.path(folder, "Signals")

# Simulate estimates ----------------------------------------------------
simulationSettings <- createSimulationSettings(
  exposureOutcomeSettings = c(
    lapply(rep(1000, 900), createExposureOutcomeSettings, logRrMean = 0, logRrSd = 0),
    lapply(rep(1000, 100), createExposureOutcomeSettings, logRrMean = log(2), logRrSd = 0.25)
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
                        iterations = 100)

# Create discovery system settings --------------------------------------
discoverySystemSettings <- createDiscoverySystemSettings(
  alpha = c(0.05, 0.25, 0.5),
  useCalibratedMaxSprt = FALSE,
  useMaxSprt = FALSE)

# Run discovery system simulations --------------------------------------
runDiscoverySystemIterations(simulationsFolder = simulationsFolder,
                             signalsFolder = signalsFolder,
                             discoverySystemSettings = discoverySystemSettings,
                             threads = maxCores,
                             useCaching = TRUE)

# Compute confusion matrices --------------------------------------------
confusionMatrices <- evaluateIterations(signalsFolder = signalsFolder)

saveRDS(confusionMatrices, file.path(folder, "confusionMatrices.rds"))

library(dplyr)
x <- confusionMatrices %>%
  group_by(label, alpha) %>%
  summarize(mean(type1), mean(type2))
