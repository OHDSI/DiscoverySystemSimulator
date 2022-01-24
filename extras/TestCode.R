library(DiscoverySystemSimulator)

simulationSettings <- createSimulationSettings()
simulation <- simulateDiscoverySystem(simulationSettings)

outputFolder <- "s:/temp/Simulations"
maxCores <- 20

ParallelLogger::addDefaultFileLogger(file.path(outputFolder, "discoverySystemLog.txt"))
ParallelLogger::addDefaultErrorReportLogger(file.path(outputFolder, "errorReport.txt"))

runSimulationIterations(outputFolder = outputFolder,
                        threads = maxCores)


discoverySystemSettings <- createDiscoverySystemSettings()
runDiscoverySystemIterations(outputFolder = outputFolder,
                             discoverySystemSettings = discoverySystemSettings,
                             threads = maxCores)



simulation <- readRDS(file.path(outputFolder, "Simulation_i1.rds"))
discoverySystemSettings <- createDiscoverySystemSettings()

signals <- runDiscoverySystem(simulation, discoverySystemSettings)

ParallelLogger::unregisterLogger("DEFAULT_FILE_LOGGER")
ParallelLogger::unregisterLogger("DEFAULT_ERRORREPORT_LOGGER")
saveRDS(signals, file.path(outputFolder, "Signals_i1.rds"))



signals <- readRDS(file.path(outputFolder, "Signals_i1.rds"))
computeConfusionMatrix(signals = signals,
                       simulationSettings = simulationSettings)
