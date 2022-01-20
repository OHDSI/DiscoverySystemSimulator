library(DiscoverySystemSimulator)

simulationSettings <- createSimulationSettings()
simulation <- simulateDiscoverySystem(simulationSettings)

outputFolder <- "s:/temp/Simulations"
maxCores <- 20

runSimulationIterations(outputFolder = outputFolder,
                        threads = maxCores)
