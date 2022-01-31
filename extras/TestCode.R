library(DiscoverySystemSimulator)


maxCores <- 25


simulationSettings <- createSimulationSettings()
simulationsFolder <- "s:/DiscoverSytemSimulations/Simulations"
runSimulationIterations(simulationsFolder = simulationsFolder,
                        threads = maxCores)

discoverySystemSettings <- createDiscoverySystemSettings()
discoverySystemSettings$alpha <- 0.50
signalsFolder <- "s:/DiscoverSytemSimulations/SignalsAlpha50"
runDiscoverySystemIterations(simulationsFolder = simulationsFolder,
                             signalsFolder = signalsFolder,
                             discoverySystemSettings = discoverySystemSettings,
                             threads = maxCores)

confusionMatrices <- evaluateIterations(signalsFolder = signalsFolder)

confusionMatrices <- evaluateIterations(signalsFolder = signalsFolder, level = "across looks")
mean(confusionMatrices$type1[confusionMatrices$label == "Calibrated MaxSPRT"] ) * 10000
mean(confusionMatrices$type1[confusionMatrices$label == "MaxSPRT"]) * 10000
mean(confusionMatrices$type1[confusionMatrices$label == "P"]) * 10000
alpha <- 0.5
maxSprtAlpha <- alpha /
  length(simulationSettings$exposureOutcomeSettings) /
  length(simulationSettings$databaseSettings) /
  length(simulationSettings$methodSettings) /
  length(simulationSettings$timeAtRiskSettings)
maxSprtAlpha * 10000
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

level = "across looks"
simulationSettings <- readRDS(file.path(signalsFolder, "simulationSettings.rds"))
signalFiles <- list.files(path = signalsFolder, pattern = "Signals_.*.rds", full.names = TRUE)

doEvaluation <- function(signalFile) {
  iteration <- as.numeric(gsub("^.*_i", "", gsub(".rds", "", signalFile)))
  signals <- readRDS(signalFile) %>%
    filter(.data$methodId == 1)
  matrix <- computeConfusionMatrix(signals = signals,
                                   simulationSettings = simulationSettings,
                                   level = level) %>%
    mutate(iteration = !!iteration) %>%
    return()
}
confusionMatrices <-map_dfr(signalFiles, doEvaluation)



x <- c(1, 0)
y <- c(6, 0)
logTime <- log(c(16000, 32000))

cyclopsData <- Cyclops::createCyclopsData(y ~ x + offset(logTime), modelType = "pr")
fit <- Cyclops::fitCyclopsModel(cyclopsData)
logRr <- coef(fit)["x"]
ci <- confint(fit, parm = "x")
llNull <- Cyclops::getCyclopsProfileLogLikelihood(
  object = fit,
  parm = "x",
  x = 0
)$value
llr <- fit$log_likelihood - llNull

profile <- Cyclops::getCyclopsProfileLogLikelihood(
  object = fit,
  parm = "x",
  bounds = c(log(0.1), 12)
)
plot(profile$point, profile$value)
cumulativeCount <- 6
cumulativeObserved <- 6
p <- 16000 / 32000
llr = dbinom(cumulativeObserved, cumulativeCount, cumulativeObserved / cumulativeCount, log = TRUE) -
  dbinom(cumulativeObserved, cumulativeCount, p, log = TRUE)

