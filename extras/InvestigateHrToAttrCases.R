# Some code for evaluating converting HR to attributable cases using dumb approach
library(survival)
library(Cyclops)
library(dplyr)

settings <- list(
  n = 10000,
  treatedFraction = 0.2,
  nStrata = 10,
  minBackgroundHazard = 2e-6,
  maxBackgroundHazard = 2e-5,
  hazardRatio = 2,
  censorHazard = 0.01
)

runSimulation <- function(i, settings) {
  set.seed(i)
  population <- data.frame(rowId = 1:settings$n,
                           stratumId = round(runif(settings$n,
                                                   min = 1,
                                                   max = settings$nStrata)),
                           y = 0,
                           x = as.numeric(runif(settings$n) < settings$treatedFraction))
  population$timeToCensor <- 1 + round(rexp(n = settings$n, settings$censorHazard))
  strataBackgroundHazard <- runif(settings$nStrata,
                                  min = settings$minBackgroundHazard,
                                  max = settings$maxBackgroundHazard)
  population$backgroundHazard <- strataBackgroundHazard[population$stratumId]
  population$timeToBackgroundEvent <- 1 + round(rexp(n = settings$n, population$backgroundHazard))
  population$timeToCausalEvent <- Inf
  idx <- population$x == 1
  population$timeToCausalEvent[idx] <- 1 + round(rexp(sum(idx), population$backgroundHazard[idx] * (settings$hazardRatio - 1)))
  population$time <- pmin(population$timeToCensor, population$timeToBackgroundEvent, population$timeToCausalEvent)
  population$y <- population$time != population$timeToCensor
  population$causalY <- population$y & population$time == population$timeToCausalEvent
  exposedCases <- sum(population$y[idx])

  cyclopsData <- createCyclopsData(Surv(time, y) ~ x + strata(stratumId), data = population, modelType = "cox")
  fit <- fitCyclopsModel(cyclopsData)
  ci <- exp(confint(fit, 1))[2:3]
  # effectiveHr <- sum(population$y) / sum(population$y & ! population$causalY)
  effectiveHr <- settings$hazardRatio

  # effectiveAc <- sum(population$causalY)
  effectiveAc <- exposedCases - (exposedCases/settings$hazardRatio)
  acCi <- exposedCases - (exposedCases / ci)
  return(data.frame(
    exposedCases = exposedCases,
    hrCoverage = ci[1] < effectiveHr & effectiveHr < ci[2],
    acCoverage = acCi[1] < effectiveAc & effectiveAc < acCi[2])
  )
}
cluster <- ParallelLogger::makeCluster(20)
ParallelLogger::clusterRequire(cluster, "survival")
ParallelLogger::clusterRequire(cluster, "Cyclops")

results <- ParallelLogger::clusterApply(cluster, 1:1000, runSimulation, settings = settings)
bind_rows(results) %>%
  summarise(hrCoverage = mean(hrCoverage, na.rm = TRUE),
            acCoverage = mean(acCoverage, na.rm = TRUE),
            zeroExposedCases = mean(exposedCases == 0))

ParallelLogger::stopCluster(cluster)

