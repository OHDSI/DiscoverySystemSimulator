library(Cyclops)
library(dplyr)

settings <- list(
  nExposed = 100000,
  rr = 2,
  backgroundRate = 0.0001,
  alpha = 0.05,
  minAc = 1,
  sampling = FALSE
)

simulate <- function(settings) {
  if (settings$sampling) {
    exposedCases <- rbinom(1, settings$nExposed, settings$backgroundRate * settings$rr)
    unexposedCases <- rbinom(1, settings$nExposed, settings$backgroundRate)
  } else {
    exposedCases <- settings$nExposed * settings$backgroundRate * settings$rr
    unexposedCases <- settings$nExposed * settings$backgroundRate
  }

  data <- data.frame(outcome = rep(0, settings$nExposed * 2),
                     exposure = c(rep(0, settings$nExposed), rep(1, settings$nExposed)))
  data$outcome[1:unexposedCases] <- 1
  data$outcome[settings$nExposed + 1:exposedCases] <- 1

  cyclopsData <- createCyclopsData(outcome ~ exposure, data = data, modelType = "lr")
  fit <- fitCyclopsModel(cyclopsData)
  logRr <- coef(fit)[2]
  ci <- confint(fit, "exposure")
  seLogRr <- (ci[3] - ci[2])/(2 * qnorm(0.975))
  pRr <- 1 - pnorm(logRr/seLogRr)

  # significantRr <- pRr < settings$alpha

  if (exposedCases == 0) {
    ac <- 0
    acCi <- c(0,0)
    pAc <- 1
  } else {
    ac <- exposedCases - (exposedCases / exp(logRr))
    acCi <- exposedCases - (exposedCases / exp(ci[2:3]))
    pAc <- 1 - pnorm((logRr - log(exposedCases / (exposedCases - settings$minAc)))/seLogRr)
  }
  result <- tibble(rr = exp(logRr),
                   ci95LbRr = exp(ci[2]),
                   ci95UbRr = exp(ci[3]),
                   pRr = pRr,
                   ac = ac,
                   ci95LbAc = acCi[1],
                   ci95UbAc = acCi[2],
                   pAc = pAc)
  result
  return(result)
}

simulate(settings)

# Discovery system simulation ------------------------

simulateSystem <- function(seed) {

  set.seed(seed)
  nExposed <- round(runif(100, 1e3, 1e5))
  backgroundRate <- runif(100, 1e-4, 1e-3)
  rr <- c(rep(1, 90), rep(2, 10))
  minAc <- 1
  alpha <- 0.05

  allSettings <- list()
  for (i in 1:length(nExposed)) {
    allSettings[[i]] <- list(
      nExposed = nExposed[i],
      rr = rr[i],
      backgroundRate = backgroundRate[i],
      minAc = minAc,
      sampling = TRUE
    )
  }
  results <- lapply(allSettings, remoteSimulate)
  results <- bind_rows(results)
  results$groundTruth <- rr > 1
  results <- results %>%
    mutate(signalRr = pRr < alpha,
           signalAc = pAc < alpha) %>%
    mutate(fpSignalRr = signalRr & !groundTruth,
           fnSignalRr = !signalRr & groundTruth,
           fpSignalAc = signalAc & !groundTruth,
           fnSignalAc = !signalAc & groundTruth) %>%
    mutate(fpAcRr = ac * fpSignalRr,
           fnAcRr = ac * fnSignalRr,
           fpAcAc = ac * fpSignalAc,
           fnAcAc = ac * fnSignalAc)

  results %>%
    summarise(fpSignalsRr = sum(fpSignalRr),
              fnSignalsRr = sum(fnSignalRr),
              fpAcsRr = sum(fpAcRr),
              fnAcsRr = sum(fnAcRr),
              fpSignalsAc = sum(fpSignalAc),
              fnSignalsAc = sum(fnSignalAc),
              fpAcsAc = sum(fpAcAc),
              fnAcsAc = sum(fnAcAc)) %>%
    return()
}
cluster <- ParallelLogger::makeCluster(20)
ParallelLogger::clusterRequire(cluster, "dplyr")
ParallelLogger::clusterRequire(cluster, "Cyclops")

# load simulate function in to remove nodes:
f <- function(x, y) {remoteSimulate <<- y}
results <- ParallelLogger::clusterApply(cluster, 1:length(cluster), f, y = simulate)

results <- ParallelLogger::clusterApply(cluster, 1:100, simulateSystem)
ParallelLogger::stopCluster(cluster)

results <- bind_rows(results)
apply(results, 2, mean)
