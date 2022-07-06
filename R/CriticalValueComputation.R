# Copyright 2022 Observational Health Data Sciences and Informatics
#
# This file is part of DiscoverySystemSimulator
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

computeCriticalValuesWithCache <- function(cvsToCompute, cvCacheFile, threads = 1) {
  results <- tibble()
  useCache <- !is.null(cvCacheFile)
  if (useCache) {
    cache <- readCvCache(cvCacheFile)
    if (!is.null(cache)) {
      results <- cvsToCompute %>%
        inner_join(cache, by = colnames(cvsToCompute))
      cvsToCompute <- cvsToCompute %>%
        anti_join(cache, by = colnames(cvsToCompute))
    }
  }

  if (nrow(cvsToCompute) > 0) {
    cluster <- ParallelLogger::makeCluster(threads)
    newResults <- ParallelLogger::clusterApply(
      cluster = cluster,
      x = split(cvsToCompute, seq_len(nrow(cvsToCompute))),
      fun = computeCv)
    ParallelLogger::stopCluster(cluster)

    newResults <- bind_rows(newResults)

    if (useCache) {
      appendToCvCache(newResults, cvCacheFile)
    }
    results <- bind_rows(results, newResults)
  }
  return(results)
}

# row <- cvsToCompute[1, ]
computeCv <- function(row) {
  events <- round(1:row$looks * row$expectedEventsPerLook)
  if (length(events) > 1) {
    events[2:length(events)] <- events[2:length(events)] - events[1:(length(events)-1)]
    events <- events[events != 0]
  }
  if (length(events) == 0) {
    cv <- Inf
    cvAlpha <- 0
  } else {
    suppressMessages(
      cv <- EmpiricalCalibration::computeCvBinomial(groupSizes = events,
                                                    z = row$z,
                                                    minimumEvents = 1,
                                                    sampleSize = 1e6,
                                                    alpha = row$alpha,
                                                    nullMean = row$systematicErrorMean,
                                                    nullSd = row$systematicErrorSd)
    )
    cvAlpha <- attr(cv, "alpha")
  }
  row$cv <- cv
  row$cvAlpha <- cvAlpha
  return(row)
}

readCvCache <- function(cvCacheFile) {
  lockFile <- getLockFileName(cvCacheFile)
  lock <- filelock::lock(lockFile, exclusive = FALSE, timeout = 30000)
  on.exit(filelock::unlock(lock))
  if (file.exists(cvCacheFile)) {
    readRDS(cvCacheFile) %>%
      return()
  } else {
    return(NULL)
  }
}

appendToCvCache <- function(newResults, cvCacheFile) {
  lockFile <- getLockFileName(cvCacheFile)
  lock <- filelock::lock(lockFile, exclusive = TRUE, timeout = 30000)
  on.exit(filelock::unlock(lock))
  if (file.exists(cvCacheFile)) {
    cache <- readRDS(cvCacheFile)

    # Other process could have already added these while we were computing,
    # so make sure to avoid duplicates:
    cacheKeys <- cache %>%
      select(-.data$cv,
             -.data$cvAlpha)
    newResults <- newResults %>%
      anti_join(cacheKeys, by = colnames(cacheKeys))
    if (nrow(newResults) > 0) {
      cache <- bind_rows(cache, newResults)
      saveRDS(cache, cvCacheFile)
    }
  } else {
    saveRDS(newResults, cvCacheFile)
  }
}

getLockFileName <- function(cvCacheFile) {
  return(paste(cvCacheFile, "lock", sep = "."))
}
