# Copyright 2023 Observational Health Data Sciences and Informatics
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

# group = groups[[361]]
computeCv <- function(group) {
  events <- group %>%
    arrange(.data$lookId) %>%
    mutate(expectedEvents = round(cumsum(expectedEvents))) %>%
    pull(.data$expectedEvents)
  events[2:length(events)] <- events[2:length(events)] - events[1:(length(events) - 1)]
  events <- events[events != 0]
  if (length(events) == 0) {
    cv <- Inf
    cvAlpha <- 0
  } else {
    row <- group %>%
      filter(.data$lookId == max(.data$lookId))
    suppressMessages(
      cv <- EmpiricalCalibration::computeCvBinomial(groupSizes = events,
                                                    z = row$z,
                                                    minimumEvents = 1,
                                                    sampleSize = 1e6,
                                                    alpha = row$alphaPerMethod,
                                                    nullMean = row$systematicErrorMean,
                                                    nullSd = row$systematicErrorSd)
    )
    cvAlpha <- attr(cv, "alpha")
  }
  row <- row %>%
    mutate(cv = cv,
           cvAlpha = cvAlpha) %>%
    select(-"lookId", -"z", -"expectedEvents", -"systematicErrorMean", -"systematicErrorSd")
  return(row)
}
