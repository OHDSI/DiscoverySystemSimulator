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

# group = groups[[1]]
computeCv <- function(group) {
  events <- group %>%
    tidyr::pivot_longer(names_to = "lookId", names_prefix = "look_", cols = starts_with("look_")) %>%
    arrange(.data$lookId) %>%
    mutate(expectedEvents = round(cumsum(.data$value))) %>%
    pull(.data$expectedEvents)

  if (length(events) > 1) {
    events[2:length(events)] <- events[2:length(events)] - events[1:(length(events) - 1)]
  }
  events <- events[events != 0]
  if (length(events) == 0) {
    cv <- Inf
    cvAlpha <- 0
  } else {
    suppressMessages(
      cv <- EmpiricalCalibration::computeCvBinomial(groupSizes = events,
                                                    z = group$z,
                                                    minimumEvents = 1,
                                                    sampleSize = 1e6,
                                                    alpha = group$alphaPerDatabase,
                                                    nullMean = group$systematicErrorMean,
                                                    nullSd = group$systematicErrorSd)
    )
    cvAlpha <- attr(cv, "alpha")
  }
  group <- group %>%
    mutate(cv = cv,
           cvAlpha = cvAlpha)
  return(group)
}
