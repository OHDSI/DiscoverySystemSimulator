% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunDiscovery.R
\name{createDiscoverySystemSettings}
\alias{createDiscoverySystemSettings}
\title{Create setting for the discovery system}
\usage{
createDiscoverySystemSettings(
  alpha = exp(seq(log(0.05), log(5), length.out = 10)),
  databaseIdsToIgnore = c(),
  methodIdsToIgnore = c()
)
}
\arguments{
\item{alpha}{The family-wise type 1 error rate we're aiming for. Can be a grid of values.}

\item{databaseIdsToIgnore}{A list of database IDs to ignore (first database has ID 1, second 2, etc.).}

\item{methodIdsToIgnore}{A list of method IDs to ignore (first method has ID 1, method 2, etc.).}
}
\value{
A settings object
}
\description{
Create setting for the discovery system
}
