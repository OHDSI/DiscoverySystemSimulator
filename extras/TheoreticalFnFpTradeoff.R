# Create a plot showing the theoretical trade-off between type 1 and type 2 error. Assumes there's only
# one real positive, and the rest is negative.
library(dplyr)
library(purrr)
library(ggplot2)

nEvaluated <- 75   # exposure-outcome pairs to evaluate
power <- 0.8       # 1 - type 2 error
nTarget <- 7e5     # Number of subjects in target
nComparator <- 1e7 # Number of subjects in comparator
alphas <- c(0.05, 0.5, 1, 2, 5) # Expected number of false positives
prevalences <- exp(seq(log(1e-6), log(1e-3), length.out = 100))


computeForAlpha <- function(alpha) {
  correctedAlpha <- alpha / nEvaluated # Bonferroni
  z1MinAlpha <- qnorm(1 - correctedAlpha)
  zBeta <- -qnorm(1 - power)
  pA <- nTarget / (nTarget + nComparator)
  pB <- 1 - pA
  nPrevalences <- length(prevalences)
  mdrr <- rep(0, nPrevalences)
  for (i in 1:nPrevalences) {
    totalEvents <- prevalences[i] * (nTarget + nComparator)
    mdrr[i] <- exp(sqrt((zBeta + z1MinAlpha)^2/(totalEvents * pA * pB)))
  }
  return(tibble(alpha = rep(alpha, nPrevalences),
                prevalence = prevalences,
                mdrr = mdrr))
}

plotData <- map_dfr(alphas, computeForAlpha)
plotData <- plotData %>%
  mutate(alpha = as.factor(alpha))

xBreaks <- c(1e-6, 1e-5, 1e-4, 1e-3)
xLabels <- c("1/1,000,000", "1/100,000", "1/10,000", "1/1,000")

ggplot(plotData, aes(x = prevalence, y = mdrr, group = alpha, color = alpha)) +
  geom_line(size = 1) +
  labs(color = "E(FP)") +
  scale_x_log10("Outcome incidence", breaks = xBreaks, labels = xLabels) +
  scale_y_continuous("Minimum detectable relative risk (MDRR)") +
  coord_cartesian(ylim = c(1, 5))

ggsave(sprintf("s:/temp/fnFpTradeoff%d.png", nEvaluated), width = 5, height = 4, dpi = 200)


