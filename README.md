DiscoverySystemSimulator
========================


Introduction
============
DiscoverySystemSimulator is an R package that simulates effect discovery in an active surveillance system.

Features
========
Simulates:

- Multiple databases
- Multiple exposure-outcome pairs
- Multiple methods (simulated by resampling the comparator, reflecting different ways to construct the counterfactual), each with different systematic error distributions
- Multiple times at risk
- Multiple looks over time

For each estimate, this is computed:

- Patient time in target and comparator (during time at risk)
- Outcomes observed in target and comparator (during time at risk)
- Estimated effect size, confidence interval, p-value, log likelihood ratio, and likelihood profile (Poisson regression)


Technology
============
DiscoverySystemSimulator is an R package.

System Requirements
===================
Requires R (version 4.0.0 or higher). Installation on Windows requires [RTools](https://cran.r-project.org/bin/windows/Rtools/). Libraries used in DiscoverySystemSimulator require Java.

Installation
=============
1. See the instructions [here](https://ohdsi.github.io/Hades/rSetup.html) for configuring your R environment, including RTools and Java.

2. In R, use the following commands to download and install CohortMethod:

  ```r
  install.packages("remotes")
  remotes::install_github("schuemie/DiscoverySystemSimulator")
  ```

User Documentation
==================
* Package manual: [DiscoverySystemSimulator.pdf](https://raw.githubusercontent.com/schuemie/DiscoverySystemSimulator/master/extrasDiscoverySystemSimulator.pdf)

Support
=======
* Developer questions/comments/feedback: <a href="http://forums.ohdsi.org/c/developers">OHDSI Forum</a>
* We use the <a href="https://github.com/schuemie/DiscoverySystemSimulator/issues">GitHub issue tracker</a> for all bugs/issues/enhancements

Contributing
============
Read [here](https://ohdsi.github.io/Hades/contribute.html) how you can contribute to this package.

License
=======
CohortMethod is licensed under Apache License 2.0

Development
===========
DiscoverySystemSimulator is being developed in R Studio.

### Development status

DiscoverySystemSimulator is under development. Do not use.
