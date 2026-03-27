# Disparity analyses are robust to ancestral state estimation uncertainty

Author(s): [Caleb Scutt](mailto:cnscutt1.@sheffield.ac.uk), [Natalie Cooper](https://github.com/nhcooper123), [Gavin Thomas](https://github.com/ghthomas) and [Thomas Guillerme](https://github.com/TGuillerme)

This study uses simulations to evaluate different methods of ancestral state estimation for recovering morphological disparity under varying evolutionary models, fossil sampling densities, and trait types.
This repository contains the scripts to reproduce the analyses of the paper.

---

## Data

The `Data/` folder contains the trees and simulated matrices that are used in the analyses. These are also easily reproduced using the `Analyses/` folder.
The ancestral state estimations, ordinations, disparity and LMM outputs are on dryad (add link here when).

---

## Analyses

All code used to run the analyses are contained in this repo. 
The `Analyses/` folder contains the scripts to reproduce the analyses, functionalised so that you can reproduce the analyses with however many replicates you desire.
The `Cluster/` folder contains scripts used to produce the analyses on the Sheffield HPC.
Note, that the analyses take over 1.5 CPU years to run.

* **01_tree_generation.R** generates the 300 trees (100 for each tree size level) used in both continuous and discrete analyses.
* **02_continuous_sims_ace.R** & **03_discrete_sims_ace.R** simulates continuous and discrete traits across the trees, simulates fossil sampling, runs pre/post-ordination and point/distribution ancestral state estimation and generates trait spaces.
* **05_continuous_disparity.R** & **06_discrete_disparity.R** calculates disparity errors across the estimated trait spaces.
* **07_continuous_lmm.Rmd** & **08_continuous_lmm.Rmd** run the linear mixed models on the disparity errors.

---

## Vignette

```{R}
library(treats)

source("Analyses/01_tree_generation.R")




```
