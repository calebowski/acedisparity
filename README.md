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
The `Analyses/` folder contains the scripts to reproduce the analyses, functionalised so that you can reproduce the analyses with however many replicates you choose. Ensure your working directory is set to `Analyses/`.
The `Cluster/` folder contains scripts used to produce the analyses on the Sheffield HPC.
Note, that the analyses take over 1.5 CPU years to run.

Sections 01-05 are formatted so that you can choose to run for either one replicate, or the whole shebang.

* **01_tree_generation.R** generates the 300 trees (100 for each tree size level) used in both continuous and discrete analyses.
* **02_continuous_sims_ace.R** & **03_discrete_sims_ace.R** simulates continuous and discrete traits across the trees, simulates fossil sampling, runs pre/post-ordination and point/distribution ancestral state estimation and generates trait spaces.
* **04_continuous_disparity.R** & **05_discrete_disparity.R** calculates disparity errors across the estimated continuous and discrete trait spaces.
* **06_continuous_lmm.Rmd** & **07_continuous_lmm.Rmd** runs the aggregated and weighted linear mixed models on the continuous and discrete disparity errors. It also plots the heatmaps of estimated marginal means from the three-way model.

---

## Packages & Functions

The majority of the functions used in this paper come from the [`treats`](https://github.com/TGuillerme/treats/tree/map.traits.events) & [`dispRity`](https://github.com/TGuillerme/dispRity) packages. Note that the paper uses the **latest development versions of these packages** (as of March 2026). To use the exact versions of these packages used in this study, install with this:

```{R}
remotes::install_github("TGuillerme/treats", ref = "a02738522a05eb0f4da29806011acf9cbb4a6a83")

remotes::install_github("TGuillerme/dispRity", ref = "753bddb3da93f1e067e70fb1de7a120aa73e385a")
```

There are also custom functions used for simulating fossil sampling, contained within `Functions/`.

---
<!-- 
## Vignette

```{R}
library(treats)

tree <- read.tree("Data/trees/tree_50t_001.tre")
plot(tree)



``` -->
