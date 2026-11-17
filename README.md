# When does ancestral state estimation improve disparity analyses? A simulation-based evaluation

Authors: [Caleb Scutt](mailto:cnscutt1@sheffield.ac.uk), [Natalie Cooper](https://github.com/nhcooper123), [Gavin Thomas](https://github.com/ghthomas) and [Thomas Guillerme](https://github.com/TGuillerme)

[Preprint here](https://www.biorxiv.org/content/10.64898/2026.04.21.719166v1)

This study uses simulations to evaluate different methods of ancestral state estimation for recovering morphological disparity under varying evolutionary models, fossil sampling densities, and trait types.
This repository contains the scripts to reproduce the analyses of the paper.

<!-- ---

## Data

The `Data/` folder contains the trees and simulated matrices that are used in the analyses. These are also easily reproduced using the `Analyses/` folder. 
The LMM outputs are on [Dryad](https://doi.org/10.5061/dryad.280gb5n4p)

--- -->

## Analyses

All code used to run the analyses are contained in this repo. 
The `Analyses/` folder contains the scripts to reproduce the analyses, functionalised so that you can reproduce the analyses with however many replicates you choose. Ensure your working directory is set to `Analyses/`.
Within `Analyses/sensitivity_analyses`, there are scripts to complete the sensitivity analyses (namely running the pipeline with asymmetrical trees and sampling both tips+nodes during fossil sampling simulation). These scripts are prefixed with "S" and then a number, i.e., "S02_tip_nodes_fossil_sample_continuous_sims_ace.R", with the number corresponding which script in the main analysis it replaces.
The `Cluster/` folder contains scripts used to produce the analyses on the Sheffield HPC.
Note, that the analyses take over 1.5 CPU years to run.

Sections 01-05 are formatted so that you can choose to run for either one replicate, or the whole analysis.

* **01_tree_generation.R** generates the 300 trees (100 for each tree size level) used in both continuous and discrete analyses.
* **02_continuous_sims_ace.R** & **03_discrete_sims_ace.R** simulates continuous and discrete traits across the trees, simulates fossil sampling, runs pre/post-ordination and point/distribution ancestral state estimation and generates trait spaces.
* **04_continuous_disparity.R** & **05_discrete_disparity.R** calculates disparity errors across the estimated continuous and discrete trait spaces.
* **06_continuous_lmm.Rmd** & **07_discrete_lmm.Rmd** runs the aggregated and weighted linear mixed models on the continuous and discrete disparity errors. It also plots the heatmaps of estimated marginal means from the three-way model.
* **08_continuous_plot.R** & **09_discrete_plot.R** generate the boxplots (main figures) of relative disparity errors for combinations of ASE method, fossil sampling, model/transition rate and disparity metric.

---

## Packages & Functions

The majority of the functions used in this paper come from the [`treats`](https://github.com/TGuillerme/treats/tree/map.traits.events) & [`dispRity`](https://github.com/TGuillerme/dispRity) packages. Note that the paper uses the **latest development versions of these packages** (as of March 2026). To use the exact versions of these packages used in this study, install with this:

```r
remotes::install_github("TGuillerme/treats", ref = "a02738522a05eb0f4da29806011acf9cbb4a6a83")

remotes::install_github("TGuillerme/dispRity", ref = "753bddb3da93f1e067e70fb1de7a120aa73e385a")
```

There are also custom functions used for simulating fossil sampling, contained within `Functions/`.

---


### Getting sampled distributions instead of point estimates from the estimations

You can also use `multi.ace` to sample distributions instead of getting point estimates.
For example, by default, when using `multi.ace` on continuous characters, the output estimation is the average most likely value for the trait at every specific nodes.

```r
set.seed(1)
## Creating a coalescent tree
my_tree <- makeNodeLabel(rcoal(15))

## Creating a matrix of continuous characters
data <- space.maker(elements = 15, dimensions = 1, distribution = rnorm,
                    elements.name = my_tree$tip.label)

## Running the estimation with details
point_estimates <- multi.ace(data, my_tree, estimation.details = list("CI95"))
```

For example for trait 1 for Node 1, the estimated trait value is `r unname(round(point_estimates$estimations[1,1], 3))` but with a confidence interval of `r unname(round(point_estimates$details[[1]]$CI95[1,1], 3))` and `r unname(round(point_estimates$details[[1]]$CI95[1,2], 3))` which is quiet the spread.
To avoid treating the ancestral estimation as a point estimate (here the average of a large spread) we can also just get a number of samples to get an estimation of the distribution instead of the point estimate using the `ml.collapse` argument and `type = "sample"`.

```r
## Getting 1000 samples from the distribution
distribution_estimates <- multi.ace(data, my_tree, ml.collapse = list(type = "sample", sample = 1000))
```

This results in 1000 matrices that are sampled by default uniformaly.
We can visualise that for Node 1 that will be randomly sampled between the estimate's 95% confidence interval highlighted above.

```{r, fig.height = 4, fig.width = 4}
## Visualising the sampling of Node 1
node_1_estimates <- unlist(lapply(distribution_estimates, function(x) x[1,1]))
hist(node_1_estimates, main = "Distribution of estimates for Node 1")
# Note that the two small bars are because of the
# last and first bins should start/end at the actual 95 CI
```

Of course, we can change the type of function to use for sampling knowing the confidence interval of our estimations.
We can specify this by providing a named list containing the following two elements:

 1. `fun` that is the function to use for the sampling (`rnorm` in the example below);
 2. `param` that contains a list of functions or variables on how to calculate the extra arguments for `fun`. For example for `rnorm`, the function needs 1) a `mean` argument (here `mean = mean`, to use the function `mean` to calculate... the mean); and 2) a `sd` argument (here `sd = approx.sd`, a function to approximate the stnadard deviation knowing the confidence interval - assuming normality of the data).

```r
## Function for approximating the standard deviation
## from the range of a confidence interval.
approx.sd <- function(CI)return(diff(range(CI))/4)

## Getting 1000 samples from the distribution using a normal one
distribution_estimates <- multi.ace(data, my_tree, ml.collapse = list(type = "sample", sample = 1000,
                            sample.fun = list(fun = rnorm,
                                              param = list(mean = mean,
                                                           sd = approx.sd))))
```

And again, we can visualise that:

```{r, fig.height = 4, fig.width = 4}
## Visualising the sampling of Node 1
node_1_estimates <- unlist(lapply(distribution_estimates, function(x) x[1,1]))
hist(node_1_estimates, main = "Distribution of estimates for Node 1")
# Which is now a neat normal distribution
```

In this example we used the default and a specified `sample.fun` argument for a matrix with only one trait.
However you can use any number of distributions for any number of traits: either just the one to be applied to all traits or a list of the same length as the number of traits with a list of `fun` and `param` for each.

To play around with different sampling strategies, but if your estimations take a lot of time, you can save your estimations as a `"multi.ace"` object that you can then recycle and reuse very easily:

```r
## Saving the estimation without getting any estimate (raw):
raw_estimates <- multi.ace(data, my_tree, output = "multi.ace")

## Sampling from the raw_estimates
sampled_estimates <- multi.ace(raw_estimates, ml.collapse = list(type = "sample", sample = 100))
```