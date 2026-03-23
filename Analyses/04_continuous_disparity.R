library(treats)
if(packageVersion("treats") != "1.1.7") {
    library(devtools)
    install_github("TGuillerme/treats", ref = "map.traits.events") ## requires this branch of treats
}

if(packageVersion("dispRity") != "1.9.8") {
    library(devtools)
    install_github("TGuillerme/dispRity")
}

dir.create("../Data/continuous/disparity", recursive = TRUE, showWarnings = FALSE)
source("../Functions/utility.R")



run.cont.disparity <- function(tree_size, replicate_id){

    ord_pre_sample <- readRDS(sprintf("../Data/continuous/ord/ord_pre_sample_%st_%03d.rds", tree_size, replicate_id))
    ord_pre_point <- readRDS(sprintf("../Data/continuous/ord/ord_pre_point_%st_%03d.rds", tree_size, replicate_id))
    ord_no_ace <- readRDS( sprintf("../Data/continuous/ord/ord_no_ace_%st_%03d.rds", tree_size, replicate_id))
    ord_true <- readRDS( sprintf("../Data/continuous/ord/ord_true_%st_%03d.rds", tree_size, replicate_id))
    post_ord_ace <- readRDS(sprintf("../Data/continuous/ace/post_ord_ace_%st_%03d.rds", tree_size, replicate_id))
    point_post_ord_ace <- readRDS(sprintf("../Data/continuous/ord/ord_post_point_%st_%03d.rds", tree_size, replicate_id))
    sample_post_ord_ace <- readRDS(sprintf("../Data/continuous/ord/ord_post_sample_%st_%03d.rds", tree_size, replicate_id))


    calc.error <- function(estimates, true_vals, metric) {
    Map(function(rate_est, rate_true) {
        lapply(rate_est, function(fossil_est) {
        # Handle both point and sample estimations
        if(is.list(fossil_est) && !is.matrix(fossil_est)) {
            # Is a sample method, nested list
            lapply(fossil_est, function(sample) {
            est <- get.disparity(dispRity(sample, metric = metric))
            (est[[1]] - rate_true[[1]]) / rate_true[[1]]
            })
        } else {
            est <- get.disparity(dispRity(fossil_est, metric = metric))
            (est[[1]] - rate_true[[1]]) / rate_true[[1]]
        }
        })
    }, estimates, true_vals)
    }

    # Define metrics
    metrics <- list(
    sum_var = c(sum, variances),
    sum_quant = c(sum, quantiles),
    pairwise = c(mean, pairwise.dist.na.rm)
    )

    results_raw <- lapply(names(metrics), function(metric_name) {
    metric <- metrics[[metric_name]]
    true_disp <- lapply(ord_true, function(rate) get.disparity(dispRity(rate, metric = metric)))
    
    list(
        pre_ord_sample = calc.error(ord_pre_sample, true_disp, metric),
        pre_ord_point = calc.error(ord_pre_point, true_disp, metric),
        no_ace = calc.error(ord_no_ace, true_disp, metric),
        post_ord_point = calc.error(point_post_ord_ace, true_disp, metric),
        post_ord_sample = calc.error(sample_post_ord_ace, true_disp, metric)
    )
    })
    names(results_raw) <- names(metrics)


    saveRDS(lapply(results_raw, `[[`, "pre_ord_sample"), sprintf("../Data/continuous/disparity/pre_ord_sample_disparity_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(lapply(results_raw, `[[`, "pre_ord_point"), sprintf("../Data/continuous/disparity/pre_ord_point_disparity_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(lapply(results_raw, `[[`, "no_ace"), sprintf("../Data/continuous/disparity/no_ace_disparity_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(lapply(results_raw, `[[`, "post_ord_point"), sprintf("../Data/continuous/disparity/post_ord_point_disparity_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(lapply(results_raw, `[[`, "post_ord_sample"), sprintf("../Data/continuous/disparity/post_ord_sample_disparity_%st_%03d.rds", tree_size, replicate_id))
}

## using same tree_size/replicate id as `02_continuous_sims_ace.R`
run.cont.disparity(tree_size = 50, replicate_id = 1)

## Full reproducibility
tree_sizes <- c(50, 100, 150)
for(tree_size in tree_sizes) {
    for(i in 1:100) {
        run.cont.disparity(tree_size, i)
    }
}
