args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
tree_size <- args[2]
job_id <- args[3]

library(dispRity)

base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/continuous/", tree_size, "/")

write.path <- function(subfolder, filename) {
  paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id))
}
source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")
# models <- c("bm", "bm_t", "ou_w", "ou_st", "ou_sh")
# fossils <- c("all", "fossil_high", "fossil_med", "fossil_low", "living")

## Load data

cat("Loading data...\n")

# Load ordination matrices
point_pre_ord_ace <- readRDS(write.path("ord", "ord_point_%03d.rds"))
sample_pre_ord_ace <- readRDS(write.path("ord", "ord_sample_%03d.rds"))
ord_no_ace <- readRDS(write.path("ord", "ord_no_ace_%03d.rds"))
ord_true <- readRDS(write.path("ord", "ord_true_%03d.rds"))
point_post_ord_ace <-  readRDS(write.path("ord", "post_ord_point_%03d.rds"))
sample_post_ord_ace <-  readRDS(write.path("ord", "post_ord_sample_%03d.rds"))


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
    pre_ord_sample = calc.error(sample_pre_ord_ace, true_disp, metric),
    pre_ord_point = calc.error(point_pre_ord_ace, true_disp, metric),
    no_ace = calc.error(ord_no_ace, true_disp, metric),
    post_ord_point = calc.error(point_post_ord_ace, true_disp, metric),
    post_ord_sample = calc.error(sample_post_ord_ace, true_disp, metric)
  )
})
names(results_raw) <- names(metrics)

# Ensure disparity directory exists
raw_disparity_dir <- paste0(base_path, "disparity")
if(!dir.exists(raw_disparity_dir)) dir.create(raw_disparity_dir, recursive = TRUE)


saveRDS(lapply(results_raw, `[[`, "pre_ord_sample"), write.path("disparity", "pre_ord_sample_%03d.rds"))
saveRDS(lapply(results_raw, `[[`, "pre_ord_point"), write.path("disparity", "pre_ord_point_%03d.rds"))
saveRDS(lapply(results_raw, `[[`, "no_ace"), write.path("disparity", "no_ace_%03d.rds"))
saveRDS(lapply(results_raw, `[[`, "post_ord_point"), write.path("disparity", "post_ord_point_%03d.rds"))
saveRDS(lapply(results_raw, `[[`, "post_ord_sample"), write.path("disparity", "post_ord_sample_%03d.rds"))

cat("Finished replicate...\n")
