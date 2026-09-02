args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
tree_size <- args[2]
job_id <- args[3]
source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")

library(dispRity)

base_path <- "/mnt/parscratch/users/bip24cns/acedisparity/revisions/discrete/"

write.path <- function(subfolder, filename) {
  paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id))
}

rates <- c("slow", "med", "fast")
fossils <- c("fossil_high", "fossil_med", "fossil_low", "living")

# LOAD DATA

cat("Loading data...\n")

# Load pre-ordination ACE
point_pre_ord_ace <- readRDS(write.path("ord", "ord_point_%03d.rds"))

sample_pre_ord_ace <- lapply(setNames(rates, rates), function(rate) {
  lapply(setNames(fossils, fossils), function(fossil) {
    lapply(1:100, function(i) {
      readRDS(paste0(base_path, "checkpoints_ord/rep", sprintf("%03d", replicate_id), 
                     "/", rate, "/", fossil, "/", i, ".rds"))
    })
  })
})

# Load other data
ord_no_ace <- readRDS(write.path("ord", "ord_no_ace_%03d.rds"))
ord_true <- readRDS(write.path("ord", "ord_true_%03d.rds"))

# Load post-ordination ACE
point_post_ord_ace <- lapply(setNames(rates, rates), function(rate) {
  lapply(setNames(fossils, fossils), function(fossil) {
    readRDS(paste0(base_path, "checkpoints_ace/rep", sprintf("%03d", replicate_id), 
                   "/", rate, "/", fossil, "/", job_id, "_point_post_ord_ace.rds"))
  })
})

sample_post_ord_ace <- lapply(setNames(rates, rates), function(rate) {
  lapply(setNames(fossils, fossils), function(fossil) {
    readRDS(paste0(base_path, "checkpoints_ace/rep", sprintf("%03d", replicate_id), 
                   "/", rate, "/", fossil, "/", job_id, "_sample_post_ord_ace.rds"))
  })
})


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

################################################################################
# CALCULATE ERRORS FOR ALL METRICS USING ALL AXES
################################################################################

cat("Calculating raw disparity errors...\n")

# Define metrics
metrics <- list(
  sum_var = c(sum, variances),
  sum_quant = c(sum, quantiles),
  pairwise = c(mean, pairwise.dist.na.rm)
)

# Fix: Calculate with correct metric names
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

raw_disparity_dir <- paste0(base_path, "disparity")
if(!dir.exists(raw_disparity_dir)) dir.create(raw_disparity_dir, recursive = TRUE)

cat("Saving raw disparity results...\n")

# Save by treatment type
saveRDS(lapply(results_raw, `[[`, "pre_ord_sample"), write.path("disparity", "pre_ord_sample_%03d.rds"))
saveRDS(lapply(results_raw, `[[`, "pre_ord_point"), write.path("disparity", "pre_ord_point_%03d.rds"))
saveRDS(lapply(results_raw, `[[`, "no_ace"), write.path("disparity", "no_ace_%03d.rds"))
saveRDS(lapply(results_raw, `[[`, "post_ord_point"), write.path("disparity", "post_ord_point_%03d.rds"))
saveRDS(lapply(results_raw, `[[`, "post_ord_sample"), write.path("disparity", "post_ord_sample_%03d.rds"))

cat("Completed raw...!\n")
