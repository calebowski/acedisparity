args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
tree_size <- args[2]
job_id <- args[3]

base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/discrete/crown/", tree_size, "/")

set.seed(100 + replicate_id)



write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}


fossil_matrices <- readRDS(write.path("matrices", "fossil_matrices_%03d.rds"))



tasks  <- expand.grid(rate = names(fossil_matrices), fossil_level = names(fossil_matrices[[1]]),  stringsAsFactors = FALSE) 

res_pre_ord_ace <- mclapply(seq_len(nrow(tasks)), function(i){ ## loop over each model combination by row
    task <- tasks[i,]
    level <- fossil_matrices[[task$rate]][[task$fossil_level]]
    tryCatch({
    multi.ace(level$matrix, level$tree, models = "ER", output = "multi.ace")
  }, error = function(e) {
    cat("ERROR:", task$fossil_level, e$message, "\n") ## error handeling
    NULL
  })
}, mc.cores = 15)

pre_ord_ace <- list()
for(i in seq_along(res_pre_ord_ace)) {
  r <- tasks$rate[i]  
  l <- tasks$fossil_level[i]
  pre_ord_ace[[r]][[l]] <- res_pre_ord_ace[[i]]
}

point_fossil_anc <- lapply(pre_ord_ace, lapply, multi.ace, threshold = "max_tiebreaker", output = "combined.matrix", verbose = TRUE)


saveRDS(point_fossil_anc, write.path("anc", "pre_ord_tiebreaker_%03d.rds"))

ord_point <- lapply(point_fossil_anc, lapply,  function(rep){
dist <- char.diff(rep, method = "mord", by.col = FALSE)
ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})
saveRDS(ord_point, write.path("ord", "ord_tiebreaker_%03d.rds"))


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

# Ensure disparity directory exists
raw_disparity_dir <- paste0(base_path, "disparity/raw")
if(!dir.exists(raw_disparity_dir)) dir.create(raw_disparity_dir, recursive = TRUE)

cat("Saving raw disparity results...\n")

# Save by treatment type
saveRDS(lapply(results_raw, `[[`, "pre_ord_sample"), write.path("disparity/raw", "pre_ord_sample_%03d.rds"))
saveRDS(lapply(results_raw, `[[`, "pre_ord_point"), write.path("disparity/raw", "pre_ord_point_%03d.rds"))
saveRDS(lapply(results_raw, `[[`, "no_ace"), write.path("disparity/raw", "no_ace_%03d.rds"))
saveRDS(lapply(results_raw, `[[`, "post_ord_point"), write.path("disparity/raw", "post_ord_point_%03d.rds"))
saveRDS(lapply(results_raw, `[[`, "post_ord_sample"), write.path("disparity/raw", "post_ord_sample_%03d.rds"))


