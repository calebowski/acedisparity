args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
# library(parallel)
library(treats)

cat("Starting replicate", replicate_id, "\n")


set.seed(100 + replicate_id)

file_path <- "/mnt/parscratch/users/bip24cns/acedisparity/continuous/50t/ord/"

ord_sample <- readRDS(paste0(file_path, sprintf("ord/sample_ord_%03d.rds", replicate_id)))
ord_point <- readRDS(paste0(file_path, sprintf("ord/point_ord_%03d.rds", replicate_id)))
ord_no_ace <- readRDS(paste0(file_path, sprintf("ord/no_ace_ord_%03d.rds", replicate_id)))
ord_true <- readRDS(paste0(file_path, sprintf("ord/true_ord_%03d.rds", replicate_id)))
point_post_ord_ace_living <- readRDS(paste0(file_path, sprintf("ord/point_post_ord_%03d.rds", replicate_id)))
sample_post_ord_ace_living <- readRDS(paste0(file_path, sprintf("ord/sample_post_ord_%03d.rds", replicate_id)))


metrics_list <- list(
  sum_var = c(sum, variances),
  sum_quant = c(sum, quantiles),
  mean_neigh = c(mean, neighbours)
)

# Calculate disparity for all metrics
for(metric_name in names(metrics_list)) {
  metric <- metrics_list[[metric_name]]
  
  # Apply to each ordination type with appropriate structure
  assign(paste0(metric_name, "_sample"), 
         lapply(ord_sample, lapply, lapply, function(mat) dispRity(mat, metric = metric)$disparity))
  
  assign(paste0(metric_name, "_point"), 
         lapply(ord_point, lapply, function(mat) dispRity(mat, metric = metric)$disparity))
  
  assign(paste0(metric_name, "_no_ace"), 
         lapply(ord_no_ace, lapply, function(mat) dispRity(mat, metric = metric)$disparity))
  
  assign(paste0(metric_name, "_true"), 
         lapply(ord_true, function(mat) dispRity(mat, metric = metric)$disparity))
  
  assign(paste0(metric_name, "_point_postord"), 
         lapply(point_post_ord_ace_living, lapply, function(mat) dispRity(mat, metric = metric)$disparity))

  assign(paste0(metric_name, "_sample_postord"), 
         lapply(sample_post_ord_ace_living, lapply, lapply, function(mat) dispRity(mat, metric = metric)$disparity))
}

cat("Disparity calculation finished...\n")

for(metric_name in names(metrics_list)) {
  for(ord_type in c("sample", "point", "point_postord", "sample_postord", "no_ace", "true")) {
    var_name <- paste0(metric_name, "_", ord_type)
    saveRDS(get(var_name), paste0(file_path, sprintf("disparity/%s_%03d.rds", var_name, replicate_id)))
  }
}

cat("Disparity files saved...\n")



est_disp <- list(
  point = mget(paste0(c("sum_var", "sum_quant", "mean_neigh"), "_point")),
  sample = mget(paste0(c("sum_var", "sum_quant", "mean_neigh"), "_sample")),
  no_ace = mget(paste0(c("sum_var", "sum_quant", "mean_neigh"), "_no_ace")),
  point_postord = mget(paste0(c("sum_var", "sum_quant", "mean_neigh"), "_point_postord")),
  sample_postord = mget(paste0(c("sum_var", "sum_quant", "mean_neigh"), "_sample_postord"))
)

true_disp <- mget(paste0(c("sum_var", "sum_quant", "mean_neigh"), "_true"))



diffs <- lapply(est_disp, function(method){
  method_result <- Map(function(est_metric, true_metric) {
    Map(function(est_model, true_model){
      true_disp <- true_model[[1]]$elements[1]
      lapply(est_model, function(fossil_lev){
        if (is.list(fossil_lev) && !is.null(fossil_lev[[1]]$elements)) {
        diff <- (fossil_lev[[1]]$elements[1] - true_disp) / true_disp ## calculate relative disparity difference 
        return(diff)
        } else {lapply(fossil_lev, function(rep) { ## for sample methods with additional nested list
          diff <- (rep[[1]]$elements[1] - true_disp) / true_disp 
          return(diff)
          })
        }
      })
    }, est_metric, true_metric)
  }, method, true_disp)
})

cat("Disparity diffs calculated...\n")


methods <- names(diffs)
for (method in methods){
  saveRDS(diffs[[method]], paste0(file_path, sprintf("disparity/diff_%s_%03d.rds", method, replicate_id)))
}

cat("Disparity diffs saved...\n")


cat("=== CHECKING WARNINGS ===\n")
warning_list <- warnings()
if(!is.null(warning_list)) {
  cat("Number of warnings:", length(warning_list), "\n")
  print(warning_list)
} else {
  cat("No warnings detected\n")
}
cat("Script completed\n")