create.matched.data <- function(true_vals, ace_vals_clean) {
  true_matched <- c()
  ace_matched <- c()
  fossil_levels <- c()  # Add this line
  
  for(method_name in names(ace_vals_clean)) {
    for(rate_name in names(ace_vals_clean[[method_name]])) {
      for(rep_idx in seq_along(ace_vals_clean[[method_name]][[rate_name]])) {
        
        true_val <- true_vals[[rate_name]][[rep_idx]]
        ace_rep_vals <- ace_vals_clean[[method_name]][[rate_name]][[rep_idx]]
        
        # Get fossil level names and repeat true values
        fossil_names <- names(ace_rep_vals)
        n_fossil_levels <- length(ace_rep_vals)
        
        true_matched <- c(true_matched, rep(true_val, n_fossil_levels))
        ace_matched <- c(ace_matched, unlist(ace_rep_vals))
        fossil_levels <- c(fossil_levels, fossil_names)  # Add this line
      }
    }
  }
  
  return(list(true = true_matched, ace = ace_matched, fossil_level = fossil_levels))
}

fix.zero.branches <- function(tree, min_length = 1e-8) {
  tree$edge.length[tree$edge.length <= 0] <- min_length
  return(tree)
}

disp.diff <- function(ace, true){
  tryCatch({
    if(is.null(ace[[1]]$elements)){
      return(lapply(ace, function(x){
        disp.diff(x, true)
      }))
    } 
    diff <- (ace[[1]]$elements[1] - true[[1]]$elements[1]) / true[[1]]$elements[1]
    return(diff)
  }, error = function(e) {
    cat("Error in disp.diff:", e$message, "\n")
    cat("ace class:", class(ace), "true class:", class(true), "\n")
    return(NA)
  })
}

# write.path <- function(subfolder, filename) {
#   return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
# }



write.path <- function(folder, subfolder, filename){
  file.path("..", "Data", "folder", "subfolder", sprintf(filename, replicate_id))
}

BM.trend.process <- function(x0 = 0, edge.length = 1,
                             Sigma = diag(length(x0)), trend = 0.1, ...) {
  drift <- trend * edge.length 
  if (edge.length < (0.01 * Sigma[1,1])){ drift <- 0}
  return(t(MASS::mvrnorm(n = 1, mu = x0 + drift, Sigma = Sigma * edge.length, ...)))
}

calc.error <- function(estimates, true_vals, metric_fun) {
  Map(function(rate_est, rate_true) {
    lapply(rate_est, function(fossil_est) {
      # Handle both single estimates and lists of samples
      if(is.list(fossil_est) && !is.data.frame(fossil_est) && length(fossil_est) > 1) {
        # This is a list of samples (100 replicates)
        lapply(fossil_est, function(sample) {
          est <- get.disparity(dispRity(sample, metric = metric_fun))[[1]]
          (est - rate_true[[1]]) / rate_true[[1]]
        })
      } else {
        # This is a single estimate
        est <- get.disparity(dispRity(fossil_est, metric = metric_fun))[[1]]
        (est - rate_true[[1]]) / rate_true[[1]]
      }
    })
  }, estimates, true_vals)
}

pairwise.dist.na.rm <- function(matrix, method = "euclidean", ...) {
 
        # Calculate distances on clean data
        distances <- vegan::vegdist(matrix, method = method, na.rm = TRUE, ...)
        return(as.vector(distances))
}