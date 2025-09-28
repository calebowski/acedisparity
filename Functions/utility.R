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
  if(is.null(ace[[1]]$elements )){
    return(lapply(ace, function(x){
      disp.diff(x, true)
    }))
  } 
  diff <- (ace[[1]]$elements[1] - true[[1]]$elements[1]) / true[[1]]$elements[1]
  return(diff)
}

write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}