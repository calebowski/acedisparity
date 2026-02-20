library(lme4)
library(lmerTest) 
library(emmeans)
library(multcomp)
library(multcompView)



tree_sizes <- c(50, 100, 150)
methods <- c("pre_ord_point", "pre_ord_sample", "post_ord_point", "post_ord_sample", "no_ace")
n_replicates <- 100 ## if doing single replicate change to 1

results <- list()
for (size in tree_sizes) {
  results[[size]] <- list()
  for (method in methods) {
    results[[size]][[method]] <- list()
    for(i in 1:n_replicates) {
      file_path <- file.path("..", "Data", "continuous", "disparity", sprintf("%s_disparity_%st_%03d.rds", method, size, i))
      if(file.exists(file_path)) {
        results[[size]][[method]][[i]] <- readRDS(file_path)
      } else {
        warning("Missing file: ", file_path)
        results[[size]][[method]][[i]] <- NULL
      }
    }
  }
}


metric_names <- c("sum_quant", "sum_var", "pairwise")

# Transpose structure: metric -> tree_size -> method -> replicate
results_relisted <- setNames(lapply(metric_names, function(metric) {
  setNames(lapply(tree_sizes, function(size) {
    setNames(lapply(methods, function(method) {
      is_sample_method <- grepl("sample", method)
      
      lapply(results[[size]][[method]], function(rep_data) {
        metric_data <- rep_data[[metric]]
        
        if(is_sample_method) {
          lapply(metric_data, function(model_data) {
            lapply(model_data, function(fossil_data) {
              if(is.list(fossil_data) && !is.data.frame(fossil_data)) {
                unlist(fossil_data)
              } else {
                fossil_data
              }
            })
          })
        } else {
          metric_data
        }
      })
    }), methods)
  }), tree_sizes)
}), metric_names)