---
title: "Continuous LMM Analysis"
output: html_document
---

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



# Build dataframe with tree_size
results_df <- do.call(rbind, lapply(names(results_relisted), function(metric_name) {
  metric_data <- results_relisted[[metric_name]]
  
  do.call(rbind, lapply(names(metric_data), function(size_name) {
    size_data <- metric_data[[size_name]]
    
    do.call(rbind, lapply(names(size_data), function(method_name) {
      method_data <- size_data[[method_name]]
      
      do.call(rbind, lapply(seq_along(method_data), function(rep_idx) {
        rep_data <- method_data[[rep_idx]]
        
        if(is.null(rep_data)) return(NULL)
        
        do.call(rbind, lapply(names(rep_data), function(model_name) {
          model_data <- rep_data[[model_name]]
          
          if(is.null(model_data)) return(NULL)
          
          data.frame(
            replicate = rep_idx,
            all = I(list(model_data$all)),
            high = I(list(model_data$fossil_high)),
            med = I(list(model_data$fossil_med)),
            low = I(list(model_data$fossil_low)),
            living = I(list(model_data$living)),
            model = model_name,
            method = method_name,
            metric = metric_name,
            tree_size = size_name,  #  Add tree size
            stringsAsFactors = FALSE
          )
        }))
      }))
    }))
  }))
}))



# Pivot to long format using base R - FIXED VERSION
results_df_long <- do.call(rbind, lapply(1:nrow(results_df), function(i) {
  row <- results_df[i, ]
  
  # Extract each fossil sampling level separately
  fossil_levels <- c("all", "high", "med", "low", "living")
  
  # Build data frame by binding each fossil level
  do.call(rbind, lapply(fossil_levels, function(fossil_level) {
    values <- unlist(row[[fossil_level]])
    if(length(values) > 0) {
      data.frame(
        tree_size = row$tree_size,
        replicate = row$replicate,
        model = row$model,
        method = row$method,
        metric = row$metric,
        fossil_sampling = fossil_level,
        error = values,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    } else {
      NULL
    }
  }))
}))


# Aggregate sample methods, keep point methods as-is
results_df_long$is_sample <- grepl("sample", results_df_long$method)

# Split data by sample vs point methods
sample_data <- results_df_long[results_df_long$is_sample, ]
point_data <- results_df_long[!results_df_long$is_sample, ]

# Aggregate only sample methods
sample_aggregated <- aggregate(error ~ tree_size + replicate + model + method + metric + fossil_sampling,
                               data = sample_data, 
                               FUN = median)
