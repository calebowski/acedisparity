library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest) 
library(emmeans)
library(multcomp)
library(multcompView)
library(xtable)
# Define job IDs for each tree size
# job_ids <- list(
#   "50t" = "8558401",
#   "100t" = "8561556",
#   "150t" = "8561712"
# )

    
job_ids <- list(
  "50t" = "8558401",
  "100t" = "8561556",  # Replace with actual job ID
  "150t" = "8561712"   # Replace with actual job ID
)

tree_sizes <- c("50t", "100t", "150t")
methods <- c("pre_ord_point_tiebreaker", "pre_ord_sample", "post_ord_point",
             "post_ord_sample", "no_ace")

results <- list()
for (size in tree_sizes) {
  results[[size]] <- list()
  job_id <- job_ids[[size]]  #  Get job ID for this tree size
  
  for (method in methods) {
    results[[size]][[method]] <- list()
    for(i in 1:100) {
      file_path <- file.path("..", "Data", "cluster", "discrete_crown", size, 
                            "disparity", "raw",
                            sprintf("%s_%s_%03d.rds", job_id, method, i))      
      if(file.exists(file_path)) {
        results[[size]][[method]][[i]] <- readRDS(file_path)
      } else {
        warning("Missing file: ", file_path)
        results[[size]][[method]][[i]] <- NULL
      }
    }
  }
}

metric_names <- c("sum_var", "sum_quant", "pairwise")

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

rm(results)
gc()

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

rm(results_relisted)
gc()
# Aggregate sample methods, keep point methods as-is
results_df_long$is_sample <- grepl("sample", results_df_long$method)

# Split data by sample vs point methods
sample_data <- results_df_long[results_df_long$is_sample, ]
point_data <- results_df_long[!results_df_long$is_sample, ]

# Aggregate only sample methods
sample_aggregated <- aggregate(error ~ tree_size + replicate + model + method + metric + fossil_sampling,
                               data = sample_data, 
                               FUN = median)

# Combine back together
results_df_long <- rbind(
  point_data[, c("tree_size", "replicate", "model", "method", "metric", "fossil_sampling", "error")],
  sample_aggregated
)

# Factor variables
results_df_long$tree_size <- factor(results_df_long$tree_size)
results_df_long$method <- factor(results_df_long$method)
results_df_long$fossil_sampling <- factor(results_df_long$fossil_sampling)
results_df_long$model <- factor(results_df_long$model)
results_df_long$metric <- factor(results_df_long$metric)

results_df_long$tree_unique_id <- interaction(results_df_long$tree_size, 
                                              results_df_long$replicate)

table(results_df_long$method)


results_df_long$abs_error <- abs(results_df_long$error)
results_df_long$log_abs_error <- log(results_df_long$abs_error + 0.001)


lmm_model <- lmer(log_abs_error ~ model * method * fossil_sampling * metric + tree_size +
                                  (1|tree_unique_id),
                    data = results_df_long, REML = TRUE)
saveRDS(lmm_model, "../Data/cluster/discrete_crown/lmm/raw/lmm_model_four_way.rds")
lmm_model <- readRDS("../Data/cluster/discrete_crown/lmm/raw/lmm_model_four_way.rds")



#################################################################################################
lmm_path <- "/home/caleb/Documents/PhD/acedisparity/Data/cluster/discrete_crown/lmm/raw/"

# Residual plots
pdf(paste0(lmm_path, "diagnostics/diagnostic_plot.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2))

# QQ plot
qqnorm(resid(lmm_model), main = "Q-Q Plot of Residuals")
qqline(resid(lmm_model), col = "red")

# Residuals vs Fitted
plot(fitted(lmm_model), resid(lmm_model),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# Scale-Location
plot(fitted(lmm_model), sqrt(abs(resid(lmm_model))),
     xlab = "Fitted Values", ylab = "Square Root of |Residuals|",
     main = "Scale-Location")

# Histogram of residuals
hist(resid(lmm_model), breaks = 50, main = "Histogram of residuals",
     xlab = "Residuals")
dev.off()


#################################################################################################

sink(paste0(lmm_path, "model_summaries_four_way.txt"))
cat("\n Model Summary \n")
print(summary(lmm_model))

cat("\n Model fit \n")
cat("AIC:", AIC(lmm_model), "\n")
cat("BIC:", BIC(lmm_model), "\n")
cat("Log-Likelihood:", logLik(lmm_model), "\n")

cat("\n Random effects variance \n")
print(VarCorr(lmm_model))

# ANOVA table type 3 due to unequal sample sizes
cat("\n ANOVA table type 3 due to unequal sample sizes \n")
print(anova(lmm_model, type = 3))
sink()

#################################################################################################

# Estimated marginal means and Tukey comparisons
cat("\n--- Rate Rankings ---\n")
emm_model <- emmeans(lmm_model, ~ model)
cld_model <- cld(emm_model, Letters = letters, alpha = 0.05)
write.csv(cld_model, paste0(lmm_path, "tables/model_multcomp.csv"))
print(xtable(cld_model), 
      include.rownames = FALSE)

cat("\n--- Method Rankings ---\n")
emm_method <- emmeans(lmm_model, ~ method)
cld_method <- cld(emm_method, Letters = letters, alpha = 0.05)
write.csv(cld_method, paste0(lmm_path, "tables/method_multcomp.csv"))
print(xtable(cld_method), 
      include.rownames = FALSE)


cat("\n--- Fossil sampling Level Rankings ---\n")
emm_fossils <- emmeans(lmm_model, ~ fossil_sampling)
cld_fossils <- cld(emm_fossils, Letters = letters, alpha = 0.05)
write.csv(cld_fossils, paste0(lmm_path, "tables/fossil_multcomp.csv"))
print(xtable(cld_fossils), include.rownames = FALSE)

cat("\n--- Metric Level Rankings ---\n")
emm_metric <- emmeans(lmm_model, ~ metric)
cld_metric <- cld(emm_metric, Letters = letters, alpha = 0.05)
write.csv(cld_metric, paste0(lmm_path, "tables/metric_multcomp.csv"))
print(xtable(cld_metric), include.rownames = FALSE)

#################################################################################################

# Rate × Method interaction
cat("\n--- Model × Method Interaction ---\n")
emm_model_method <- emmeans(lmm_model, ~ model * method)
write.csv(cld(emm_model_method, Letters = letters, alpha = 0.05), paste0(lmm_path, "tables/model_method_multcomp.csv"))

# Method × Fossil sampling interaction
cat("\n--- Method × Fossil sampling Interaction ---\n")
emm_method_fossil <- emmeans(lmm_model, ~ method * fossil_sampling)
write.csv(cld(emm_method_fossil, Letters = letters, alpha = 0.05), paste0(lmm_path, "tables/fossil_method_multcomp.csv"))
print(xtable(cld(emm_method_fossil, Letters = letters, alpha = 0.05), include.rownames = FALSE))


# Method x metric
cat("\n--- Method × Metric Interaction ---\n")
emm_method_metric <- emmeans(lmm_model, ~ method * metric)
write.csv(cld(emm_method_metric, Letters = letters, alpha = 0.05), paste0(lmm_path, "tables/method_metric_multcomp.csv"))
print(xtable(cld(emm_method_metric, Letters = letters, alpha = 0.05), include.rownames = FALSE))


#################################################################################################

# Three-way interaction
cat("\n Method x Fossil sampling x model \n")
emm_three_way <- emmeans(lmm_model, ~ method * model * fossil_sampling)
write.csv(cld(emm_three_way, Letters = letters, alpha = 0.05), paste0(lmm_path, "tables/method_model_fossil_multcomp.csv"))
#################################################################################################


# Four-way interaction
cat("\n--- Four-way Interaction ---\n")
emm_four_way <- emmeans(lmm_model, ~ model * method * fossil_sampling * metric)
# Print only a subset to avoid overwhelming output
write.csv(cld(emm_four_way, Letters = letters, alpha = 0.05), paste0(lmm_path, "tables/method_model_fossil_metric_multcomp.csv"))

by_method <- emmeans(lmm_model, ~ method | model * fossil_sampling * metric)
by_method_cld <- cld(by_method, Letters = letters, alpha = 0.05)
write.csv(cld(by_method, Letters = letters, alpha = 0.05), paste0(lmm_path, "tables/method_model_fossil_by_metric_multcomp.csv"))

#################################################################################################

# Compare ACE methods to No ACE
cat("\n Contrast ASE vs no ase \n")
ase_vs_noase <- contrast(emm_method, 
                         list("ASE vs No ASE" = c(1/4, 1/4, 1/4, 1/4, -1)))
print(ase_vs_noase)
print(confint(ase_vs_noase))

# Compare Pre-ord vs Post-ord (averaging point and distribution)
cat("\n Contrast pre-ord vs post-ord \n")
preord_vs_postord <- contrast(emm_method,
                               list("Pre-ord vs Post-ord" = c(1/2, 1/2, -1/2, -1/2, 0)))
print(preord_vs_postord)
print(confint(preord_vs_postord))

# Compare Point vs Sample (averaging pre and post)
cat("\n Contrast point estimate vs probabilistic \n")
point_vs_distribution <- contrast(emm_method,
                             list("Point vs distribution" = c(1/2, -1/2, 1/2, -1/2, 0)))
print(point_vs_distribution)
print(confint(point_vs_distribution))

write.csv(as.data.frame(summary(ase_vs_noase)), paste0(lmm_path, "tables/contrast_ase_vs_noase.csv"))
write.csv(as.data.frame(summary(preord_vs_postord)), paste0(lmm_path, "tables/contrast_preord_vs_postord.csv"))
write.csv(as.data.frame(summary(point_vs_distribution)), paste0(lmm_path, "tables/contrast_point_vs_dist.csv"))


cat("Finished LMM continuos...\n")



df_all <- as.data.frame(by_method_cld)

# 2. Filter for your metric
# Equivalent to: filter(metric == "sum_quant")


winners <- list()
metrics <- c("sum_quant", "sum_var", "pairwise")
for (metric in metrics) {
  df_sub <- df_all[df_all$metric == metric, ]

  df_sub$shared_a <- ave(grepl("a", trimws(df_sub$.group)),
                                           df_sub$model, df_sub$fossil_sampling,
                                           FUN = function(x) sum(x)>1)

  df_sub$min_score <- ave(df_sub$emmean, 
                          df_sub$model, 
                          df_sub$fossil_sampling, 
                          FUN = min)

  # Filter to keep only the Winners
  winners_df <- df_sub[df_sub$emmean == df_sub$min_score, ]
  


  winners[[metric]] <- winners_df
}






















# # saveRDS(lmm_model_final, file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "discrete", "crown","lmm", "lmm_discrete_raw_final.rds"))

# cat("Finished Raw LMM...\n")



# if(args[1] == "rmaxes") {

#     cat("Starting rmaxes LMM...\n")

#     # Load all data with different job IDs
#     rmaxes_results <- list()
#     for (size in tree_sizes) {
#     rmaxes_results[[size]] <- list()
#     job_id <- job_ids[[size]]  #  Get job ID for this tree size
    
#     for (method in methods) {
#         rmaxes_results[[size]][[method]] <- list()
#         for(i in 1:100) {
#         file_path <- file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "discrete", "crown", size,
#                             "disparity", "rm_axes",
#                                 sprintf("%s_%s_%03d.rds", job_id, method, i))
#         if(file.exists(file_path)) {
#             rmaxes_results[[size]][[method]][[i]] <- readRDS(file_path)
#         } else {
#             warning("Missing file: ", file_path)
#             rmaxes_results[[size]][[method]][[i]] <- NULL
#         }
#         }
#     }
#     }

#     metric_names <- names(rmaxes_results[["50t"]][[1]][[1]])  # Get from first available


#     # Transpose structure: metric -> tree_size -> method -> replicate
#     rmaxes_results_relisted <- setNames(lapply(metric_names, function(metric) {
#     setNames(lapply(tree_sizes, function(size) {
#         setNames(lapply(methods, function(method) {
#         is_sample_method <- grepl("sample", method)
        
#         lapply(rmaxes_results[[size]][[method]], function(rep_data) {
#             metric_data <- rep_data[[metric]]
            
#             if(is_sample_method) {
#             lapply(metric_data, function(rate_data) {
#                 lapply(rate_data, function(fossil_data) {
#                 if(is.list(fossil_data) && !is.data.frame(fossil_data)) {
#                     unlist(fossil_data)
#                 } else {
#                     fossil_data
#                 }
#                 })
#             })
#             } else {
#             metric_data
#             }
#         })
#         }), methods)
#     }), tree_sizes)
#     }), metric_names)

#     # Build dataframe with tree_size
#     rmaxes_results_df <- do.call(rbind, lapply(names(rmaxes_results_relisted), function(metric_name) {
#     metric_data <- rmaxes_results_relisted[[metric_name]]
    
#     do.call(rbind, lapply(names(metric_data), function(size_name) {
#         size_data <- metric_data[[size_name]]
        
#         do.call(rbind, lapply(names(size_data), function(method_name) {
#         method_data <- size_data[[method_name]]
        
#         do.call(rbind, lapply(seq_along(method_data), function(rep_idx) {
#             rep_data <- method_data[[rep_idx]]
            
#             if(is.null(rep_data)) return(NULL)
            
#             do.call(rbind, lapply(names(rep_data), function(rate_name) {
#             rate_data <- rep_data[[rate_name]]
            
#             if(is.null(rate_data)) return(NULL)
            
#             data.frame(
#                 replicate = rep_idx,
#                 all = I(list(rate_data$all)),
#                 high = I(list(rate_data$fossil_high)),
#                 mid = I(list(rate_data$fossil_med)),
#                 low = I(list(rate_data$fossil_low)),
#                 living = I(list(rate_data$living)),
#                 rate = rate_name,
#                 method = method_name,
#                 metric = metric_name,
#                 tree_size = size_name,  #  Add tree size
#                 stringsAsFactors = FALSE
#             )
#             }))
#         }))
#         }))
#     }))
#     }))

#     # Pivot to long format
#     rmaxes_results_df_long <- rmaxes_results_df %>%
#     pivot_longer(
#         cols = c("living", "low", "mid", "high", "all"),
#         names_to = "preservation_level",
#         values_to = "error"
#     ) %>%
#     unnest(error)

#     # Factor variables
#     rmaxes_results_df_long$tree_size <- factor(
#     rmaxes_results_df_long$tree_size,
#     levels = c("50t", "100t", "150t"),
#     labels = c("50 taxa", "100 taxa", "150 taxa")
#     )

#     rmaxes_results_df_long$method <- factor(
#     rmaxes_results_df_long$method,
#     levels = c("pre_ord_point_tiebreaker", "pre_ord_sample", 
#                 "post_ord_point", "post_ord_sample", "no_ace"),
#     labels = c("Pre-ord ASE\n(point + tiebreaker)",
#                 "Pre-ord ASE\n(dist)", "Post-ord ASE\n(point)", 
#                 "Post-ord ASE\n(dist)", "No ASE")
#     )

#     rmaxes_results_df_long$preservation_level <- factor(
#     rmaxes_results_df_long$preservation_level,
#     levels = c("living", "low", "mid", "high", "all"),
#     labels = c("0%", "5%", "15%", "50%", "100%")
#     )

#     rmaxes_results_df_long$rate <- factor(
#     rmaxes_results_df_long$rate,
#     levels = c("slow", "med", "fast"),
#     labels = c("Slow", "Medium", "Fast")
#     )

#     rmaxes_results_df_long$metric <- factor(
#     rmaxes_results_df_long$metric,
#     levels = c("sum_var", "sum_quant", "pairwise"),
#     labels = c("Sum of Variances", "Sum of Quantiles", "Mean Pairwise Distance")
#     )

#     rmaxes_results_df_long$abs_error <- abs(rmaxes_results_df_long$error)
#     rmaxes_results_df_long$log_abs_error <- log(rmaxes_results_df_long$abs_error + 0.001)


#     lmm_model <- lmer(log_abs_error ~ rate * method * preservation_level * metric + 
#                         (1|replicate) + (1|tree_size), 
#                         data = rmaxes_results_df_long,
#                         REML = FALSE)


#     saveRDS(lmm_model, file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "discrete", "crown", "lmm", "lmm_discrete_rmaxes.rds"))

# }

# if(args[1] == "combined_raw_rmaxes") {
#     cat("Starting combined raw and rmaxes LMM...\n")

#     # Load all data with different job IDs
#     raw_results <- list()
#     for (size in tree_sizes) {
#     raw_results[[size]] <- list()
#     job_id <- job_ids[[size]]  #  Get job ID for this tree size
    
#     for (method in methods) {
#         raw_results[[size]][[method]] <- list()
#         for(i in 1:100) {
#         file_path <- file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "discrete", "crown", size,
#                             "disparity", "raw",
#                                 sprintf("%s_%s_%03d.rds", job_id, method, i))
#         if(file.exists(file_path)) {
#             raw_results[[size]][[method]][[i]] <- readRDS(file_path)
#         } else {
#             warning("Missing file: ", file_path)
#             raw_results[[size]][[method]][[i]] <- NULL
#         }
#         }
#     }
#     }

#     metric_names <- names(raw_results[["50t"]][[1]][[1]])  # Get from first available


#     # Transpose structure: metric -> tree_size -> method -> replicate
#     raw_results_relisted <- setNames(lapply(metric_names, function(metric) {
#     setNames(lapply(tree_sizes, function(size) {
#         setNames(lapply(methods, function(method) {
#         is_sample_method <- grepl("sample", method)
        
#         lapply(raw_results[[size]][[method]], function(rep_data) {
#             metric_data <- rep_data[[metric]]
            
#             if(is_sample_method) {
#             lapply(metric_data, function(rate_data) {
#                 lapply(rate_data, function(fossil_data) {
#                 if(is.list(fossil_data) && !is.data.frame(fossil_data)) {
#                     unlist(fossil_data)
#                 } else {
#                     fossil_data
#                 }
#                 })
#             })
#             } else {
#             metric_data
#             }
#         })
#         }), methods)
#     }), tree_sizes)
#     }), metric_names)

#     rm(raw_results)
#     gc()


#     # Build dataframe with tree_size
#     raw_results_df <- do.call(rbind, lapply(names(raw_results_relisted), function(metric_name) {
#     metric_data <- raw_results_relisted[[metric_name]]
    
#     do.call(rbind, lapply(names(metric_data), function(size_name) {
#         size_data <- metric_data[[size_name]]
        
#         do.call(rbind, lapply(names(size_data), function(method_name) {
#         method_data <- size_data[[method_name]]
        
#         do.call(rbind, lapply(seq_along(method_data), function(rep_idx) {
#             rep_data <- method_data[[rep_idx]]
            
#             if(is.null(rep_data)) return(NULL)
            
#             do.call(rbind, lapply(names(rep_data), function(rate_name) {
#             rate_data <- rep_data[[rate_name]]
            
#             if(is.null(rate_data)) return(NULL)
            
#             data.frame(
#                 replicate = rep_idx,
#                 all = I(list(rate_data$all)),
#                 high = I(list(rate_data$fossil_high)),
#                 mid = I(list(rate_data$fossil_med)),
#                 low = I(list(rate_data$fossil_low)),
#                 living = I(list(rate_data$living)),
#                 rate = rate_name,
#                 method = method_name,
#                 metric = metric_name,
#                 tree_size = size_name,  #  Add tree size
#                 stringsAsFactors = FALSE
#             )
#             }))
#         }))
#         }))
#     }))
#     }))

#     rm(raw_results_relisted)
#     gc()

#     # Pivot to long format
#     raw_results_df_long <- raw_results_df %>%
#     pivot_longer(
#         cols = c("living", "low", "mid", "high", "all"),
#         names_to = "preservation_level",
#         values_to = "error"
#     ) %>%
#     unnest(error)

#     rm(raw_results_df)
#     gc()

#     # Factor variables
#     raw_results_df_long$tree_size <- factor(
#     raw_results_df_long$tree_size,
#     levels = c("50t", "100t", "150t"),
#     labels = c("50 taxa", "100 taxa", "150 taxa")
#     )

#     raw_results_df_long$method <- factor(
#     raw_results_df_long$method,
#     levels = c("pre_ord_point_tiebreaker", "pre_ord_sample", 
#                 "post_ord_point", "post_ord_sample", "no_ace"),
#     labels = c("Pre-ord ASE\n(point + tiebreaker)",
#                 "Pre-ord ASE\n(dist)", "Post-ord ASE\n(point)", 
#                 "Post-ord ASE\n(dist)", "No ASE")
#     )

#     raw_results_df_long$preservation_level <- factor(
#     raw_results_df_long$preservation_level,
#     levels = c("living", "low", "mid", "high", "all"),
#     labels = c("0%", "5%", "15%", "50%", "100%")
#     )

#     raw_results_df_long$rate <- factor(
#     raw_results_df_long$rate,
#     levels = c("slow", "med", "fast"),
#     labels = c("Slow", "Medium", "Fast")
#     )

#     raw_results_df_long$metric <- factor(
#     raw_results_df_long$metric,
#     levels = c("sum_var", "sum_quant", "pairwise"),
#     labels = c("Sum of Variances", "Sum of Quantiles", "Mean Pairwise Distance")
#     )

#     # Load all data with different job IDs
#     rmaxes_results <- list()
#     for (size in tree_sizes) {
#     rmaxes_results[[size]] <- list()
#     job_id <- job_ids[[size]]  #  Get job ID for this tree size
    
#     for (method in methods) {
#         rmaxes_results[[size]][[method]] <- list()
#         for(i in 1:100) {
#         file_path <- file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "discrete", "crown", size,
#                             "disparity", "rm_axes",
#                                 sprintf("%s_%s_%03d.rds", job_id, method, i))
#         if(file.exists(file_path)) {
#             rmaxes_results[[size]][[method]][[i]] <- readRDS(file_path)
#         } else {
#             warning("Missing file: ", file_path)
#             rmaxes_results[[size]][[method]][[i]] <- NULL
#         }
#         }
#     }
#     }

#     metric_names <- names(rmaxes_results[["50t"]][[1]][[1]])  # Get from first available


#     # Transpose structure: metric -> tree_size -> method -> replicate
#     rmaxes_results_relisted <- setNames(lapply(metric_names, function(metric) {
#     setNames(lapply(tree_sizes, function(size) {
#         setNames(lapply(methods, function(method) {
#         is_sample_method <- grepl("sample", method)
        
#         lapply(rmaxes_results[[size]][[method]], function(rep_data) {
#             metric_data <- rep_data[[metric]]
            
#             if(is_sample_method) {
#             lapply(metric_data, function(rate_data) {
#                 lapply(rate_data, function(fossil_data) {
#                 if(is.list(fossil_data) && !is.data.frame(fossil_data)) {
#                     unlist(fossil_data)
#                 } else {
#                     fossil_data
#                 }
#                 })
#             })
#             } else {
#             metric_data
#             }
#         })
#         }), methods)
#     }), tree_sizes)
#     }), metric_names)

#     rm(rmaxes_results)
#     gc()

#     # Build dataframe with tree_size
#     rmaxes_results_df <- do.call(rbind, lapply(names(rmaxes_results_relisted), function(metric_name) {
#     metric_data <- rmaxes_results_relisted[[metric_name]]
    
#     do.call(rbind, lapply(names(metric_data), function(size_name) {
#         size_data <- metric_data[[size_name]]
        
#         do.call(rbind, lapply(names(size_data), function(method_name) {
#         method_data <- size_data[[method_name]]
        
#         do.call(rbind, lapply(seq_along(method_data), function(rep_idx) {
#             rep_data <- method_data[[rep_idx]]
            
#             if(is.null(rep_data)) return(NULL)
            
#             do.call(rbind, lapply(names(rep_data), function(rate_name) {
#             rate_data <- rep_data[[rate_name]]
            
#             if(is.null(rate_data)) return(NULL)
            
#             data.frame(
#                 replicate = rep_idx,
#                 all = I(list(rate_data$all)),
#                 high = I(list(rate_data$fossil_high)),
#                 mid = I(list(rate_data$fossil_med)),
#                 low = I(list(rate_data$fossil_low)),
#                 living = I(list(rate_data$living)),
#                 rate = rate_name,
#                 method = method_name,
#                 metric = metric_name,
#                 tree_size = size_name,  #  Add tree size
#                 stringsAsFactors = FALSE
#             )
#             }))
#         }))
#         }))
#     }))
#     }))

#     rm(rmaxes_results_relisted)
#     gc()

#     # Pivot to long format
#     rmaxes_results_df_long <- rmaxes_results_df %>%
#     pivot_longer(
#         cols = c("living", "low", "mid", "high", "all"),
#         names_to = "preservation_level",
#         values_to = "error"
#     ) %>%
#     unnest(error)

#     rm(rmaxes_results_df)
#     gc()

#     # Factor variables
#     rmaxes_results_df_long$tree_size <- factor(
#     rmaxes_results_df_long$tree_size,
#     levels = c("50t", "100t", "150t"),
#     labels = c("50 taxa", "100 taxa", "150 taxa")
#     )

#     rmaxes_results_df_long$method <- factor(
#     rmaxes_results_df_long$method,
#     levels = c("pre_ord_point_tiebreaker",  "pre_ord_sample", 
#                 "post_ord_point", "post_ord_sample", "no_ace"),
#     labels = c("Pre-ord ASE\n(point + tiebreaker)",
#                 "Pre-ord ASE\n(dist)", "Post-ord ASE\n(point)", 
#                 "Post-ord ASE\n(dist)", "No ASE")
#     )

#     rmaxes_results_df_long$preservation_level <- factor(
#     rmaxes_results_df_long$preservation_level,
#     levels = c("living", "low", "mid", "high", "all"),
#     labels = c("0%", "5%", "15%", "50%", "100%")
#     )

#     rmaxes_results_df_long$rate <- factor(
#     rmaxes_results_df_long$rate,
#     levels = c("slow", "med", "fast"),
#     labels = c("Slow", "Medium", "Fast")
#     )

#     rmaxes_results_df_long$metric <- factor(
#     rmaxes_results_df_long$metric,
#     levels = c("sum_var", "sum_quant", "pairwise"),
#     labels = c("Sum of Variances", "Sum of Quantiles", "Mean Pairwise Distance")
#     )

#     raw_results_df_long$axes_treatment <- "Raw"
#     rmaxes_results_df_long$axes_treatment <- "Axes Removed"

#     # Now combine
#     combined_rm_raw_df <- rbind(raw_results_df_long, rmaxes_results_df_long)

#     rm(raw_results_df_long)
#     rm(rmaxes_results_df_long)
#     gc()

#     # Optionally make it a factor with specific order
#     combined_rm_raw_df$axes_treatment <- factor(
#     combined_rm_raw_df$axes_treatment,
#     levels = c("Raw", "Axes Removed")
#     )

#     combined_rm_raw_df$abs_error <- abs(combined_rm_raw_df$error)
#     combined_rm_raw_df$log_abs_error <- log(combined_rm_raw_df$abs_error + 0.001)


#     lmm_model <- lmer(log_abs_error ~ rate * method * preservation_level * metric + 
#                         axes_treatment +
#                         axes_treatment:method + axes_treatment:preservation_level +
#                         (1|replicate) + (1|tree_size), 
#                         data = combined_rm_raw_df,
#                         REML = FALSE)


#     saveRDS(lmm_model, file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "discrete", "crown", "lmm", "lmm_discrete_combined_raw_rmaxes.rds"))
#     cat("Finished combined raw and rmaxes LMM...\n")

# }



