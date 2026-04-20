cat("Starting R script...\n")


cat("Loading packages...\n")
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest) 
library(emmeans)
library(multcomp)
library(multcompView)
library(xtable)
library(ggplot2)



    
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
      file_path <- file.path("/mnt/parscratch/users/bip24cns/acedisparity/discrete/crown", size, 
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



counts <- ave(
  rep(1, nrow(results_df_long)),
  results_df_long$tree_size,            # A vector of 1s to count
  results_df_long$replicate,            # Grouping var 1
  results_df_long$method,
  results_df_long$metric,
  results_df_long$fossil_sampling, 
  results_df_long$model,                  # Grouping var 2
  FUN = length                   # Function to count rows
)

results_df_long$weight <- 1 / counts

checks <- aggregate(weight ~ replicate + method + tree_size + fossil_sampling + model + metric, data = results_df_long, FUN = sum) ## should all equal 1


if(!any(is.na(checks$weight == 1))){cat("All weights add up to one...\n")}

rm(results_relisted)
rm(results_df)

gc()

# Combine back together
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


lmm_model_weighted <- lmer(log_abs_error ~ model * method * fossil_sampling * metric + tree_size +
                                  (1|tree_unique_id), weights = weight,
                    data = results_df_long, REML = TRUE)




lmm_path <- "/mnt/parscratch/users/bip24cns/acedisparity/discrete/crown/lmm/"

saveRDS(lmm_model_weighted, paste0(lmm_path, "lmm_4_way_weighted.rds"))

# Residual plots
pdf(paste0(lmm_path, "diagnostics/diagnostic_plot.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2))

# QQ plot
qqnorm(resid(lmm_model_weighted), main = "Q-Q Plot of Residuals")
qqline(resid(lmm_model_weighted), col = "red")

# Residuals vs Fitted
plot(fitted(lmm_model_weighted), resid(lmm_model_weighted),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# Scale-Location
plot(fitted(lmm_model_weighted), sqrt(abs(resid(lmm_model_weighted))),
     xlab = "Fitted Values", ylab = "Square Root of |Residuals|",
     main = "Scale-Location")

# Histogram of residuals
hist(resid(lmm_model_weighted), breaks = 50, main = "Histogram of residuals",
     xlab = "Residuals")
dev.off()


sink(paste0(lmm_path, "model_summaries_four_way.txt"))
cat("\n Model Summary \n")
print(summary(lmm_model_weighted))

cat("\n Model fit \n")
cat("AIC:", AIC(lmm_model_weighted), "\n")
cat("BIC:", BIC(lmm_model_weighted), "\n")
cat("Log-Likelihood:", logLik(lmm_model_weighted), "\n")

cat("\n Random effects variance \n")
print(VarCorr(lmm_model_weighted))

# ANOVA table type 3 due to unequal sample sizes
cat("\n ANOVA table type 3 due to unequal sample sizes \n")
print(anova(lmm_model_weighted, type = 3))
sink()

#################################################################################################

# Estimated marginal means and Tukey comparisons
cat("\n--- Rate Rankings ---\n")
emm_model <- emmeans(lmm_model_weighted, ~ model)
cld_model <- cld(emm_model, Letters = letters, alpha = 0.05)
write.csv(cld_model, paste0(lmm_path, "tables/model_multcomp.csv"))
print(xtable(cld_model), 
      include.rownames = FALSE)

cat("\n--- Method Rankings ---\n")
emm_method <- emmeans(lmm_model_weighted, ~ method)
cld_method <- cld(emm_method, Letters = letters, alpha = 0.05)
write.csv(cld_method, paste0(lmm_path, "tables/method_multcomp.csv"))
print(xtable(cld_method), 
      include.rownames = FALSE)


cat("\n--- Fossil sampling Level Rankings ---\n")
emm_fossils <- emmeans(lmm_model_weighted, ~ fossil_sampling)
cld_fossils <- cld(emm_fossils, Letters = letters, alpha = 0.05)
write.csv(cld_fossils, paste0(lmm_path, "tables/fossil_multcomp.csv"))
print(xtable(cld_fossils), include.rownames = FALSE)

cat("\n--- Metric Level Rankings ---\n")
emm_metric <- emmeans(lmm_model_weighted, ~ metric)
cld_metric <- cld(emm_metric, Letters = letters, alpha = 0.05)
write.csv(cld_metric, paste0(lmm_path, "tables/metric_multcomp.csv"))
print(xtable(cld_metric), include.rownames = FALSE)

#################################################################################################

# Rate × Method interaction
cat("\n--- Model × Method Interaction ---\n")
emm_model_method <- emmeans(lmm_model_weighted, ~ model * method)
write.csv(cld(emm_model_method, Letters = letters, alpha = 0.05), paste0(lmm_path, "tables/model_method_multcomp.csv"))

# Method × Fossil sampling interaction
cat("\n--- Method × Fossil sampling Interaction ---\n")
emm_method_fossil <- emmeans(lmm_model_weighted, ~ method * fossil_sampling)
write.csv(cld(emm_method_fossil, Letters = letters, alpha = 0.05), paste0(lmm_path, "tables/fossil_method_multcomp.csv"))
print(xtable(cld(emm_method_fossil, Letters = letters, alpha = 0.05), include.rownames = FALSE))


# Method x metric
cat("\n--- Method × Metric Interaction ---\n")
emm_method_metric <- emmeans(lmm_model_weighted, ~ method * metric)
write.csv(cld(emm_method_metric, Letters = letters, alpha = 0.05), paste0(lmm_path, "tables/method_metric_multcomp.csv"))
print(xtable(cld(emm_method_metric, Letters = letters, alpha = 0.05), include.rownames = FALSE))


#################################################################################################

# Three-way interaction
cat("\n Method x Fossil sampling x model \n")
emm_three_way <- emmeans(lmm_model_weighted, ~ method * model * fossil_sampling)
saveRDS(emm_three_way, paste0(lmm_path, "three_way_emm.rds"))

write.csv(cld(emm_three_way, Letters = letters, alpha = 0.05), paste0(lmm_path, "tables/method_model_fossil_multcomp.csv"))
print(xtable(cld(emm_three_way, Letters = letters, alpha = 0.05), include.rownames = FALSE))

#################################################################################################


# Four-way interaction
cat("\n--- Four-way Interaction ---\n")
emm_four_way <- emmeans(lmm_model_weighted, ~ model * method * fossil_sampling * metric)
# Print only a subset to avoid overwhelming output
write.csv(cld(emm_four_way, Letters = letters, alpha = 0.05), paste0(lmm_path, "tables/method_model_fossil_metric_multcomp.csv"))

by_method <- emmeans(lmm_model_weighted, ~ method | model * fossil_sampling * metric)
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

