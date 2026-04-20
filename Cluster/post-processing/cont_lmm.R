library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest) 
library(emmeans)
library(multcomp)
library(multcompView)
library(xtable)
library(tatoo)
library(partykit)

job_ids <- list(
  "50t" = "8690335",
  "100t" = "8690474",  # Replace with actual job ID
  "150t" = "8690995"   # Replace with actual job ID
)

tree_sizes <- c("50t", "100t", "150t")
methods <- c("pre_ord_point", "pre_ord_sample", "post_ord_point", "post_ord_sample", "no_ace")

# # Load all data with different job IDs
# results <- list()
# for (size in tree_sizes) {
#   results[[size]] <- list()
#   job_id <- job_ids[[size]]  #  Get job ID for this tree size
  
#   for (method in methods) {
#     results[[size]][[method]] <- list()
#     for(i in 1:100) {
#       file_path <- file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "continuous", size,
#                             "disparity",
#                             sprintf("%s_%s_%03d.rds", job_id, method, i))
#       if(file.exists(file_path)) {
#         results[[size]][[method]][[i]] <- readRDS(file_path)
#       } else {
#         warning("Missing file: ", file_path)
#         results[[size]][[method]][[i]] <- NULL
#       }
#     }
#   }
# }

results <- list()
for (size in tree_sizes) {
  results[[size]] <- list()
  job_id <- job_ids[[size]]  #  Get job ID for this tree size
  
  for (method in methods) {
    results[[size]][[method]] <- list()
    for(i in 1:100) {
      file_path <- file.path("..", "Data", "continuous", size, 
                            "disparity",
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

rm(results_relisted)
gc()


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




results_df_long_filtered <- results_df_long[results_df_long$model %in% c("bm", "bm_t", "ou_st"),]
table(results_df_long_filtered$method)

results_df_long_filtered$tree_unique_id <- interaction(results_df_long_filtered$tree_size, 
                                              results_df_long_filtered$replicate)

# 2. Unique Simulation ID (Accounts for the specific trait generation)
# results_df_long_filtered$sim_unique_id <- interaction(results_df_long_filtered$tree_unique_id, 
#                                              results_df_long_filtered$model)

rm(results_df_long)
gc()
# Transform for normality (log of absolute values)
results_df_long_filtered$abs_error <- abs(results_df_long_filtered$error)
results_df_long_filtered$log_abs_error <- log(results_df_long_filtered$abs_error + 0.001)



#  What questions do we want to know: 
#   - Fossil sampling x method
#   - Metric x method
#   - Model x method
#   - FOssil sampling x metric x method x model
  
#  Linear Mixed Model
lmm_model <- lmer(log_abs_error ~ model * method * fossil_sampling * metric + tree_size +
                                  (1|tree_unique_id),
                    data = results_df_long_filtered, REML = TRUE)


# saveRDS(lmm_model, file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "continuous","lmm", "lmm_continuous_3_model.rds"))

# filtered_model <- lmerTest::step(lmm_model) ## backwards step selection to find best model

# got_model <- lmerTest::get_model(filtered_model)
# final_formula <- formula(got_model)

# lmm_model_final <- lmer(final_formula, data = results_df_long_filtered)

saveRDS(lmm_model, "../Data/cluster/continuous/lmm/lmm_model_four_way.rds")
lmm_model <- readRDS("../Data/continuous/lmm/lmm_model_four_way.rds")


#################################################################################################
lmm_path <- "/home/caleb/Documents/PhD/acedisparity/Data/cluster/continuous/lmm/"

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
print(xtable((anova(lmm_model, type = 3))))
print((anova(lmm_model, type = 3)))

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

# Model × Method interaction
cat("\n--- Model × Method Interaction ---\n")
emm_model_method <- emmeans(lmm_model, ~ model * method)
write.csv(cld(emm_model_method, Letters = letters, alpha = 0.05), paste0(lmm_path, "tables/model_method_multcomp.csv"))
print(xtable(cld(emm_model_method, Letters = letters, alpha = 0.05)))

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
print(xtable(cld(emm_three_way, Letters = letters, alpha = 0.05), include.rownames = FALSE))

#################################################################################################


# Four-way interaction
cat("\n--- Four-way Interaction ---\n")
emm_four_way <- emmeans(lmm_model, ~ model * method * fossil_sampling * metric)
# Print only a subset to avoid overwhelming output
write.csv(cld(emm_four_way, Letters = letters, alpha = 0.05), paste0(lmm_path, "tables/method_model_fossil_metric_multcomp.csv"))

by_method <- emmeans(lmm_model, ~ method | model * fossil_sampling * metric)
write.csv(cld(by_method, Letters = letters, alpha = 0.05), paste0(lmm_path, "tables/method_model_fossil_by_metric_multcomp.csv"))
by_method_cld <- cld(by_method, Letters = letters, alpha = 0.05)
print(xtable(by_method_cld, include.rownames = FALSE))

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


#################################################################################################

fossil_method_model <- read.csv("../Data/cluster/continuous/lmm/tables/method_model_fossil_metric_multcomp.csv")


chunks <- split(fossil_method_model, list(fossil_method_model$model, fossil_method_model$fossil_sampling), drop = TRUE) ## visually inspect it


fossil_method_model_metric <- read.csv("../Data/cluster/continuous/lmm/tables/method_model_fossil_by_metric_multcomp.csv")
chunks <- split(fossil_method_model_metric, list(fossil_method_model_metric$metric, fossil_method_model_metric$model, fossil_method_model_metric$fossil_sampling), drop = TRUE) ## visually inspect it


chunks_nested <- lapply(split(fossil_method_model_metric, fossil_method_model_metric$metric), function(metric_data) {
  lapply(split(metric_data, metric_data$model), function(model_data) {
    lapply(split(model_data, model_data$fossil_sampling), function(fossil_data) {
      fossil_data
    })
  })
})


#########################################################################################################
## HEATMAP
library(ggplot2)
library(dplyr)

# 1. Prepare the Data
# Convert emmeans object to a standard dataframe
# 1. Prepare the data
plot_df <- as.data.frame(emm_three_way)
emm_three_way <- emmeans(lmm_model, ~ method * model * fossil_sampling, type = "response")
plot_df <- as.data.frame(emm_three_way)


# 2. Clean up factors
plot_df$fossil_sampling <- factor(plot_df$fossil_sampling, 
                                   levels = c("living", "low", "med", "high", "all"),labels = c("0%", "5%", "15%", "50%", "100%"))

plot_df$method <- factor(plot_df$method, 
                                   levels = c("pre_ord_point", "pre_ord_sample", "post_ord_point", "post_ord_sample", "no_ace"),
                                   labels = c("Pre-ord point", "Pre-ord dist", "Post-ord point", "Post-ord dist", "No ASE"))

plot_df$model <- factor(plot_df$model, 
                                   levels = c("bm", "bm_t", "ou_st" ),
                                   labels = c("BM", "BM + trend", "OU"))

plot_df$e_emmean <- exp(plot_df$emmean)

my_breaks <- quantile(plot_df$e_emmean, probs = seq(0, 1, 0.25))


p <- ggplot(plot_df, aes(x = fossil_sampling, y = method, fill = e_emmean)) +
    geom_tile(color = "white", linewidth = 1) +
    # geom_text(aes(label = sprintf("%.2f", emmean)), color = "white", size = 3) +
    facet_wrap(~ model, nrow = 1) +
    scale_fill_viridis_c(option = "rocket", direction = -1, 
                         name = "EMM Relative\nDisparity Error", values = scales::rescale(my_breaks)) +
    # scale_fill_distiller("RdYlBu", direction = +1) +
    # labs(title = paste("Metric:", metric_name)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 14),
      strip.text = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 20, face = "bold", color = "black"),
      ) +
      labs(x = "Fossil Sampling", y = "Estimation method") +
      coord_equal() +
      scale_y_discrete(limits = rev)

ggsave("../../Manuscript/draft/figures/cont_heatmap.pdf", p, "pdf", width = 10, height = 4, units = "in")






# ############################################################################################
# # 1. Prepare Data
# df_all <- as.data.frame(by_method_cld)

# # 2. Filter for your metric
# # Equivalent to: filter(metric == "sum_quant")


# winners <- list()
# metrics <- c("sum_quant", "sum_var", "pairwise")
# for (metric in metrics) {
#   df_sub <- df_all[df_all$metric == metric, ]

#   df_sub$shared_a <- ave(grepl("a", trimws(df_sub$.group)),
#                                            df_sub$model, df_sub$fossil_sampling,
#                                            FUN = function(x) sum(x)>1)

#   df_sub$min_score <- ave(df_sub$emmean, 
#                           df_sub$model, 
#                           df_sub$fossil_sampling, 
#                           FUN = min)

#   # Filter to keep only the Winners
#   winners_df <- df_sub[df_sub$emmean == df_sub$min_score, ]
  


#   winners[[metric]] <- winners_df
# }

# # --- SETUP PALETTE ---
# # Define colors for your methods manually (Clean, professional look)
# # Adjust names to match your exact method names
# method_colors <- c("pre_ord_sample" = "#377EB8",  # Blue
#                    "post_ord_sample" = "#4DAF4A", # Green
#                    "no_ace" = "#E41A1C",          # Red
#                    "pre_ord_point" = "#984EA3",   # Purple
#                    "post_ord_point" = "#FF7F00")  # Orange

# # Map methods to colors in the dataframe
# winners_df$color <- method_colors[as.character(winners_df$method)]
# # Fallback to grey if a method name doesn't match the list above
# winners_df$color[is.na(winners_df$color)] <- "grey"

# # --- DRAW PLOT ---
# # 1. Setup Canvas (type="n" means 'none', just setup axes)
# par(mar = c(5, 5, 4, 7)) # Increase right margin for legend
# plot(x = range(winners_df$x_coord), 
#      y = range(winners_df$y_coord), 
#      type = "n", 
#      axes = FALSE, 
#      xlab = "Fossil Availability", 
#      ylab = "Evolutionary Model",
#      main = "Optimal Ancestral Reconstruction Strategy",
#      xlim = c(0.5, 5.5), 
#      ylim = c(0.5, length(unique(winners_df$model)) + 0.5))

# # 2. Draw Tiles (rectangles)
# # We draw a box centered on the x,y coordinate
# rect(xleft = winners_df$x_coord - 0.5, 
#      ybottom = winners_df$y_coord - 0.5, 
#      xright = winners_df$x_coord + 0.5, 
#      ytop = winners_df$y_coord + 0.5, 
#      col = winners_df$color, 
#      border = "white", 
#      lwd = 2) # White grid lines

# # 3. Add Text Labels (Method Names)
# text(x = winners_df$x_coord, 
#      y = winners_df$y_coord, 
#      labels = winners_df$method, 
#      col = "white", 
#      font = 2, # Bold
#      cex = 0.8) # Slightly smaller text

# # 4. Add "Statistical Tie" Points
# # Check if .group is NOT "a" (using trimws to remove spaces like "a ")
# ties <- winners_df[trimws(winners_df$.group) != "a", ]

# if(nrow(ties) > 0) {
#   points(x = ties$x_coord, 
#          y = ties$y_coord, 
#          pch = 21, # Circle with outline
#          bg = "white", 
#          col = "black", 
#          cex = 1.5)
# }

# # 5. Custom Axes
# ggplot(plot_df, aes(x = fossil_sampling, y = emmean, color = method, group = method)) +
#   # Add error bars (assuming you have SE)
#   geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
#   # Add points and lines
#   geom_point(size = 2) +
#   geom_line(size = 1) +
#   # Facet by Metric (rows) and Model (cols)
#   facet_grid(metric ~ model) +
#   # Flip the axis if 'emmean' is error (so lower is visually lower)
#   # If negative values are 'good', keep as is. If positive values are error:
#   scale_y_continuous(trans = "reverse") + 
#   labs(
#     y = "Log Disparity Error (Lower is Better)",
#     x = "Fossil Sampling Intensity",
#     title = "Method Performance across Evolutionary Contexts"
#   ) +
#   theme_bw()
# ###########################################################################
# axis(1, at = 1:5, labels = levels(winners_df$fossil_sampling), tick = FALSE)
# axis(2, at = 1:length(levels(winners_df$model)), labels = levels(winners_df$model), tick = FALSE, las = 2)

# # 6. Add Legend
# legend("right", 
#        legend = c("Stat. Tie"), 
#        pch = 21, 
#        pt.bg = "white",
#        col = "black",
#        bty = "n",
#        inset = c(-0.15, 0), # Move outside plot
#        xpd = TRUE)



## weighted

lmm_path <- "/home/caleb/Documents/PhD/acedisparity/Cluster/Data/continuous/lmm/weighted/"

lmm_four_way <- readRDS(paste0(lmm_path, "lmm_4_way_weighted.rds"))
summary(lmm_four_way,ddf = "Satterthwaite")
emm_three_way <- readRDS(paste0(lmm_path, "three_way_emm.rds"))


print(anova(lmm_four_way, type = 3))

library(ggplot2)
library(dplyr)

# 1. Prepare the Data
# Convert emmeans object to a standard dataframe
# 1. Prepare the data
plot_df <- as.data.frame(emm_three_way)

# 2. Clean up factors
plot_df$fossil_sampling <- factor(plot_df$fossil_sampling, 
                                   levels = c("living", "low", "med", "high", "all"),labels = c("0%", "5%", "15%", "50%", "100%"))

plot_df$method <- factor(plot_df$method, 
                                   levels = c("pre_ord_point", "pre_ord_sample", "post_ord_point", "post_ord_sample", "no_ace"),
                                   labels = c("Pre-ord point", "Pre-ord dist", "Post-ord point", "Post-ord dist", "No ASE"))

plot_df$model <- factor(plot_df$model, 
                                   levels = c("bm", "bm_t", "ou_st" ),
                                   labels = c("BM", "BM + trend", "OU"))

plot_df$e_emmean <- exp(plot_df$emmean)


my_breaks <- quantile(plot_df$e_emmean, probs = seq(0, 1, 0.25))

p <- ggplot(plot_df, aes(x = fossil_sampling, y = method, fill = e_emmean)) +
    geom_tile(color = "white", linewidth = 1) +
    # geom_text(aes(label = sprintf("%.2f", emmean)), color = "white", size = 3) +
    facet_wrap(~ model, nrow = 1) +
    scale_fill_viridis_c(option = "rocket", direction = -1, 
                         name = "Emmean Relative\nDisparity Error", values = scales::rescale(my_breaks)) +
    # scale_fill_distiller("RdYlBu", direction = +1) +
    # labs(title = paste("Metric:", metric_name)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 14),
      strip.text = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 20, face = "bold", color = "black"),
      ) +
      labs(x = "Fossil Sampling", y = "Estimation method") +
      coord_equal() +
      scale_y_discrete(limits = rev)

ggsave("../Manuscript/draft/figures/cont_heatmap_weighted.pdf", p, "pdf", width = 10, height = 4, units = "in")
