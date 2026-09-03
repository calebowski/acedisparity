library(ggplot2)
library(tidyr)
library(patchwork)
library(scales)
# Read and rearrange in one step
rate_names <- c("slow", "med", "fast")


methods <- c("pre_ord_point", "pre_ord_sample", "post_ord_point", 
             "post_ord_sample", "no_ace")

# Load all data with different job IDs
results <- list()


for (method in methods) {
  results[[method]] <- list()
  for(i in 1:100) {
    file_path <- file.path("..", "Data",  "revisions", "wagner", "discrete",
                          "disparity",
                          sprintf("11429393_%s_%03d.rds", method, i))
    if(file.exists(file_path)) {
      results[[method]][[i]] <- readRDS(file_path)
    } else {
      warning("Missing file: ", file_path)
      results[[method]][[i]] <- NULL
    }
  }
}


metric_names <- names(results[[1]][[1]])  # Get from first available

# Transpose structure: metric -> tree_size -> method -> replicate
results_relisted <- setNames(lapply(metric_names, function(metric) {
    setNames(lapply(methods, function(method) {
      is_sample_method <- grepl("sample", method)
      
      lapply(results[[method]], function(rep_data) {
        metric_data <- rep_data[[metric]]
        
        if(is_sample_method) {
          lapply(metric_data, function(rate_data) {
            lapply(rate_data, function(fossil_data) {
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
}), metric_names)

# Build dataframe with tree_size
results_df <- do.call(rbind, lapply(names(results_relisted), function(metric_name) {
  metric_data <- results_relisted[[metric_name]]
  do.call(rbind, lapply(names(metric_data), function(method_name) {
    method_data <- metric_data[[method_name]]
    
    do.call(rbind, lapply(seq_along(method_data), function(rep_idx) {
      rep_data <- method_data[[rep_idx]]
      
      if(is.null(rep_data)) return(NULL)
      
      do.call(rbind, lapply(names(rep_data), function(rate_name) {
        rate_data <- rep_data[[rate_name]]
        
        if(is.null(rate_data)) return(NULL)
        
        data.frame(
          replicate = rep_idx,
          high = I(list(rate_data$fossil_high)),
          mid = I(list(rate_data$fossil_med)),
          low = I(list(rate_data$fossil_low)),
          living = I(list(rate_data$living)),
          rate = rate_name,
          method = method_name,
          metric = metric_name,
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))
  }))

# Pivot to long format
results_df_long <- results_df %>%
  pivot_longer(
    cols = c("living", "low", "mid", "high"),
    names_to = "preservation_level",
    values_to = "error"
  ) %>%
  unnest(error)


results_df_long$method <- factor(
  results_df_long$method,
  levels = c("pre_ord_point", "pre_ord_sample", 
             "post_ord_point", "post_ord_sample", "no_ace"),
  labels = c("Pre-ord ASE\n(point)",
             "Pre-ord ASE\n(dist)", "Post-ord ASE\n(point)", 
             "Post-ord ASE\n(dist)", "Sampled taxa only")
)

results_df_long$preservation_level <- factor(
  results_df_long$preservation_level,
  levels = c("living", "low", "mid", "high"),
  labels = c("0", "5", "15", "50")
)

results_df_long$rate <- factor(
  results_df_long$rate,
  levels = c("slow", "med", "fast"),
  labels = c("Slow", "Medium", "Fast")
)

results_df_long$metric <- factor(
  results_df_long$metric,
  levels = c("sum_var", "sum_quant", "pairwise"),
  labels = c("Sum of Variances", "Sum of Quantiles", "Mean Pairwise Distance")
)

## do boxplot

colours <- c("Slow" = "#8805A8", "Medium" = "#00B945", "Fast" = "#FFB600")

boxplot_plot_all <- ggplot(results_df_long, 
                           aes(x = preservation_level, y = error, fill = rate, colour = rate)) +
  geom_boxplot(
    alpha = 0.7,
    # outlier.alpha = 0.1,
    outlier.shape = NA,  # Hide outliers for cleaner plot
    position = position_dodge(width = 0.8)
  ) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.8, linetype = "dashed") +
  facet_grid(metric ~ method) +  #  Only facet by metric and method
  labs(
    x = "Fossil Sampling (%)",
    y = "Relative Disparity Error",
    fill = "Transition Rate"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = colours,
    labels = c("Slow", "Medium", "Fast")
  ) +
  scale_colour_manual(
    values = colours,
    labels = c("Slow", "Medium", "Fast"),
    guide = "none"
  ) +
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15),
    strip.text = element_text(size = 20, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.3, "cm")
  ) +
  scale_y_continuous(
    breaks = c(-1, -0.5, 0, 0.50, 1.0),
    labels = c("-1.0", "-0.5", "0", "0.50", "1.0")
  ) +
  coord_cartesian(ylim = c(-1, 1))

ggsave("../../Manuscript/draft/figures/discrete_boxplot_nodes_tips_sampled.png", boxplot_plot_all, , width = 18, height = 14, dpi = 700, units = "in", bg = "white")
