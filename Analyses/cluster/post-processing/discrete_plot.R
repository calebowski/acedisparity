## read in the disp diffs


# Read and rearrange in one step
rate_names <- c("slow", "med", "fast")


# r <- readRDS("../Data/cluster/discrete/raw/diff_sumvar_strict_001.rds")

sum_var_diff_strict <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/diff_sumvar_strict_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

sum_var_diff_sample <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/diff_sumvar_sample_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

sum_var_diff_sample_median <- lapply(sum_var_diff_sample, lapply, lapply, median)

sum_var_diff_rel <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/diff_sumvar_rel_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

sum_var_diff_no_ace <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/diff_sumvar_no_ace_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

sum_var_diff_point_postord <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/diff_sumvar_point_postord_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)


sum_var_diff_sample_postord <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/diff_sumvar_sample_postord_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)


saveRDS(sum_var_diff_strict, "../Data/cluster/discrete/dispdiffs/sumvar_diff_strict.rds")
saveRDS(sum_var_diff_sample, "../Data/cluster/discrete/dispdiffs/sumvar_diff_sample.rds")
saveRDS(sum_var_diff_rel, "../Data/cluster/discrete/dispdiffs/sumvar_diff_rel.rds")
saveRDS(sum_var_diff_no_ace, "../Data/cluster/discrete/dispdiffs/sumvar_diff_no_ace.rds")
saveRDS(sum_var_diff_point_postord, "../Data/cluster/discrete/dispdiffs/sumvar_diff_point_postord.rds")
saveRDS(sum_var_diff_sample_postord, "../Data/cluster/discrete/dispdiffs/sumvar_diff_sample_postord.rds")

###################### with relative

library(ggplot2)
library(tidyr)

library(scales)
combined_data <- do.call(rbind, lapply(c("sum_var_diff_strict", "sum_var_diff_rel", "sum_var_diff_sample", "sum_var_diff_point_postord", "sum_var_diff_sample_postord", "sum_var_diff_no_ace"), function(dataset_name) {
  dataset <- get(dataset_name)  # Dynamically get the dataset by name
  do.call(rbind, lapply(names(dataset), function(level) {
    do.call(rbind, lapply(seq_along(dataset[[level]]), function(replicate_index) {
      data <- dataset[[level]][[replicate_index]]
      df <- data.frame(
        replicate = replicate_index,
        all = data$all,
        high = data$fossil_high,
        mid = data$fossil_med,
        low = data$fossil_low,
        living = data$living,
        level = level,
        threshold_method = dataset_name  # Add a column for treatment
      )
      df_long <- pivot_longer(df, cols = c("living", "low", "mid", "high", "all"),
                              names_to = "preservation_level", values_to = "distance")
      df_long$preservation_level <- factor(df_long$preservation_level, levels = c("living", "low", "mid", "high", "all"))
      return(df_long)
    }))
  }))
}))

combined_data$threshold_method <- factor(combined_data$threshold_method, levels = c("sum_var_diff_strict", "sum_var_diff_rel", "sum_var_diff_sample", "sum_var_diff_point_postord", "sum_var_diff_sample_postord", "sum_var_diff_no_ace"), labels = c("Strict majority rule", "Relative majority rule", "No majority rule (sample)", "Post-ord ACE (point)", "Post-ord ACE (sample)", "No ACE"))

# log10_reverse_trans <- trans_new(
#   name = "log10_reverse",
#   transform = function(x) -log10(x),  # Apply log10 and reverse
#   inverse = function(x) 10^(-x)      # Inverse of the transformation
# )

colors <- c("slow" = "#8805A8", "med" = "#00B945", "fast" = "#FFB600")
# Generate the combined plot

p <- ggplot(combined_data, aes(x = preservation_level, y = distance, color = level, fill = level)) +
  geom_boxplot(
    alpha = 0.7,
    outlier.size = 0.5,
    outlier.alpha = 0.6,
    position = position_dodge(width = 0.8)
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6, linetype = "dashed") +
  facet_wrap(~ threshold_method, ncol = 6, nrow = 1) +
  labs(
    x = "Fossil sampling",
    y = "Relative disparity error",
    color = "Transition Rate",
    fill = "Transition Rate"
  ) +
  theme_minimal() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) + 
  theme(
    axis.text = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.ticks = element_line(size = 1),
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 1),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  coord_trans(y = "pseudo_log") +
  scale_y_reverse()

ggsave("../Manuscript/figures/sumvar.png", p, width = 12, height = 8, dpi = 700, bg = "white")


### without rel


combined_data <- do.call(rbind, lapply(c("sum_var_diff_strict",  "sum_var_diff_sample", "sum_var_diff_point_postord", "sum_var_diff_sample_postord", "sum_var_diff_no_ace"), function(dataset_name) {
  dataset <- get(dataset_name)  # Dynamically get the dataset by name
  do.call(rbind, lapply(names(dataset), function(level) {
    do.call(rbind, lapply(seq_along(dataset[[level]]), function(replicate_index) {
      data <- dataset[[level]][[replicate_index]]
      df <- data.frame(
        replicate = replicate_index,
        all = data$all,
        high = data$fossil_high,
        mid = data$fossil_med,
        low = data$fossil_low,
        living = data$living,
        level = level,
        threshold_method = dataset_name  # Add a column for treatment
      )
      df_long <- pivot_longer(df, cols = c("living", "low", "mid", "high", "all"),
                              names_to = "preservation_level", values_to = "distance")
      df_long$preservation_level <- factor(df_long$preservation_level, levels = c("living", "low", "mid", "high", "all"))
      return(df_long)
    }))
  }))
}))

combined_data$threshold_method <- factor(combined_data$threshold_method, levels = c("sum_var_diff_strict",  "sum_var_diff_sample", "sum_var_diff_point_postord", "sum_var_diff_sample_postord", "sum_var_diff_no_ace"), labels = c("Strict majority rule", "No majority rule (sample)", "Post-ord ACE (point)", "Post-ord ACE (sample)", "No ACE"))

q <- ggplot(combined_data, aes(x = preservation_level, y = distance, color = level, fill = level)) +
  geom_boxplot(
    alpha = 0.7,
    outlier.size = 0.5,
    outlier.alpha = 0.6,
    position = position_dodge(width = 0.8)
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6, linetype = "dashed") +
  facet_wrap(~ threshold_method, ncol = 5, nrow = 1) +
  labs(
    x = "Fossil sampling",
    y = "Relative disparity error",
    color = "Transition Rate",
    fill = "Transition Rate"
  ) +
  theme_minimal() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) + 
  theme(
    axis.text = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.ticks = element_line(size = 1),
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 1),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  coord_trans(y = "pseudo_log") +
  scale_y_reverse()

ggsave("../Manuscript/figures/sumvar_norel.png", q, width = 12, height = 8, dpi = 700, bg = "white")

##############################################################################################################




neighbours_diff_strict <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/diff_neighbours_strict_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

neighbours_diff_sample <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/diff_neighbours_sample_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

neighbours_diff_sample_median <- lapply(neighbours_diff_sample, lapply, lapply, median)

neighbours_diff_rel <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/diff_neighbours_rel_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

neighbours_diff_no_ace <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/diff_neighbours_no_ace_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

neighbours_diff_point_postord <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/diff_neighbours_point_postord_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)


neighbours_diff_sample_postord <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/diff_neighbours_sample_postord_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)


saveRDS(neighbours_diff_strict, "../Data/cluster/discrete/dispdiffs/neighbours_diff_strict.rds")
saveRDS(neighbours_diff_sample, "../Data/cluster/discrete/dispdiffs/neighbours_diff_sample.rds")
saveRDS(neighbours_diff_rel, "../Data/cluster/discrete/dispdiffs/neighbours_diff_rel.rds")
saveRDS(neighbours_diff_no_ace, "../Data/cluster/discrete/dispdiffs/neighbours_diff_no_ace.rds")
saveRDS(neighbours_diff_point_postord, "../Data/cluster/discrete/dispdiffs/neighbours_diff_point_postord.rds")
saveRDS(neighbours_diff_sample_postord, "../Data/cluster/discrete/dispdiffs/neighbours_diff_sample_postord.rds")


library(ggplot2)
library(tidyr)

library(scales)
combined_data <- do.call(rbind, lapply(c("neighbours_diff_strict", "neighbours_diff_rel", "neighbours_diff_sample", "neighbours_diff_point_postord", "neighbours_diff_sample_postord", "neighbours_diff_no_ace"), function(dataset_name) {
  dataset <- get(dataset_name)  # Dynamically get the dataset by name
  do.call(rbind, lapply(names(dataset), function(level) {
    do.call(rbind, lapply(seq_along(dataset[[level]]), function(replicate_index) {
      data <- dataset[[level]][[replicate_index]]
      df <- data.frame(
        replicate = replicate_index,
        all = data$all,
        high = data$fossil_high,
        mid = data$fossil_med,
        low = data$fossil_low,
        living = data$living,
        level = level,
        threshold_method = dataset_name  # Add a column for treatment
      )
      df_long <- pivot_longer(df, cols = c("living", "low", "mid", "high", "all"),
                              names_to = "preservation_level", values_to = "distance")
      df_long$preservation_level <- factor(df_long$preservation_level, levels = c("living", "low", "mid", "high", "all"))
      return(df_long)
    }))
  }))
}))

combined_data$threshold_method <- factor(combined_data$threshold_method, levels = c("neighbours_diff_strict", "neighbours_diff_rel", "neighbours_diff_sample", "neighbours_diff_point_postord", "neighbours_diff_sample_postord", "neighbours_diff_no_ace"), labels = c("Strict majority rule", "Relative majority rule", "No majority rule (sample)", "Post-ord ACE (point)", "Post-ord ACE (sample)", "No ACE"))

# log10_reverse_trans <- trans_new(
#   name = "log10_reverse",
#   transform = function(x) -log10(x),  # Apply log10 and reverse
#   inverse = function(x) 10^(-x)      # Inverse of the transformation
# )

colors <- c("slow" = "#8805A8", "med" = "#00B945", "fast" = "#FFB600")
# Generate the combined plot

p <- ggplot(combined_data, aes(x = preservation_level, y = distance, color = level, fill = level)) +
  geom_boxplot(
    alpha = 0.7,
    outlier.size = 0.5,
    outlier.alpha = 0.6,
    position = position_dodge(width = 0.8)
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6, linetype = "dashed") +
  facet_wrap(~ threshold_method, ncol = 6, nrow = 1) +
  labs(
    x = "Fossil sampling",
    y = "Relative disparity error",
    color = "Transition Rate",
    fill = "Transition Rate"
  ) +
  theme_minimal() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) + 
  theme(
    axis.text = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.ticks = element_line(size = 1),
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 1),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  coord_trans(y = "pseudo_log") +
  scale_y_reverse()

ggsave("../Manuscript/figures/neighbours.png", p, width = 12, height = 8, dpi = 700, bg = "white")


### without rel


combined_data <- do.call(rbind, lapply(c("neighbours_diff_strict",  "neighbours_diff_sample", "neighbours_diff_point_postord", "neighbours_diff_sample_postord", "neighbours_diff_no_ace"), function(dataset_name) {
  dataset <- get(dataset_name)  # Dynamically get the dataset by name
  do.call(rbind, lapply(names(dataset), function(level) {
    do.call(rbind, lapply(seq_along(dataset[[level]]), function(replicate_index) {
      data <- dataset[[level]][[replicate_index]]
      df <- data.frame(
        replicate = replicate_index,
        all = data$all,
        high = data$fossil_high,
        mid = data$fossil_med,
        low = data$fossil_low,
        living = data$living,
        level = level,
        threshold_method = dataset_name  # Add a column for treatment
      )
      df_long <- pivot_longer(df, cols = c("living", "low", "mid", "high", "all"),
                              names_to = "preservation_level", values_to = "distance")
      df_long$preservation_level <- factor(df_long$preservation_level, levels = c("living", "low", "mid", "high", "all"))
      return(df_long)
    }))
  }))
}))

combined_data$threshold_method <- factor(combined_data$threshold_method, levels = c("neighbours_diff_strict",  "neighbours_diff_sample", "neighbours_diff_point_postord", "neighbours_diff_sample_postord", "neighbours_diff_no_ace"), labels = c("Strict majority rule", "No majority rule (sample)", "Post-ord ACE (point)", "Post-ord ACE (sample)", "No ACE"))

q <- ggplot(combined_data, aes(x = preservation_level, y = distance, color = level, fill = level)) +
  geom_boxplot(
    alpha = 0.7,
    outlier.size = 0.5,
    outlier.alpha = 0.6,
    position = position_dodge(width = 0.8)
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6, linetype = "dashed") +
  facet_wrap(~ threshold_method, ncol = 5, nrow = 1) +
  labs(
    x = "Fossil sampling",
    y = "Relative disparity error",
    color = "Transition Rate",
    fill = "Transition Rate"
  ) +
  theme_minimal() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) + 
  theme(
    axis.text = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.ticks = element_line(size = 1),
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 1),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  coord_trans(y = "pseudo_log") +
  scale_y_reverse()


ggsave("../Manuscript/figures/neighbours_withoutrel.png", q, width = 12, height = 8, dpi = 700, bg = "white")

##############################################################################################################


##' rescaled

rate_names <- c("slow", "med", "fast")


sum_var_diff_strict <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/rescaled/diff_sumvar_strict_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

sum_var_diff_sample <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/rescaled/diff_sumvar_sample_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

sum_var_diff_sample_median <- lapply(sum_var_diff_sample, lapply, lapply, median)

sum_var_diff_rel <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/rescaled/diff_sumvar_rel_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

sum_var_diff_no_ace <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/rescaled/diff_sumvar_no_ace_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

sum_var_diff_point_postord <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/rescaled/diff_sumvar_point_postord_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)


sum_var_diff_sample_postord <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/raw/rescaled/diff_sumvar_sample_postord_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)



combined_data <- do.call(rbind, lapply(c("sum_var_diff_strict", "sum_var_diff_rel", "sum_var_diff_sample", "sum_var_diff_point_postord", "sum_var_diff_sample_postord", "sum_var_diff_no_ace"), function(dataset_name) {
  dataset <- get(dataset_name)  # Dynamically get the dataset by name
  do.call(rbind, lapply(names(dataset), function(level) {
    do.call(rbind, lapply(seq_along(dataset[[level]]), function(replicate_index) {
      data <- dataset[[level]][[replicate_index]]
      df <- data.frame(
        replicate = replicate_index,
        all = data$all,
        high = data$fossil_high,
        mid = data$fossil_med,
        low = data$fossil_low,
        living = data$living,
        level = level,
        threshold_method = dataset_name  # Add a column for treatment
      )
      df_long <- pivot_longer(df, cols = c("living", "low", "mid", "high", "all"),
                              names_to = "preservation_level", values_to = "distance")
      df_long$preservation_level <- factor(df_long$preservation_level, levels = c("living", "low", "mid", "high", "all"))
      return(df_long)
    }))
  }))
}))

combined_data$threshold_method <- factor(combined_data$threshold_method, levels = c("sum_var_diff_strict", "sum_var_diff_rel", "sum_var_diff_sample", "sum_var_diff_point_postord", "sum_var_diff_sample_postord", "sum_var_diff_no_ace"), labels = c("Strict majority rule", "Relative majority rule", "No majority rule (sample)", "Post-ord ACE (point)", "Post-ord ACE (sample)", "No ACE"))

# log10_reverse_trans <- trans_new(
#   name = "log10_reverse",
#   transform = function(x) -log10(x),  # Apply log10 and reverse
#   inverse = function(x) 10^(-x)      # Inverse of the transformation
# )

colors <- c("slow" = "#8805A8", "med" = "#00B945", "fast" = "#FFB600")
# Generate the combined plot

p <- ggplot(combined_data, aes(x = preservation_level, y = distance, color = level, fill = level)) +
  geom_boxplot(
    alpha = 0.7,
    outlier.size = 0.5,
    outlier.alpha = 0.6,
    position = position_dodge(width = 0.8)
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6, linetype = "dashed") +
  facet_wrap(~ threshold_method, ncol = 6, nrow = 1) +
  labs(
    x = "Fossil sampling",
    y = "Relative disparity error",
    color = "Transition Rate",
    fill = "Transition Rate"
  ) +
  theme_minimal() +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) + 
  theme(
    axis.text = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.ticks = element_line(size = 1),
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 1),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  coord_trans(y = "pseudo_log") +
  scale_y_reverse()



#################################################################################################################
## ANOVA


combined_data$abs_distance <- abs(combined_data$distance)


model <- aov(abs_distance ~ level * threshold_method * preservation_level, data = combined_data)
summary(model)

par(mfrow = c(2,2))
plot(model)  # Residual plots
shapiro.test(residuals(model))



## model assumptions violated

combined_data$log_abs_distance <- log(combined_data$abs_distance + 0.001)
model_log <- aov(log_abs_distance ~ level * threshold_method * preservation_level, data = combined_data)
plot(model_log)  # Residual plots



# Use aggregate function (base R equivalent of group_by + summarise)
combined_data_balanced <- aggregate(
  distance ~ replicate + level + threshold_method + preservation_level, 
  data = combined_data, 
  FUN = median, 
  na.rm = TRUE
)


combined_data_balanced$abs_distance <- abs(combined_data_balanced$distance)
combined_data_balanced$log_abs_distance <- log(combined_data_balanced$abs_distance + 0.001)


model_balanced <- aov(log_abs_distance ~ level * threshold_method * preservation_level, 
                     data = combined_data_balanced)

par(mfrow = c(2,2))
plot(model_balanced)
shapiro.test(residuals(model_balanced))


library(agricolae)

# Tukey HSD test for each factor
tukey_rate <- HSD.test(model_balanced, "level", group = TRUE)
tukey_method <- HSD.test(model_balanced, "threshold_method", group = TRUE)
tukey_preservation <- HSD.test(model_balanced, "preservation_level", group = TRUE)

# View rankings with statistical groups
print("Rate rankings:")
tukey_rate$groups

print("Method rankings:")
tukey_method$groups

print("Preservation rankings:")
tukey_preservation$groups


rate_method_means <- emmeans(model_balanced, ~ level * threshold_method)
summary(rate_method_means)

# Rank by performance with statistical significance
pairs(rate_method_means, adjust = "tukey")
cld(rate_method_means, alpha = 0.05, Letters = letters)

# Same for other interactions
rate_preservation_means <- emmeans(model_balanced, ~ level * preservation_level)
method_preservation_means <- emmeans(model_balanced, ~ threshold_method * preservation_level)
summary(method_preservation_means)
cld(method_preservation_means, alpha = 0.05, Letters = letters)

