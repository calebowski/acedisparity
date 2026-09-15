library(ggplot2)

data_dir <- "/home/caleb/Documents/PhD/acedisparity/Cluster/Data/revisions/wagner2/discrete/dist_disparity"

methods <- c(
  ace = "ACE reconstructed",
  sampled_only = "Sampled taxa only"
)

file_patterns <- c(
  ace = "^11532243_ace_dist_disparity_errors_[0-9]{3}\\.rds$",
  sampled_only = "^11532243_sampled_only_dist_disparity_errors_[0-9]{3}\\.rds$"
)

load_method_data <- function(method_name, file_pattern) {
  files <- list.files(
    data_dir,
    pattern = file_pattern,
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No files found for method: ", method_name)
  }

  do.call(rbind, lapply(seq_along(files), function(replicate_id) {
    error_data <- readRDS(files[replicate_id])

    do.call(rbind, lapply(names(error_data), function(rate_name) {
      rate_data <- error_data[[rate_name]]

      do.call(rbind, lapply(names(rate_data), function(fossil_level) {
        data.frame(
          replicate = replicate_id,
          rate = rate_name,
          preservation_level = fossil_level,
          error = as.numeric(rate_data[[fossil_level]]),
          method = method_name,
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))
}

results_df <- do.call(
  rbind,
  Map(load_method_data, names(methods), file_patterns)
)

results_df$method <- factor(
  results_df$method,
  levels = names(methods),
  labels = unname(methods)
)

results_df$rate <- factor(
  results_df$rate,
  levels = c("slow", "med", "fast"),
  labels = c("Slow", "Medium", "Fast")
)

results_df$preservation_level <- factor(
  results_df$preservation_level,
  levels = c( "fossil_low", "fossil_med", "fossil_high", "all"),
  labels = c("5", "15", "50", "100")
)

colours <- c(
  "Slow" = "#8805A8",
  "Medium" = "#00B945",
  "Fast" = "#FFB600"
)

dist_disparity_plot <- ggplot(
  results_df,
  aes(
    x = preservation_level,
    y = error,
    fill = rate,
    colour = rate
  )
) +
  geom_boxplot(
    alpha = 0.7,
    outlier.shape = NA,
    position = position_dodge(width = 0.8)
  ) +
  geom_hline(
    yintercept = 0,
    colour = "black",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  facet_wrap(~ method, ncol = 1) +
  labs(
    x = "Fossil Sampling (%)",
    y = "Absolute Mean Pairwise Distance Error",
    fill = "Transition Rate"
  ) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours, guide = "none") +
#   scale_y_continuous(
#     breaks = c(-1, -0.5, 0, 0.5, 1),
#     labels = c("-1.0", "-0.5", "0", "0.5", "1.0")
#   ) +
#   coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16, colour = "black"),
    axis.title = element_text(size = 20, face = "bold", colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15),
    strip.text = element_text(size = 20, face = "bold"),
    strip.background = element_rect(
      fill = "gray95",
      colour = "black",
      linewidth = 0.3
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.4
    ),
    panel.grid.major.y = element_line(
      colour = "gray90",
      linewidth = 0.3
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# ggsave(
#   "/home/caleb/Documents/PhD/acedisparity/Manuscript/draft/figures/wagner2_dist_disparity_boxplot.png",
#   dist_disparity_plot,
#   width = 14,
#   height = 10,
#   dpi = 700,
#   units = "in",
#   bg = "white"
# )