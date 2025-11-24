args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])  # 1-1500
tree_size <- args[2]
job_id <- as.numeric(args[3])
offset <- as.numeric(args[4])
library(treats)
library(dispRity)
library(ape)

task_id <- offset + task_id ## for sending in two batches (batch 1: 1:1000, batch 2: 1001-1500)

base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/discrete/crown/", tree_size, "/")

cat("Starting task", task_id, "\n")

# Create all 1500 ACE runs: 100 reps × 3 rates × 5 levels
rates <- c("fast", "med", "slow")
levels <- c("all", "fossil_high", "fossil_low", "fossil_med", "living")

all_ace_runs <- expand.grid(
  fossil_level = levels,  # Changes FASTEST
  rate = rates,
  replicate = 1:100,      # Changes SLOWEST
  stringsAsFactors = FALSE
)

# Get this specific task
run <- all_ace_runs[task_id, ]
set.seed(100 + run$replicate)

cat("Task", task_id, "- Replicate:", run$replicate, 
    "Rate:", run$rate, "Level:", run$fossil_level, "\n")

write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", 
                sprintf(filename, run$replicate)))
}

# Load data for this replicate
cat("Loading data...\n")
fossil_matrices <- readRDS(write.path("matrices", "fossil_matrices_%03d.rds"))

# Extract specific rate/level
level_data <- fossil_matrices[[run$rate]][[run$fossil_level]]
tree_data <- level_data$tree
matrix_data <- level_data$matrix

# Extract tips only
tip_matrix <- matrix_data[grepl("^t", rownames(matrix_data)), ]

cat("Performing ordination...\n")
# Ordinate
dist <- char.diff(tip_matrix, method = "mord", by.col = FALSE)
ord_matrix <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points

cat("Running ACE on", ncol(ord_matrix), "dimensions...\n")
start_time <- Sys.time()

# Run ACE
ace_result <- tryCatch({
  multi.ace(ord_matrix, tree_data, models = "BM", output = "multi.ace")
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  return(NULL)
})

elapsed <- difftime(Sys.time(), start_time, units = "hours")
cat("ACE completed in", round(elapsed, 2), "hours\n")

# Create save directory
save_dir <- paste0(base_path, "checkpoints_ace/rep", sprintf("%03d", run$replicate), 
                   "/", run$rate, "/", run$fossil_level)
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# Save ACE result with job_id
ace_file <- paste0(save_dir, "/", job_id, "_post_ord_ace.rds")
saveRDS(ace_result, ace_file)
cat("Saved ACE result:", ace_file, "\n")

cat("Generating point estimates...\n")

# Point estimates
point_post_ord_ace <- tryCatch({
  multi.ace(ace_result, output = "combined.matrix")
}, error = function(e) {
  cat("ERROR in point estimates:", e$message, "\n")
  return(NULL)
})

cat("Generating 100 samples...\n")

# Sample estimates
trait_normal <- list(
  fun = rnorm, 
  param = list(mean = mean, sd = function(x) diff(range(x)) / 4)
)

sample_post_ord_ace <- tryCatch({
  multi.ace(ace_result, sample = 100, sample.fun = trait_normal, output = "combined.matrix")
}, error = function(e) {
  cat("ERROR in sample estimates:", e$message, "\n")
  return(NULL)
})

# Save point and sample estimates with job_id
point_file <- paste0(save_dir, "/", job_id, "_point_post_ord_ace.rds")
sample_file <- paste0(save_dir, "/", job_id, "_sample_post_ord_ace.rds")

saveRDS(point_post_ord_ace, point_file)
saveRDS(sample_post_ord_ace, sample_file)

cat("Saved point estimates:", point_file, "\n")
cat("Saved sample estimates:", sample_file, "\n")

cat("Task", task_id, "complete!\n")