args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])  # 1-1500
tree_size <- args[2]
job_id <- as.numeric(args[3])
library(treats)
library(dispRity)
library(ape)

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

cat("Running ACE (this may take a while for 'all' level)...\n")
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

# Save checkpoint
checkpoint_dir <- paste0(base_path, "checkpoints_ace/rep", 
                        sprintf("%03d", run$replicate), "/", run$rate)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

checkpoint_file <- paste0(checkpoint_dir, "/", run$fossil_level, ".rds")
saveRDS(ace_result, checkpoint_file)

cat("Saved:", checkpoint_file, "\n")
cat("Task", task_id, "complete!\n")