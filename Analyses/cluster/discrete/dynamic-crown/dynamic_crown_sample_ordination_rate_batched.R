args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])  # 1-1000
tree_size <- args[2]
job_id <- as.numeric(args[3])
library(parallel)
library(treats)

base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/discrete/crown/", tree_size, "/")

cat("Starting task", task_id, "\n")

# ✅ Simple: List all 150,000 ordinations we need to do
rates <- c("fast", "med", "slow")
levels <- c("all", "fossil_high", "fossil_low", "fossil_med", "living")

all_ordinations <- expand.grid(
  sample = 1:100,         # ✅ Changes FASTEST (innermost loop)
  fossil_level = levels,  
  rate = rates,           
  replicate = 1:100,      # ✅ Changes SLOWEST (outermost loop)
  stringsAsFactors = FALSE
)

# ✅ Simple: Each of 1000 tasks does 150 ordinations
ordinations_per_task <- 150
start <- (task_id - 1) * ordinations_per_task + 1
end <- min(task_id * ordinations_per_task, nrow(all_ordinations))

cat("Task", task_id, "will process ordinations", start, "to", end, "\n")

# ✅ Simple: Loop through my assigned ordinations
current_replicate <- -1
sample_data <- NULL

for(ord_num in start:end) {
  
  # Get what this ordination needs
  ord <- all_ordinations[ord_num, ]
  
  # Load replicate data (only if we haven't already)
  if(ord$replicate != current_replicate) {
    cat("Loading replicate", ord$replicate, "...\n")
    file <- paste0(base_path, "anc/", job_id, "_pre_ord_sample_", 
                   sprintf("%03d", ord$replicate), ".rds")
    sample_data <- readRDS(file)
    current_replicate <- ord$replicate
  }
  
  # Do the ordination
  matrix <- sample_data[[ord$rate]][[ord$fossil_level]][[ord$sample]]
  dist <- char.diff(matrix, method = "mord", by.col = FALSE)
  ord_result <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
  
  # Save it
  save_dir <- paste0(base_path, "checkpoints/rep", sprintf("%03d", ord$replicate), 
                     "/", ord$rate, "/", ord$fossil_level)
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  save_file <- paste0(save_dir, "/", ord$sample, ".rds")
  saveRDS(ord_result, save_file)
  
  # Progress update every 10 ordinations
  if((ord_num - start + 1) %% 10 == 0) {
    cat("  Done", ord_num - start + 1, "of", ordinations_per_task, "\n")
  }
}

cat("Task", task_id, "complete!\n")