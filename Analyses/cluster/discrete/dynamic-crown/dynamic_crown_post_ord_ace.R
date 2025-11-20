args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
tree_size <- args[2]
job_id <- as.numeric(args[3])
library(parallel)
library(treats)
library(dispRity)
library(ape)

base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/discrete/crown/", tree_size, "/")
# job_id <- 8558401

write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}

# Checkpoint helpers
checkpoint_dir <- paste0(base_path, "checkpoints/")
if(!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

checkpoint_file <- paste0(checkpoint_dir, job_id, "_post_ord_progress_", sprintf("%03d.rds", replicate_id))

fossil_matrices <- readRDS(write.path("matrices", "fossil_matrices_%03d.rds"))
fossil_trees <- lapply(fossil_matrices, lapply, function(level){
  tree <- level$tree
})

tip_mat <- lapply(fossil_matrices, lapply, function(x){
    mat <- x$matrix
    tip_mat <- mat[grepl("^t", rownames(mat)), ]
    return(tip_mat)
})


ord_fossil_tips <- lapply(tip_mat, lapply, function(x){
  dist <- char.diff(x, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})


# Extract living labels
# Split each matrix into 4 partitions
split_matrices <- lapply(ord_fossil_tips, function(rate_list) {
  lapply(rate_list, function(mat) {
    n_dims <- ncol(mat)
    partition_size <- ceiling(n_dims / 6)
    
    # Create 4 partitions
    partitions <- list()
    for(i in 1:6) {
      start_col <- (i - 1) * partition_size + 1
      end_col <- min(i * partition_size, n_dims)
      partitions[[i]] <- mat[, start_col:end_col, drop = FALSE]
    }
    
    return(partitions)
  })
})

# Create tasks: 60 tasks (3 rates × 5 levels × 6 partitions)
tasks_post_ord <- expand.grid(
  rate = names(ord_fossil_tips), 
  fossil_level = names(ord_fossil_tips[[1]]), 
  partition = 1:6, 
  stringsAsFactors = FALSE
)

tasks_post_ord <- tasks_post_ord[order(tasks_post_ord$fossil_level, 
                                        tasks_post_ord$rate, 
                                        tasks_post_ord$partition), ]

cat("Total tasks:", nrow(tasks_post_ord), "\n")

# Check for existing checkpoint
if(file.exists(checkpoint_file)) {
  cat("Loading checkpoint...\n")
  checkpoint <- readRDS(checkpoint_file)
  res_post_ord <- checkpoint$results
  completed_tasks <- checkpoint$completed
  cat("Resuming from", length(completed_tasks), "completed tasks\n")
} else {
  res_post_ord <- vector("list", nrow(tasks_post_ord))
  completed_tasks <- c()
}

remaining_tasks <- setdiff(seq_len(nrow(tasks_post_ord)), completed_tasks)
cat("Running", length(remaining_tasks), "remaining tasks on 18 cores with batches of 18...\n")
if(length(remaining_tasks) > 0) {
    batch_size <- 18  # Process 18 tasks at a time (5 batches)
    n_batches <- ceiling(length(remaining_tasks) / batch_size)
    
    for(batch_idx in 1:n_batches) {
        start_idx <- (batch_idx - 1) * batch_size + 1
        end_idx <- min(batch_idx * batch_size, length(remaining_tasks))
        batch_tasks <- remaining_tasks[start_idx:end_idx]
        
        cat("Batch", batch_idx, "/", n_batches, ": processing", length(batch_tasks), "tasks\n")
        
        # Run batch
        res_batch <- mclapply(batch_tasks, function(i) {
            task <- tasks_post_ord[i, ]
            matrix_partition <- split_matrices[[task$rate]][[task$fossil_level]][[task$partition]]
            fossil_tree <- fossil_trees[[task$rate]][[task$fossil_level]]
            
            result <- tryCatch({
                multi.ace(matrix_partition, fossil_tree, models = "BM", output = "multi.ace")
            }, error = function(e) {
                cat("ERROR:", task$rate, task$fossil_level, "partition", task$partition, e$message, "\n")
                NULL
            })
            
            return(result)
        }, mc.cores = 18)
        
        # Store batch results
        for(j in seq_along(batch_tasks)) {
            res_post_ord[[batch_tasks[j]]] <- res_batch[[j]]
        }
        
        # Update completed list
        completed_tasks <- c(completed_tasks, batch_tasks)
        
        # Save checkpoint after EACH batch
        saveRDS(list(results = res_post_ord, completed = completed_tasks), checkpoint_file)
        cat("Checkpoint saved:", length(completed_tasks), "/", nrow(tasks_post_ord), "completed\n")
    }
}

# Reconstruct: combine 4 partitions back into full multi.ace objects
post_ord_ace <- list()

for(rate in names(ord_fossil_tips)) {
  post_ord_ace[[rate]] <- list()
  
  for(fossil_level in names(ord_fossil_tips[[1]])) {
    # Find indices for all 4 partitions of this rate×fossil_level
    partition_indices <- which(
      tasks_post_ord$rate == rate & 
      tasks_post_ord$fossil_level == fossil_level
    )
    
    # Get the 4 partition results (should be in order: partition 1,2,3,4)
    partitions <- lapply(partition_indices, function(idx) res_post_ord[[idx]])
    
    # Check if all 4 partitions completed
    if(any(sapply(partitions, is.null))) {
      cat("WARNING: Missing partitions for", rate, fossil_level, "\n")
      post_ord_ace[[rate]][[fossil_level]] <- NULL
      next
    }
    
    # Combine the 4 partitions into one multi.ace object
    # Start with partition 1
    combined <- partitions[[1]]
    
    # Add dimensions from partitions 2, 3, 4
    for(p in 2:6) {
      # Combine ace estimates (column-bind trait dimensions)
      combined$ace <- cbind(combined$ace, partitions[[p]]$ace)
      combined$CI95 <- cbind(combined$CI95, partitions[[p]]$CI95)
    }
    
    post_ord_ace[[rate]][[fossil_level]] <- combined
  }
}

saveRDS(post_ord_ace, write.path("anc", "post_ord_ace_%03d.rds"))

cat("Starting point post-ord ACE...\n")

# Point estimates
point_post_ord_ace <- lapply(post_ord_ace, lapply, function(x) {
  if(is.null(x)) return(NULL)
  multi.ace(x, output = "combined.matrix")
})

# point_post_ord_ace_living <- Map(function(rate_anc, rate_labels) {
#     Map(function(fossil_anc, label_anc) {
#       if(is.null(fossil_anc)) return(NULL)
#       fossil_anc[label_anc, , drop = FALSE]
#   }, rate_anc, rate_labels)
# }, point_post_ord_ace, labels)

# names(point_post_ord_ace_living) <- names(point_post_ord_ace)

cat("Starting sample post-ord ACE (100 samples)...\n")

# Sample estimates
trait_normal <- list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))

sample_post_ord_ace <- lapply(post_ord_ace, lapply, function(x) {
  if(is.null(x)) return(NULL)
  multi.ace(x, sample = 100, sample.fun = trait_normal, output = "combined.matrix")
})

# sample_post_ord_ace_living <- Map(function(rate_anc, rate_labels) {
#     Map(function(fossil_anc, label_anc) {
#       if(is.null(fossil_anc)) return(NULL)
#       lapply(fossil_anc, function(rep) {
#         rep[label_anc, , drop = FALSE]
#       })
#   }, rate_anc, rate_labels)
# }, sample_post_ord_ace, labels)

# names(sample_post_ord_ace_living) <- names(sample_post_ord_ace)

# Save final outputs
saveRDS(point_post_ord_ace, write.path("ord", "post_ord_point_%03d.rds"))
saveRDS(sample_post_ord_ace, write.path("ord", "post_ord_sample_%03d.rds"))

# Clean up checkpoint
if(file.exists(checkpoint_file)) file.remove(checkpoint_file)

cat("Finished replicate", replicate_id, "\n")