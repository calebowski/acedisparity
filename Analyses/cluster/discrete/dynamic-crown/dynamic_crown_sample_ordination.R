args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
tree_size <- args[2]
job_id <- as.numeric(args[3])
library(parallel)
library(treats)

base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/discrete/crown/", tree_size, "/")

write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}

# Checkpoint helpers
checkpoint_dir <- paste0(base_path, "checkpoints/")
if(!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

checkpoint_file <- paste0(checkpoint_dir, job_id, "_ord_progress_", sprintf("%03d.rds", replicate_id))

fossil_matrices <- readRDS(write.path("matrices", "fossil_matrices_%03d.rds"))
sample_fossil_anc <- readRDS(write.path("anc", "pre_ord_sample_%03d.rds"))
relative_fossil_anc <- readRDS(write.path("anc", "pre_ord_rel_%03d.rds"))
point_fossil_anc <- readRDS(write.path("anc", "pre_ord_point_%03d.rds"))
matrices <- readRDS(write.path("matrices", "matrices_%03d.rds"))

no_ace_matrix <- lapply(point_fossil_anc, lapply, function(matrix){matrix[!grepl("^n", rownames(matrix)),]}) # remove nodes


ord_rel <- lapply(relative_fossil_anc, lapply, function(rep){
  dist <- char.diff(rep, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
}) ## ordinate full matrix, prune out certain datapoints later.

ord_point <- lapply(point_fossil_anc, lapply,  function(rep){
  dist <- char.diff(rep, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})

ord_true <- lapply(matrices, function(rep){
  dist <- char.diff(rep, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})

ord_no_ace <- lapply(no_ace_matrix, lapply, function(rep){
  dist <- char.diff(rep, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})

saveRDS(ord_no_ace, write.path("ord", "ord_no_ace_%03d.rds"))
saveRDS(ord_true, write.path("ord", "ord_true_%03d.rds"))
saveRDS(ord_rel, write.path("ord", "ord_rel_%03d.rds"))
saveRDS(ord_point, write.path("ord", "ord_point_%03d.rds"))



# extract.living.nodes <- function(fossils) { 
# #   basal_node <- fossils$living$tree$node.label[1]
#   living_nodes <- fossils$living$tree$node.label
#   labels <- lapply(fossils, function(level){
#     tree <- level$tree
#     add_nodes <- c(tree$tip.label, living_nodes)
#     return(add_nodes)
#   })
#   return(labels)
# }



tasks_sample_ord <- expand.grid(
  rate = names(sample_fossil_anc), 
  fossil_level = names(sample_fossil_anc[[1]]), 
  rep_n = 1:100,  # inner samples
  stringsAsFactors = FALSE
)

tasks_sample_ord <- tasks_sample_ord[order(tasks_sample_ord$fossil_level, 
                                        tasks_sample_ord$rate, 
                                        tasks_sample_ord$rep_n), ]

cat("Total tasks:", nrow(tasks_sample_ord), "\n")

# Check for existing checkpoint
if(file.exists(checkpoint_file)) {
  cat("Loading checkpoint...\n")
  checkpoint <- readRDS(checkpoint_file)
  res_ord_sample <- checkpoint$results
  completed_tasks <- checkpoint$completed
  cat("Resuming from", length(completed_tasks), "completed tasks\n")
} else {
  res_ord_sample <- vector("list", nrow(tasks_sample_ord))
  completed_tasks <- c()
}

# Find remaining tasks
remaining_tasks <- setdiff(seq_len(nrow(tasks_sample_ord)), completed_tasks)
cat("Running", length(remaining_tasks), "remaining tasks on 64 cores\n")

if(length(remaining_tasks) > 0) {
  # Process in batches with periodic checkpointing
  batch_size <- 300  # Save every 300 tasks (corresponds to each fossil level)
  n_batches <- ceiling(length(remaining_tasks) / batch_size)
  
  for(batch_idx in 1:n_batches) {
    start_idx <- (batch_idx - 1) * batch_size + 1
    end_idx <- min(batch_idx * batch_size, length(remaining_tasks))
    batch_tasks <- remaining_tasks[start_idx:end_idx]
    
    cat("Batch", batch_idx, "/", n_batches, ": processing", length(batch_tasks), "tasks\n")
    
    # Run batch
    res_batch <- mclapply(batch_tasks, function(i){
        task <- tasks_sample_ord[i, ]
        
        # Get single replicate matrix 
        rep_matrix <- sample_fossil_anc[[task$rate]][[task$fossil_level]][[task$rep_n]]
        
        # Single ordination
        dist <- char.diff(rep_matrix, method = "mord", by.col = FALSE)
        ord_rep <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
        
        return(ord_rep)
    }, mc.cores = 64)
    
    # Store batch results
    for(j in seq_along(batch_tasks)) {
      res_ord_sample[[batch_tasks[j]]] <- res_batch[[j]]
    }
    
    # Update completed list
    completed_tasks <- c(completed_tasks, batch_tasks)
    
    # Save checkpoint
    saveRDS(list(results = res_ord_sample, completed = completed_tasks), checkpoint_file)
    cat("Checkpoint saved:", length(completed_tasks), "/", nrow(tasks_sample_ord), "completed\n")
  }
}

# Reconstruct nested structure: rate -> fossil_level -> rep_n
ord_sample <- list()
for(i in seq_along(res_ord_sample)) {
  rate <- tasks_sample_ord$rate[i]
  fossil_level <- tasks_sample_ord$fossil_level[i]
  rep_n <- tasks_sample_ord$rep_n[i]
  
  # Initialize nested lists if needed
  if(is.null(ord_sample[[rate]])) ord_sample[[rate]] <- list()
  if(is.null(ord_sample[[rate]][[fossil_level]])) ord_sample[[rate]][[fossil_level]] <- list()
  
  # Store result
  ord_sample[[rate]][[fossil_level]][[rep_n]] <- res_ord_sample[[i]]
}

saveRDS(ord_sample, write.path("ord", "ord_sample_%03d.rds"))

# Clean up checkpoint after successful completion
if(file.exists(checkpoint_file)) file.remove(checkpoint_file)

cat("Finished replicate", replicate_id, "\n")