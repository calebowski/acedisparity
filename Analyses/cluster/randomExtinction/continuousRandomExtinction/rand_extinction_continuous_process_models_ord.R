args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
model_name <- args[2] 
tree_size <- args[3]

library(treats)
library(parallel)
library(MASS)


cat("Starting replicate", replicate_id, "\n")
set.seed(100 + replicate_id)

base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/continuous/", tree_size, "t/")
job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")

write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}

source("/users/bip24cns/acedisparity/discrete/scripts/utility.R") ## read in core utility functions
source("/users/bip24cns/acedisparity/randomExtinction/scripts/find.extinction.time.R")


tree <- read.tree(write.path("trees", "tree_%03d.tre")) ## read in tree generated from <generate_tree.R> script
tree <- fix.zero.branches(tree)
tree <- set.root.time(tree)

BM.trend.process <- function(x0 = 0, edge.length = 1, Sigma = diag(length(x0)), trend = 0.1, ...) {
      # Square root gives more trend than log but less than linear
      drift <- trend * log(edge.length + 1)
      if(edge.length < (0.01 * Sigma[1,1])){drift <- 0} ## if edge length is smaller than 1% of variance, drift is 0.
      return(t(MASS::mvrnorm(n = 1, mu = x0 + drift, Sigma = Sigma * edge.length, ...)))
}

traits <- switch(model_name,
  "bm" = make.traits(process = BM.process, n = 100, 
                     process.args = list(Sigma = diag(0.25, 100))),
  
  "bm_t" = make.traits(process = BM.trend.process, n = 100, 
                      process.args = list(Sigma = diag(0.25, 100), trend = 0.3)),
  
  "ou_w" = make.traits(process = OU.process, n = 100, 
                      process.args = list(alpha = (log(2) / (tree_height / 2)), 
                                        Sigma = diag(0.25, 100))),
  
  "ou_st" = make.traits(process = OU.process, n = 100, 
                       process.args = list(alpha = (log(2) / (tree_height / 10)), 
                                         Sigma = diag(0.25, 100))),
  
  "ou_sh" = make.traits(process = OU.process, n = 100, 
                       process.args = list(optimum = 2, 
                                         alpha = (log(2) / (tree_height / 2)), 
                                         Sigma = diag(0.25, 100))),
  
  stop("Unknown model: ", model_name)
)


mat <- map.traits(traits, tree)$data ## simulate traits onto tree

cat("Traits simulated for", model_name, "\n")

saveRDS(mat, write.path("matrices", paste0(model_name, "_matrix_%03d.rds")))


cat("Trait matrices saved for replicate", replicate_id, model_name, "\n")

source("/users/bip24cns/acedisparity/discrete/scripts/fossil.pres.R")

living <- remove.fossil(mat, trees = tree, type = "continuous")
fossilised_high <- fossil.pres(mat,  trees = tree, preservation = 0.5, type = "continuous")
all_fossil <- fossil.pres(mat,  trees = tree, preservation = 1.0, type = "continuous")
fossilised_med <- fossil.pres(mat,  trees = tree, preservation = 0.15, type = "continuous")
fossilised_low <- fossil.pres(mat,  trees = tree, preservation = 0.05, type = "continuous")


fossil_matrices <- list(
      all = all_fossil,
      fossil_high = fossilised_high,
      fossil_med = fossilised_med,
      fossil_low = fossilised_low,
      living = living
)

# names(fossil_matrices) <- names(mat)

fossil_trees <- lapply(fossil_matrices, function(level){
  tree <- level$tree
})


saveRDS(fossil_matrices, write.path("matrices", paste0(model_name, "_fossil_matrices_%03d.rds")))
saveRDS(fossil_trees, write.path("trees", paste0(model_name, "_fossil_trees_%03d.rds")))

tasks  <- expand.grid(fossil_level = names(fossil_matrices), stringsAsFactors = FALSE) 

cl <- makeCluster(5) ## make 25 core cluster (one for each task)
clusterEvalQ(cl, library(treats))
clusterExport(cl, c("fossil_matrices", "tasks"))

res <- parLapply(cl, seq_len(nrow(tasks)), function(i){ ## loop over each model combination by row
    task <- tasks[i,]
    level <- fossil_matrices[[task]]
    tryCatch({
    multi.ace(level$matrix, level$tree, models = "BM", output = "multi.ace")
  }, error = function(e) {
    cat("ERROR:", task$fossil_level, e$message, "\n") ## error handeling
    NULL
  })
})

stopCluster(cl)

fossil_anc <- list()
for(i in seq_along(res)) {
  l <- tasks$fossil_level[i]
  fossil_anc[[l]] <- res[[i]]
}