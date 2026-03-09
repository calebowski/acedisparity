args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
tree_size <- args[2]
job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")

library(treats)
library(parallel)
library(MASS)

cat("Starting replicate", replicate_id, "\n")


base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/continuous/", tree_size, "/")
source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")

write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}

# Test 1: Read tree
# cat("\n=== STEP 1: Reading tree ===\n")
# tree_file <- write.path("trees", "tree_%03d.tre")
# cat("Tree file:", tree_file, "\n")

# tryCatch({
#   tree <- read.tree(tree_file)
#   cat("SUCCESS: Read tree with", length(tree$tip.label), "tips\n")
# }, error = function(e) {
#   cat("ERROR reading tree:", e$message, "\n")
#   quit(status = 1)
# })

cat("Ensuring directories exist...\n")
dirs <- c("matrices", "trees", "anc", "ord")
for(dir in dirs) {
  full_path <- paste0(base_path, dir)
  if(!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE)
    cat("Created:", full_path, "\n")
  }
}


cat("\n=== STEP 1: Generating tree ===\n")
set.seed(100 + replicate_id)  # Set seed for reproducible tree generation

tree <- read.tree(paste0("/mnt/parscratch/users/bip24cns/acedisparity/trees/overallDisparity/tree_", tree_size, sprintf("_%03d.tre", replicate_id)))

# tree <- drop.singles(tree)
# tree <- fix.zero.branches(tree)
# tree <- set.root.time(tree)

ages <- tree.age(tree) # get tip ages
extant <- ages$element[ages$ages == 0]
living_tree <- keep.tip(tree, extant) ## root at mrca - alike to crown group analyses.

crown_tree  <- extract.clade(tree, living_tree$node.label[1]) ## this is the crown tree



tree_height <- max(node.depth.edgelength(crown_tree))



cat("Tree generated for replicate", replicate_id, "\n")
bm  <-  make.traits(process = BM.process, n = 20, 
                     process.args = list(Sigma = diag(0.25, 20)))
  
bm_t <-  make.traits(process = BM.trend.process, n = 20, 
                      process.args = list(Sigma = diag(0.25, 20), trend = 0.3))
  
ou_w <-  make.traits(process = OU.process, n = 20, 
                      process.args = list(alpha = (log(2) / (tree_height * 0.75)),
                                        Sigma = diag(0.25, 20)))
  
ou_st <-  make.traits(process = OU.process, n = 20, 
                       process.args = list(alpha = (log(2) / (tree_height * 0.25)),
                                         Sigma = diag(0.25, 20)))
  
ou_sh <-  make.traits(process = OU.process, n = 20, 
                       process.args = list(optimum = 4, 
                                         alpha = (log(2) / (tree_height * 0.75)),
                                         Sigma = diag(0.25, 20)))

# traits <- list(bm = bm, bm_t = bm_t)
traits <- list(bm = bm, bm_t = bm_t, ou_w = ou_w, ou_st = ou_st, ou_sh = ou_sh)


# mat <- map.traits(traits, tree)$data

matrices <- lapply(traits, function(x){map.traits(x, crown_tree)$data})

cat("Traits simulated...\n")

# matrices <- list(bm = bm_matrices, bm_t = bm_matrices_trend, ou_w = ou_w_matrices, ou_st = ou_s_matrices, ou_sh = ou_shift_matrices)

saveRDS(matrices, write.path("matrices", "matrices_%03d.rds"))

cat("Trait matrices saved...\n")

source("/users/bip24cns/acedisparity/discrete/scripts/fossil.pres.R")
# set_seed <- 100 + replicate_id
living <- lapply(matrices, remove.fossil, trees = crown_tree, type = "continuous")
fossilised_high <- lapply(matrices, fossil.pres, trees = crown_tree, preservation = 0.5, type = "continuous", seed = set_seed)
all_fossil <- lapply(matrices, fossil.pres, trees = crown_tree, preservation = 1.0, type = "continuous", seed = set_seed)
fossilised_med <- lapply(matrices, fossil.pres, trees = crown_tree, preservation = 0.15, type = "continuous", seed = set_seed)
fossilised_low <- lapply(matrices, fossil.pres, trees = crown_tree, preservation = 0.05, type = "continuous", seed = set_seed)

fossil_matrices <- lapply(names(matrices), function(level) {
    list(
      all = all_fossil[[level]],
      fossil_high = fossilised_high[[level]],
      fossil_med = fossilised_med[[level]],
      fossil_low = fossilised_low[[level]],
      living = living[[level]]
    )
})


# Assign names to the outer list
names(fossil_matrices) <- names(matrices)


fossil_trees <- lapply(fossil_matrices, lapply, function(level){
  tree <- level$tree
})
cat("Fossil matrices completed...\n")

saveRDS(fossil_matrices, write.path("matrices", "fossil_matrices_%03d.rds"))

tasks  <- expand.grid(model = names(fossil_matrices), fossil_level = names(fossil_matrices[[1]]),  stringsAsFactors = FALSE) 

start_time <- Sys.time()
res_pre_ord_ace <- lapply(seq_len(nrow(tasks)), function(i){ ## loop over each model combination by row
    task <- tasks[i,]
    level <- fossil_matrices[[task$model]][[task$fossil_level]]
    tryCatch({
    multi.ace(level$matrix, level$tree, models = "BM", output = "multi.ace")
  }, error = function(e) {
    cat("ERROR:", task$fossil_level, e$message, "\n") ## error handeling
    NULL
  })
})

elapsed <- difftime(Sys.time(), start_time, units = "mins")

fossil_anc <- list()
for(i in seq_along(res_pre_ord_ace)) {
  m <- tasks$model[i] 
  l <- tasks$fossil_level[i]
  fossil_anc[[m]][[l]] <- res_pre_ord_ace[[i]]
}

cat("Starting point ancestral state estimation...\n")
point_anc <- lapply(fossil_anc, lapply, multi.ace, output = "combined.matrix")  ## takes median of confidence interval

trait_normal  <-  list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4))) ## samples with normal distribution

cat("Starting distribution ancestral state estimation...\n")
sample_anc <- lapply(fossil_anc, lapply, multi.ace, output = "combined.matrix", ml.collapse = list(type = "sample", sample = 100, sample.fun = trait_normal)) ## 100 matrices replicated, 

cat("Ancestral states estimated...\n")

saveRDS(fossil_anc, write.path("anc", "pre_ord_ace_%03d.rds"))
saveRDS(point_anc, write.path("anc", "pre_ord_point_%03d.rds"))
saveRDS(sample_anc, write.path("anc", "pre_ord_sample_%03d.rds"))
cat("Ancestral states saved...\n")


cat("=== STARTING ord_sample ===\n")
  ord_sample <- lapply(sample_anc, lapply, lapply, function(mat){
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })


cat("=== STARTING ord_point ===\n")
  ord_point <- lapply(point_anc, lapply, function(mat){
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })


cat("=== STARTING ord_no_ace ===\n")
  ord_no_ace <- lapply(fossil_matrices, lapply, function(mat){
    prcomp(mat$matrix, scale = FALSE, center = TRUE)$x
  })


cat("=== STARTING ord_true ===\n")
ord_true <- lapply(matrices, function(mat) {
  prcomp(mat, scale = FALSE, center = TRUE)$x
})


cat("Ordinations completed...\n")


cat("=== STARTING ordination of fossil tips for post ordiation ace ===\n")
# ord_fossil_tips <- lapply(fossil_matrices, function(x) {
#     mat <- x$matrix
#     prcomp(mat, scale = FALSE, center = TRUE)$x
#   })

cat("=== STARTING ancestral statte estimation of post ordination ===\n")

tasks_post_ord <- expand.grid(model = names(ord_no_ace), fossil_level = names(ord_no_ace[[1]]), stringsAsFactors = FALSE)

start_time <- Sys.time()
res_post_ord_ace <- lapply(seq_len(nrow(tasks_post_ord)), function(i){ ## loop over each model combination by row
    task <- tasks_post_ord[i,]
    level_matrix <- ord_no_ace[[task$model]][[task$fossil_level]]
    level_tree <- fossil_trees[[task$model]][[task$fossil_level]]
    tryCatch({
    multi.ace(level_matrix, level_tree, models = "BM", output = "multi.ace")
  }, error = function(e) {
    cat("ERROR:", task$fossil_level, e$message, "\n") ## error handeling
    NULL
  })
})

elapsed <- difftime(Sys.time(), start_time, units = "mins")


post_ord_ace  <- list()
for(i in seq_along(res_post_ord_ace)) {
    m <- tasks_post_ord$model[i]
    l <- tasks_post_ord$fossil_level[i]
    post_ord_ace[[m]][[l]] <- res_post_ord_ace[[i]]
}

cat("=== STARTING point estimation post ord ace ===\n")

point_post_ord_ace <- lapply(post_ord_ace, lapply, multi.ace, output = "combined.matrix")

trait_normal <- list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))

sample_post_ord_ace <- lapply(post_ord_ace, lapply, multi.ace, ml.collapse = list(type = "sample", sample = 100, sample.fun = trait_normal), output = "combined.matrix")

cat("Post ordination ace completed...\n")


saveRDS(post_ord_ace, write.path("anc", "post_ord_ace_%03d.rds"))
saveRDS(ord_sample, write.path("ord", "ord_sample_%03d.rds"))
saveRDS(ord_point, write.path("ord", "ord_point_%03d.rds"))
saveRDS(ord_no_ace, write.path("ord", "ord_no_ace_%03d.rds"))
saveRDS(ord_true, write.path("ord", "ord_true_%03d.rds"))
saveRDS(point_post_ord_ace, write.path("ord", "post_ord_point_%03d.rds"))
saveRDS(sample_post_ord_ace, write.path("ord", "post_ord_sample_%03d.rds"))


cat("Finished replicate", replicate_id, "at", Sys.time(), "\n")
