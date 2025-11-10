args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
tree_size <- args[2]
library(parallel)
library(treats)
library(readr)

cat("Starting replicate", replicate_id, "\n")

set.seed(100 + replicate_id)

base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/discrete/", tree_size, "/")
job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")


write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}

# bd_params <- make.bd.params(speciation = 1, extinction = 0.7)


# stop_rule <- list(max.living = readr::parse_number(tree_size))
# # set.seed(123)
# tree <- treats(stop.rule = stop_rule, bd.params = bd_params, null.error = 100)

# b_d_est <- crude.bd.est(tree, "estimate")

# metadata_df <- data.frame(
#   replicate_id = replicate_id,
#   tree_size = length(tree$tip.label),
#   seed = 100 + replicate_id,
#   speciation = b_d_est$call$speciation,
#   extinction = b_d_est$call$extinction
# )
# metadata_df
# write.csv(metadata_df, 
#           write.path("metadata","metadata_%03d.csv"),
#           row.names = FALSE)


# ## write the saved trees
# write.tree(tree, write.path("trees", "tree_%03d.tre"))

tree <- read.tree(paste0("/mnt/parscratch/users/bip24cns/acedisparity/trees/overallDisparity/tree_", tree_size, sprintf("_%03d.tre", replicate_id)))

tree <- set.root.time(tree)

cat("Mapping traits...\n")


slow_binary_transitions <- matrix(c(
    0.99, 0.01,
    0.01, 0.99
), nrow = 2, byrow = TRUE)

slow_multi_transitions <- matrix(c(
    0.99, 0.005, 0.005,
    0.005, 0.99, 0.005,
    0.005, 0.005, 0.99
), nrow = 3, byrow = TRUE)

slow_binary <- treats::make.traits(process = discrete.process, n = 85, process.args = list(transitions = slow_binary_transitions))
slow_multi <- treats::make.traits(process = discrete.process, n = 15, process.args = list(transitions = slow_multi_transitions))


med_binary_transitions <- matrix(c(
  0.90, 0.1,
  0.1, 0.9
), nrow = 2, byrow = TRUE)


med_multi_transitions <- matrix(c(
  0.9, 0.05, 0.05,
  0.05, 0.9, 0.05,
  0.05, 0.05, 0.9
), nrow = 3, byrow = TRUE)

med_binary <- treats::make.traits(process = discrete.process, n = 85, process.args = list(transitions = med_binary_transitions))
med_multi <- treats::make.traits(process = discrete.process, n = 15, process.args = list(transitions = med_multi_transitions))

fast_binary_transitions <- matrix(c(
    0.1, 0.9,
    0.9, 0.1
), nrow = 2, byrow = TRUE)

fast_multi_transitions <- matrix(c(
    0.1, 0.45, 0.45, 
    0.45, 0.1, 0.45,
    0.45, 0.45, 0.1 
), nrow = 3, byrow = TRUE)

fast_binary <- treats::make.traits(process = discrete.process, n = 85, process.args = list(transitions = fast_binary_transitions))
fast_multi <- treats::make.traits(process = discrete.process, n = 15, process.args = list(transitions = fast_multi_transitions))



## map the traits

# Create list of trait objects
trait_sets <- list(
  slow = list(binary = slow_binary, multi = slow_multi),
  med = list(binary = med_binary, multi = med_multi),
  fast = list(binary = fast_binary, multi = fast_multi)
)

# Map all traits and combine
matrices <- lapply(trait_sets, function(traits) {
  mapped_binary <- (map.traits(traits$binary, tree))$data
  mapped_multi <- (map.traits(traits$multi, tree))$data
  cbind(mapped_binary, mapped_multi)
})

saveRDS(matrices, write.path("matrices", "matrices_%03d.rds"))

cat("Trait matrices saved for replicate", replicate_id, "\n")

source("/users/bip24cns/acedisparity/discrete/scripts/fossil.pres.R")
set_seed <- 100 + replicate_id
living <- lapply(matrices, remove.fossil, trees = tree, type = "discrete")
fossilised_high <- lapply(matrices, fossil.pres, trees = tree, preservation = 0.5, type = "discrete", seed = set_seed)
all_fossil <- lapply(matrices, fossil.pres, trees = tree, preservation = 1.0, type = "discrete", seed = set_seed)
fossilised_med <- lapply(matrices, fossil.pres, trees = tree, preservation = 0.15, type = "discrete", seed = set_seed)
fossilised_low <- lapply(matrices, fossil.pres, trees = tree, preservation = 0.05, type = "discrete", seed = set_seed)


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

fossil_trees <- lapply(fossil_matrices, lapply,  function(level){
  tree <- level$tree
})

cat("Fossil matrices created\n")
saveRDS(fossil_matrices, write.path("matrices", "fossil_matrices_%03d.rds"))



########################################################################################################################
## ANC STATES
tasks  <- expand.grid(rate = names(fossil_matrices), fossil_level = names(fossil_matrices[[1]]),  stringsAsFactors = FALSE) 
cl <- makeCluster(15)
clusterEvalQ(cl, library(treats))
clusterExport(cl, c("fossil_matrices", "tasks"))

# n_cores  <- (detectCores() - 1)

res_pre_ord_ace <- parLapply(cl, seq_len(nrow(tasks)), function(i){ ## loop over each model combination by row
    task <- tasks[i,]
    level <- fossil_matrices[[task$rate]][[task$fossil_level]]
    tryCatch({
    multi.ace(level$matrix, level$tree, models = "ER", output = "multi.ace")
  }, error = function(e) {
    cat("ERROR:", task$fossil_level, e$message, "\n") ## error handeling
    NULL
  })
})
stopCluster(cl)

pre_ord_ace <- list()
for(i in seq_along(res_pre_ord_ace)) {
  r <- tasks$rate[i]  
  l <- tasks$fossil_level[i]
  pre_ord_ace[[r]][[l]] <- res_pre_ord_ace[[i]]
}

# fossil_anc <- lapply(fossil_matrices, lapply, anc.states)

cat("Ancestral states estimated\n")

sample_fossil_anc <- lapply(pre_ord_ace, lapply,  multi.ace, sample = 100)

point_fossil_anc <- lapply(pre_ord_ace, lapply, multi.ace, threshold = FALSE, output = "combined.matrix", verbose = TRUE)

relative_fossil_anc <- lapply(pre_ord_ace, lapply,  multi.ace, output = "combined.matrix", verbose = TRUE)

saveRDS(pre_ord_ace, write.path("anc", "pre_ord_anc_%03d.rds"))
saveRDS(sample_fossil_anc, write.path("anc", "pre_ord_sample_%03d.rds"))
saveRDS(relative_fossil_anc, write.path("anc", "pre_ord_rel_%03d.rds"))
saveRDS(point_fossil_anc, write.path("anc", "pre_ord_point_%03d.rds"))

########################################################################################################################

## this function finds the most basal node label from the living sampling level, then extracts the tree from other fossil sampling levels from this node.
extract.living <- function(fossils) { 
  basal_node <- fossils$living$tree$node.label[1]
  living_nodes <- fossils$living$tree$node.label
  labels <- lapply(fossils, function(level){
    tree <- level$tree
    descendents <- (extract.clade(tree, node = basal_node))$tip.label
    add_nodes <- c(descendents, living_nodes)
  })
  return(labels)
}

labels <- lapply(fossil_matrices, extract.living)

point_living <- Map(function(rate_anc, rate_labels) {
    Map(function(fossil_anc, label_anc) {
      fossil_anc[label_anc, , drop = FALSE]
  }, rate_anc, rate_labels)
}, point_fossil_anc, labels)

names(point_living) <- names(point_fossil_anc)

rel_living <- Map(function(rate_anc, rate_labels) {
    Map(function(fossil_anc, label_anc) {
      fossil_anc[label_anc, , drop = FALSE]
  }, rate_anc, rate_labels)
}, relative_fossil_anc, labels)

names(rel_living) <- names(relative_fossil_anc)

sample_living <- Map(function(rate_anc, rate_labels) {
    Map(function(fossil_anc, label_anc) {
      lapply(fossil_anc, function(rep) {
        rep[label_anc, , drop = FALSE]
      })
  }, rate_anc, rate_labels)
}, sample_fossil_anc, labels)

names(sample_living) <- names(sample_fossil_anc)

no_ace_living <- lapply(point_living, lapply, function(matrix){
  no_node <- matrix[!grepl("^n", rownames(matrix)),]
})

true_living <-  lapply(seq_along(matrices), function(rate) {
  labs <- labels[[rate]]$all
  true_living <- matrices[[rate]][labs, , drop = FALSE]
})
names(true_living) <- names(matrices)


########################################################################################################################


tasks_sample_ord <- expand.grid(rate = names(sample_living), fossil_level = names(sample_living[[1]]), stringsAsFactors = FALSE)

cl <- makeCluster(15)
clusterEvalQ(cl, library(treats))
clusterExport(cl, c("sample_living",  "tasks_sample_ord"))

res_ord_sample <- parLapply(cl, seq_len(nrow(tasks_sample_ord)), function(i){
    task <- tasks_sample_ord[i, ]
    mat <- sample_living[[task$rate]][[task$fossil_level]]
    ord <- lapply(mat, function(rep){
        dist <- char.diff(rep, method = "mord", by.col = FALSE)
        ord_rep <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
        return(ord_rep)
    })
    return(ord)
})
stopCluster(cl)

ord_sample <- list()
for(i in seq_along(res_ord_sample)) {
  rate <- tasks_sample_ord$rate[i]
  fossil_level <- tasks_sample_ord$fossil_level[i]
  ord_sample[[rate]][[fossil_level]] <- res_ord_sample[[i]]
}


ord_rel <- lapply(rel_living, lapply, function(rep){
  dist <- char.diff(rep, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})

ord_point <- lapply(point_living, lapply,  function(rep){
  dist <- char.diff(rep, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})

ord_true <- lapply(true_living, function(rep){
  dist <- char.diff(rep, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})

ord_no_ace <- lapply(no_ace_living, lapply, function(rep){
  dist <- char.diff(rep, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})

saveRDS(ord_no_ace, write.path("ord", "ord_no_ace_%03d.rds"))
saveRDS(ord_true, write.path("ord", "ord_true_%03d.rds"))
saveRDS(ord_rel, write.path("ord", "ord_rel_%03d.rds"))
saveRDS(ord_point, write.path("ord", "ord_point_%03d.rds"))
saveRDS(ord_sample, write.path("ord", "ord_sample_%03d.rds"))

cat("ordinations calculated\n")

cat("Running post ord ace", replicate_id, "\n")

## post ordination ace


ord_fossil_tips <- lapply(fossil_matrices, lapply, function(x){
  mat <- x$matrix 
  dist <- char.diff(mat, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points ## ordinate the fossil tips from earlier
})

tasks_post_ord <- expand.grid(rate = names(ord_fossil_tips), fossil_level = names(ord_fossil_tips[[1]]), stringsAsFactors = FALSE)

cl <- makeCluster(15)
clusterEvalQ(cl, library(treats))
clusterExport(cl, c("ord_fossil_tips", "fossil_trees", "tasks_post_ord"))

res_post_ord <- parLapply(cl, seq_len(nrow(tasks_post_ord)), function(i){
    task <- tasks_post_ord[i,]
    ord_matrix <- ord_fossil_tips[[task$rate]][[task$fossil_level]]
    fossil_tree <- fossil_trees[[task$rate]][[task$fossil_level]]
    tryCatch({
        multi.ace(ord_matrix, fossil_tree, models = "BM", output = "multi.ace")
    }, error = function(e) {
        cat("ERROR in level:", task, e$message, "\n")
    })
})
stopCluster(cl)

post_ord_ace <- list()
for(i in seq_along(res_post_ord)) {
  r <- tasks_post_ord$rate[i]
  l <- tasks_post_ord$fossil_level[i]
  post_ord_ace[[r]][[l]] <- res_post_ord[[i]]
}


# post_ord_ace <- Map(function(rate_matrix, rate_tree){
#   Map(function(fossil_matrix, fossil_tree){
#     multi.ace(fossil_matrix, fossil_tree, models = "BM", output = "multi.ace")
#   }, rate_matrix, rate_tree)
# }, ord_fossil_tips, fossil_trees)

saveRDS(post_ord_ace, write.path("anc", "post_ord_ace_%03d.rds"))


# trait_normal = list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))
point_post_ord_ace <- lapply(post_ord_ace, lapply,  multi.ace, output = "combined.matrix")

point_post_ord_ace_living <- Map(function(rate_anc, rate_labels) {
    Map(function(fossil_anc, label_anc) {
      fossil_anc[label_anc, , drop = FALSE]
  }, rate_anc, rate_labels)
}, point_post_ord_ace, labels)

names(point_post_ord_ace_living) <- names(point_post_ord_ace)

trait_normal <- list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))

sample_post_ord_ace <- lapply(post_ord_ace, lapply,  multi.ace, sample = 100, sample.fun = trait_normal, output = "combined.matrix")

sample_post_ord_ace_living <- Map(function(rate_anc, rate_labels) {
    Map(function(fossil_anc, label_anc) {
      lapply(fossil_anc, function(rep) {
        rep[label_anc, , drop = FALSE]
      })
  }, rate_anc, rate_labels)
}, sample_post_ord_ace, labels)
names(sample_post_ord_ace_living) <- names(sample_post_ord_ace)


cat("Finished post ord ace", replicate_id, "\n")


saveRDS(point_post_ord_ace_living, write.path("ord", "post_ord_point_%03d.rds"))
saveRDS(sample_post_ord_ace_living, write.path("ord", "post_ord_sample_%03d.rds"))

cat("Finished replicate", replicate_id, "\n")
