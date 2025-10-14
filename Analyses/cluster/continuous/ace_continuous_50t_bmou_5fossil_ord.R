args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
# library(parallel)
library(treats)
library(parallel)
library(MASS)

cat("Starting replicate", replicate_id, "\n")

set.seed(100 + replicate_id)

write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}
base_path <- "/mnt/parscratch/users/bip24cns/acedisparity/continuous/50t/"
job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")

bd_params <- make.bd.params(speciation = 1, extinction = 0.7)

stop_rule <- list(max.living = 50)


tree <- treats(stop.rule = stop_rule, bd.params = bd_params, null.error = 100)
tree <- drop.singles(tree)
cat("Tree simulated...\n")
b_d_est <- crude.bd.est(tree, "estimate") # estimate the birth death params
write.tree(tree, write.path("trees", "tree_%03d.tre")) ## write the saved tre


metadata_df <- data.frame(
  replicate_id = replicate_id,
  tree_size = length(tree$tip.label),
  seed = 100 + replicate_id,
  speciation = b_d_est$call$speciation,
  extinction = b_d_est$call$extinction
)


write.csv(metadata_df, 
          write.path("metadata","metadata_%03d.csv"),
          row.names = FALSE)

tree_height <- max(node.depth.edgelength(tree))

BM.trend.process <- function(x0 = 0, edge.length = 1, Sigma = diag(length(x0)), trend = 0.1, ...) {
      # Square root gives more trend than log but less than linear
      drift <- trend * log(edge.length + 1)
      if(edge.length < (0.01 * Sigma[1,1])){drift <- 0} ## if edge length is smaller than 1% of variance, drift is 0.
      return(t(MASS::mvrnorm(n = 1, mu = x0 + drift, Sigma = Sigma * edge.length, ...)))
}




# Sigma_matrix <- diag(0.25, 100) ## make diagonal variance-covariance matrix
bm <- make.traits(process = BM.process, n = 100, process.args = list(Sigma = diag(0.25, 100)))
bm_trend <- make.traits(process = BM.trend.process, n = 100, process.args = list(Sigma = diag(0.25, 100), trend = 0.3))
ou_strong <- make.traits(process = OU.process, n = 100, process.args = list(alpha = (log(2) / (tree_height / 10)), Sigma = diag(0.25, 100)))
ou_weak <- make.traits(process = OU.process, n = 100, process.args = list(alpha = (log(2) / 5), Sigma = diag(0.25, 100)))
ou_shift <- make.traits(process = OU.process, n = 100, process.args = list(optimum = 2, alpha = (log(2) / (tree_height / 5)), Sigma = diag(0.25, 100)))




# bm_switch <- make.events(target = "traits", condition = age.condition(tree_height / 2), modification = traits.update(process = BM.process))
bm_matrices <- map.traits(bm, tree)$data
bm_matrices_trend <- map.traits(bm_trend, tree)$data
ou_s_matrices <- map.traits(ou_strong, tree)$data
ou_w_matrices <- map.traits(ou_weak, tree)$data
ou_shift_matrices <- map.traits(ou_shift, tree)$data

# ou_s_bm_matrices <- map.traits(ou_strong, tree, events = bm_switch)
cat("Traits simulated...\n")

matrices <- list(bm = bm_matrices, bm_t = bm_matrices_trend, ou_w = ou_w_matrices, ou_st = ou_s_matrices, ou_sh = ou_shift_matrices)

saveRDS(matrices, write.path("matrices", "matrices_%03d.rds"))

cat("Trait matrices saved for replicate", replicate_id, "\n")

source("/users/bip24cns/acedisparity/discrete/scripts/fossil.pres.R")


living <- lapply(matrices, remove.fossil, trees = tree, type = "continuous")
fossilised_high <- lapply(matrices, fossil.pres, trees = tree, preservation = 0.5, type = "continuous")
all_fossil <- lapply(matrices, fossil.pres, trees = tree, preservation = 1.0, type = "continuous")
fossilised_med <- lapply(matrices, fossil.pres, trees = tree, preservation = 0.15, type = "continuous")
fossilised_low <- lapply(matrices, fossil.pres, trees = tree, preservation = 0.05, type = "continuous")




fossil_matrices <- lapply(names(matrices), function(level) {
    list(
      all = all_fossil[[level]],
      fossil_high = fossilised_high[[level]],
      fossil_med = fossilised_med[[level]],
      fossil_low = fossilised_low[[level]],
      living = living[[level]]
    )
})
names(fossil_matrices) <- names(matrices)

fossil_trees <- lapply(fossil_matrices, lapply,  function(level){
  tree <- level$tree
})
cat("Fossil matrices completed...\n")

saveRDS(fossil_matrices, write.path("matrices", "fossil_matrices_%03d.rds"))
saveRDS(fossil_trees, write.path("trees", "fossil_trees_%03d.rds"))


## ANC STATES
# n_cores  <- 18
# anc.states <- function(x) {
#   # Run multi.ace for each tree
#   anc_states <- multi.ace(data = x$matrix, 
#                           tree = x$tree, 
#                           models = "ML", 
#                           output = "multi.ace"
#                           # parallel = n_cores,
#                           # verbose = TRUE
#                           )
# return(anc_states)}

## begin cluster
# Parallel across both models AND fossil levels
library(foreach)
library(doParallel)

cl <- makeCluster(10)  # Moderate number
registerDoParallel(cl)
clusterEvalQ(cl, library(treats))
clusterExport(cl, "fossil_matrices")

fossil_anc <- foreach(model_name = names(fossil_matrices), .packages = "treats") %dopar% {
  fossil_level <- fossil_matrices[[model_name]]
  fossil_names <- names(fossil_level)
  
  output <- lapply(fossil_level, function(x) {
    tryCatch({
      multi.ace(x$matrix, x$tree, models = "BM", output = "multi.ace")
    }, error = function(e) {
      cat("ERROR in multi.ace output:", e$message, "\n")
      return(NULL)
    })
  })
  names(output) <- fossil_names
  return(output)
}
names(fossil_anc) <- names(fossil_matrices)

stopCluster(cl)

cat("=== STARTING point_anc ===\n")
  point_anc <- lapply(fossil_anc, lapply, multi.ace, output = "combined.matrix")

trait_normal  <-  list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))


sample_anc <- lapply(fossil_anc, lapply, multi.ace, output = "combined.matrix", sample = 100, sample.fun = trait_normal)


cat("Ancestral states estimated...\n")

saveRDS(fossil_anc, write.path("anc", "fossil_anc_%03d.rds"))
saveRDS(point_anc, write.path("anc", "point_anc_%03d.rds"))
saveRDS(sample_anc, write.path("anc", "sample_anc_%03d.rds"))
cat("Ancestral states saved...\n")


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
}, point_anc, labels)
names(point_living) <- names(point_anc)


sample_living <- Map(function(rate_anc, rate_labels) {
    Map(function(fossil_anc, label_anc) {
      lapply(fossil_anc, function(rep) {
        rep[label_anc, , drop = FALSE]
      })
  }, rate_anc, rate_labels)
}, sample_anc, labels)

names(sample_living) <- names(sample_anc)



no_ace_living <- lapply(point_living, lapply, function(matrix){
  no_node <- matrix[!grepl("^n", rownames(matrix)),]
})


true_living <-  lapply(seq_along(matrices), function(rate) {
  labs <- labels[[rate]]$all
  true_living <- matrices[[rate]][labs, , drop = FALSE]
})
names(true_living) <- names(matrices)

cat("Living matrices completed...\n")



cat("=== STARTING ord_sample ===\n"); 
  ord_sample <- lapply(sample_living, lapply, lapply, function(mat){
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })


cat("=== STARTING ord_point ===\n")
  ord_point <- lapply(point_living, lapply, function(mat){
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })


cat("=== STARTING ord_no_ace ===\n")
  ord_no_ace <- lapply(no_ace_living, lapply, function(mat){
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })


cat("=== STARTING ord_true ===\n"); 
  ord_true <- lapply(true_living,  function(mat){
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })


cat("Ordinations completed...\n")


cat("=== STARTING ord_fossil_tips ===\n")
ord_fossil_tips <- lapply(fossil_matrices, lapply, function(x){
    mat <- x$matrix
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })

cat("=== STARTING post_ord_ace ===\n")
# post_ord_ace <- Map(function(rate_matrix, rate_tree){
#     Map(function(fossil_matrix, fossil_tree){
#       multi.ace(fossil_matrix, fossil_tree, models = "ML", output = "multi.ace")
#     }, rate_matrix, rate_tree)
#   }, ord_fossil_tips, fossil_trees)



cl <- makeCluster(10)
registerDoParallel(cl)
clusterEvalQ(cl, library(treats))
clusterExport(cl, c("ord_fossil_tips", "fossil_trees"))

post_ord_ace <- foreach(model_name = names(ord_fossil_tips), .packages = "treats") %dopar% {
  rate_matrix <- ord_fossil_tips[[model_name]]
  rate_tree <- fossil_trees[[model_name]]
  
  # Map across fossil levels for this model
  Map(function(fossil_matrix, fossil_tree) {
    multi.ace(fossil_matrix, fossil_tree, models = "BM", output = "multi.ace")
  }, rate_matrix, rate_tree)
}
names(post_ord_ace) <- names(ord_fossil_tips)

stopCluster(cl)

cat("=== STARTING point_post_ord_ace ===\n");
point_post_ord_ace <- lapply(post_ord_ace, lapply,  multi.ace, output = "combined.matrix")


point_post_ord_ace_living <- Map(function(rate_anc, rate_labels) {
    Map(function(fossil_anc, label_anc) {
      fossil_anc[label_anc, , drop = FALSE]
  }, rate_anc, rate_labels)
}, point_post_ord_ace, labels)

names(point_post_ord_ace_living) <- names(point_post_ord_ace)

trait_normal = list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))


sample_post_ord_ace <- lapply(post_ord_ace, lapply,  multi.ace, sample = 100, sample.fun = trait_normal, output = "combined.matrix")


sample_post_ord_ace_living <- Map(function(rate_anc, rate_labels) {
    Map(function(fossil_anc, label_anc) {
      lapply(fossil_anc, function(rep) {
        rep[label_anc, , drop = FALSE]
      })
  }, rate_anc, rate_labels)
}, sample_post_ord_ace, labels)
cat("Post ordination ace completed...\n")

saveRDS(ord_sample, write.path("ord", "ord_sample_%03d.rds"))
saveRDS(ord_point, write.path("ord", "ord_point_%03d.rds"))
saveRDS(ord_no_ace, write.path("ord", "ord_no_ace_%03d.rds"))
saveRDS(ord_true, write.path("ord", "ord_true_%03d.rds"))
saveRDS(point_post_ord_ace_living, write.path("ord", "post_ord_point_%03d.rds"))
saveRDS(sample_post_ord_ace_living, write.path("ord", "post_ord_sample_%03d.rds"))

# cat("=== CHECKING WARNINGS ===\n")
# warning_list <- warnings()
# if(!is.null(warning_list)) {
#   cat("Number of warnings:", length(warning_list), "\n")
#   print(warning_list)
# } else {
#   cat("No warnings detected\n")
# }
cat("Script completed\n")