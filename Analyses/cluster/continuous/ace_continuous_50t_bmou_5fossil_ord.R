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
job_id <- 8368504

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
ou_weak <- make.traits(process = OU.process, n = 100, process.args = list(alpha = (log(2) / (tree_height / 2)), Sigma = diag(0.25, 100)))
ou_shift <- make.traits(process = OU.process, n = 100, process.args = list(optimum = 2, alpha = (log(2) / (tree_height / 2)), Sigma = diag(0.25, 100)))


bm_matrices <- map.traits(bm, tree)$data
bm_matrices_trend <- map.traits(bm_trend, tree)$data
ou_s_matrices <- map.traits(ou_strong, tree)$data
ou_w_matrices <- map.traits(ou_weak, tree)$data
ou_shift_matrices <- map.traits(ou_shift, tree)$data

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


tasks  <- expand.grid(model = names(fossil_matrices), fossil_level = names(fossil_matrices[[1]]), stringsAsFactors = FALSE) ## flatten to give 25 rows, model combinations for each row

cl <- makeCluster(25) ## make 25 core cluster (one for each task)
clusterEvalQ(cl, library(treats))
clusterExport(cl, c("fossil_matrices", "tasks"))


res <- parLapply(cl, seq_len(nrow(tasks)), function(i){ ## loop over each model combination by row
    task <- tasks[i,]
    model_combination <- fossil_matrices[[task$model]][[task$fossil_level]]
    tryCatch({
    multi.ace(model_combination$matrix, model_combination$tree, models = "BM", output = "multi.ace")
  }, error = function(e) {
    cat("ERROR:", task$model, task$fossil_level, e$message, "\n") ## error handeling
    NULL
  })
})

stopCluster(cl)

fossil_anc <- list(bm = list(), bm_t = list(), ou_w = list(), ou_st = list(), ou_sh = list())
for(i in seq_along(res)) {
  m <- tasks$model[i]
  l <- tasks$fossil_level[i]
  fossil_anc[[m]][[l]] <- res[[i]]
}

cat("Starting point ancestral state estimation...\n")
point_anc <- lapply(fossil_anc, lapply, multi.ace, output = "combined.matrix")  ## takes median of confidence interval

trait_normal  <-  list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4))) ## samples with normal distribution

cat("Starting distribution ancestral state estimation...\n")
sample_anc <- lapply(fossil_anc, lapply, multi.ace, output = "combined.matrix", sample = 100, sample.fun = trait_normal) ## 100 matrices replicated, 


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


cat("=== STARTING ordination of fossil tips for post ordiation ace ===\n")
ord_fossil_tips <- lapply(fossil_matrices, lapply, function(x) {
    mat <- x$matrix
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })

cat("=== STARTING ancestral statte estimation of post ordination ===\n")

tasks_post_ord <- expand.grid(model = names(ord_fossil_tips), fossil_level = names(ord_fossil_tips[[1]]), stringsAsFactors = FALSE)

cl <- makeCluster(25)
clusterEvalQ(cl, library(treats))
clusterExport(cl, c("ord_fossil_tips", "fossil_trees", "tasks_post_ord"))

res_post_ord <- parLapply(cl, seq_len(nrow(tasks_post_ord)), function(i){
    task <- tasks_post_ord[i,]
    ord_matrix <- ord_fossil_tips[[task$model]][[task$fossil_level]]
    fossil_tree <- fossil_trees[[task$model]][[task$fossil_level]]
    tryCatch({
        multi.ace(ord_matrix, fossil_tree, models = "BM", output = "multi.ace")
    }, error = function(e) {
        cat("ERROR:", task$model, task$fossil_level, e$message, "\n")
    })
})
stopCluster(cl)

post_ord_ace  <- list(bm = list(), bm_t = list(), ou_w = list(), ou_st = list(), ou_sh = list())
for(i in seq_along(res_post_ord)){
    m <- tasks_post_ord$model[i]
    l <- tasks_post_ord$fossil_level[i]
    post_ord_ace[[m]][[l]] <- res_post_ord[[i]]
}

cat("=== STARTING point estimation post ord ace ===\n");
point_post_ord_ace <- lapply(post_ord_ace, lapply,  multi.ace, output = "combined.matrix")


point_post_ord_ace_living <- Map(function(rate_anc, rate_labels) {
    Map(function(fossil_anc, label_anc) {
      fossil_anc[label_anc, , drop = FALSE]
  }, rate_anc, rate_labels)
}, point_post_ord_ace, labels)

names(point_post_ord_ace_living) <- names(point_post_ord_ace)

trait_normal <-  list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))


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


cat("Script completed\n")