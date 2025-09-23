args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
# library(parallel)
library(treats)

cat("Starting replicate", replicate_id, "\n")

set.seed(100 + replicate_id)

bd_params <- make.bd.params(speciation = 1, extinction = 0.7)

stop_rule <- list(max.living = 50)

tree <- treats(stop.rule = stop_rule, bd.params = bd_params, null.error = 100)

cat("Tree simulated...\n")


file_path <- c("/mnt/parscratch/users/bip24cns/acedisparity/continuous/")

tree_height <- max(node.depth.edgelength(tree))
## write the saved trees
write.tree(tree, paste0(file_path, sprintf("trees/continuous_tree_%03d.tre", replicate_id)))


bm <- make.traits(process = BM.process, n = 100)
ou_strong <- make.traits(process = OU.process, n = 100, process.args = list(alpha = (log(2) / (tree_height / 10))))
ou_weak <- make.traits(process = OU.process, n = 100, process.args = list(alpha = (log(2) / tree_height)))



# bm_switch <- make.events(target = "traits", condition = age.condition(tree_height / 2), modification = traits.update(process = BM.process))
bm_matrices <- map.traits(bm, tree)$data
ou_s_matrices <- map.traits(ou_strong, tree)$data
ou_w_matrices <- map.traits(ou_weak, tree)$data
# ou_s_bm_matrices <- map.traits(ou_strong, tree, events = bm_switch)
cat("Traits simulated...\n")

matrices <- list(bm = bm_matrices, ou_w = ou_w_matrices, ou_s = ou_s_matrices)

saveRDS(matrices, paste0(file_path, sprintf("matrices/matrices_%03d.rds", replicate_id)))

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

saveRDS(fossil_matrices, paste0(file_path, sprintf("matrices/fossil_matrices_%03d.rds", replicate_id)))
saveRDS(fossil_trees, paste0(file_path, sprintf("trees/fossil_trees_%03d.rds", replicate_id)))


## ANC STATES
# n_cores  <- 18
anc.states <- function(x) {
  # Run multi.ace for each tree
  anc_states <- multi.ace(data = x$matrix, 
                          tree = x$tree, 
                          models = "ML", 
                          output = "multi.ace"
                          # parallel = n_cores,
                          # verbose = TRUE
                          )
return(anc_states)}

fossil_anc <- lapply(fossil_matrices, lapply, anc.states)

point_anc <- lapply(fossil_anc, lapply, multi.ace, output = "combined.matrix")
trait_normal = list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))

sample_anc <- lapply(fossil_anc, lapply, multi.ace, output = "combined.matrix", sample = 100, sample.fun = trait_normal)
cat("Ancestral states estimated...\n")

saveRDS(fossil_anc, paste0(file_path, sprintf("anc/fossil_anc_%03d.rds", replicate_id)))
saveRDS(point_anc, paste0(file_path, sprintf("anc/point_anc_%03d.rds", replicate_id)))
saveRDS(sample_anc, paste0(file_path, sprintf("anc/sample_anc_%03d.rds", replicate_id)))
cat("Ancestral states saved...\n")


# sample_anc <- readRDS("../Data/cluster/continuous/anc/sample_anc_005.rds")
# point_anc <- readRDS("../Data/cluster/continuous/anc/point_anc_005.rds")

# matrices <- readRDS("../Data/cluster/continuous/matrices/matrices_005.rds")
# fossil_matrices <- readRDS("../Data/cluster/continuous/matrices/fossil_matrices_005.rds")
# tree <- read.tree("../Data/cluster/continuous/trees/continuous_tree_005.tre")
# fossil_trees <- readRDS("../Data/cluster/continuous/trees/fossil_trees_005.rds")

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



ord_sample <- lapply(sample_living, lapply, lapply, function(mat){
  prcomp(mat, scale = FALSE, center = TRUE)$x
})

ord_point <- lapply(point_living, lapply, function(mat){
  prcomp(mat, scale = FALSE, center = TRUE)$x
})

ord_no_ace <- lapply(no_ace_living, lapply, function(mat){
  prcomp(mat, scale = FALSE, center = TRUE)$x
})

ord_true <- lapply(true_living,  function(mat){
  prcomp(mat, scale = FALSE, center = TRUE)$x
})

cat("Ordinations completed...\n")


ord_fossil_tips <- lapply(fossil_matrices, lapply, function(x){
  mat <- x$matrix
  prcomp(mat, scale = FALSE, center = TRUE)$x
})



post_ord_ace <- Map(function(rate_matrix, rate_tree){
  Map(function(fossil_matrix, fossil_tree){
    multi.ace(fossil_matrix, fossil_tree, models = "ML", output = "multi.ace")
  }, rate_matrix, rate_tree)
}, ord_fossil_tips, fossil_trees)

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

saveRDS(ord_sample, paste0(file_path, sprintf("ord/sample_ord_%03d.rds", replicate_id)))
saveRDS(ord_point, paste0(file_path, sprintf("ord/point_ord_%03d.rds", replicate_id)))
saveRDS(ord_no_ace, paste0(file_path, sprintf("ord/no_ace_ord_%03d.rds", replicate_id)))
saveRDS(ord_true, paste0(file_path, sprintf("ord/true_ord_%03d.rds", replicate_id)))
saveRDS(point_post_ord_ace_living, paste0(file_path, sprintf("ord/point_post_ord_%03d.rds", replicate_id)))
saveRDS(sample_post_ord_ace_living, paste0(file_path, sprintf("ord/sample_post_ord_%03d.rds", replicate_id)))


metrics_list <- list(
  sum_var = c(sum, variances),
  sum_quant = c(sum, quantiles),
  mean_neigh = c(mean, neighbours)
)

# Calculate disparity for all metrics
for(metric_name in names(metrics_list)) {
  metric <- metrics_list[[metric_name]]
  
  # Apply to each ordination type with appropriate structure
  assign(paste0(metric_name, "_sample"), 
         lapply(ord_sample, lapply, lapply, function(mat) dispRity(mat, metric = metric)$disparity))
  
  assign(paste0(metric_name, "_point"), 
         lapply(ord_point, lapply, function(mat) dispRity(mat, metric = metric)$disparity))
  
  assign(paste0(metric_name, "_no_ace"), 
         lapply(ord_no_ace, lapply, function(mat) dispRity(mat, metric = metric)$disparity))
  
  assign(paste0(metric_name, "_true"), 
         lapply(ord_true, function(mat) dispRity(mat, metric = metric)$disparity))
  
  assign(paste0(metric_name, "_point_postord"), 
         lapply(point_post_ord_ace_living, lapply, function(mat) dispRity(mat, metric = metric)$disparity))

  assign(paste0(metric_name, "_sample_postord"), 
         lapply(sample_post_ord_ace_living, lapply, lapply, function(mat) dispRity(mat, metric = metric)$disparity))
}

cat("Disparity calculation finished...\n")

for(metric_name in names(metrics_list)) {
  for(ord_type in c("sample", "point", "point_postord", "sample_postord", "no_ace", "true")) {
    var_name <- paste0(metric_name, "_", ord_type)
    saveRDS(get(var_name), paste0(file_path, sprintf("disparity/%s_%03d.rds", var_name, replicate_id)))
  }
}

cat("Disparity files saved...\n")



est_disp <- list(
  point = mget(paste0(c("sum_var", "sum_quant", "mean_neigh"), "_point")),
  sample = mget(paste0(c("sum_var", "sum_quant", "mean_neigh"), "_sample")),
  no_ace = mget(paste0(c("sum_var", "sum_quant", "mean_neigh"), "_no_ace")),
  point_postord = mget(paste0(c("sum_var", "sum_quant", "mean_neigh"), "_point_postord")),
  sample_postord = mget(paste0(c("sum_var", "sum_quant", "mean_neigh"), "_sample_postord"))
)

true_disp <- mget(paste0(c("sum_var", "sum_quant", "mean_neigh"), "_true"))



diffs <- lapply(est_disp, function(method){
  method_result <- Map(function(est_metric, true_metric) {
    Map(function(est_model, true_model){
      true_disp <- true_model[[1]]$elements[1]
      lapply(est_model, function(fossil_lev){
        if (is.list(fossil_lev) && !is.null(fossil_lev[[1]]$elements)) {
        diff <- (fossil_lev[[1]]$elements[1] - true_disp) / true_disp ## calculate relative disparity difference 
        return(diff)
        } else {lapply(fossil_lev, function(rep) { ## for sample methods with additional nested list
          diff <- (rep[[1]]$elements[1] - true_disp) / true_disp 
          return(diff)
          })
        }
      })
    }, est_metric, true_metric)
  }, method, true_disp)
})

cat("Disparity diffs calculated...\n")


methods <- names(diffs)
for (method in methods){
  saveRDS(diffs[[method]], paste0(file_path, sprintf("disparity/diff_%s_%03d.rds", method, replicate_id)))
}

cat("Disparity diffs saved...\n")
