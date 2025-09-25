args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
# library(parallel)
library(treats)

cat("Starting replicate", replicate_id, "\n")

set.seed(100 + replicate_id)

bd_params <- make.bd.params(speciation = 1, extinction = 0.7)

stop_rule <- list(max.living = 50)

tree <- treats(stop.rule = stop_rule, bd.params = bd_params, null.error = 100)
tree <- drop.singles(tree)
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

cat("=== STARTING fossil_anc ===\n"); flush.console()
withCallingHandlers({
  fossil_anc <- lapply(fossil_matrices, lapply, anc.states)
}, warning = function(w) {
  cat("WARNING in fossil_anc:", w$message, "\n")
})

cat("=== STARTING point_anc ===\n"); flush.console()
withCallingHandlers({
  point_anc <- lapply(fossil_anc, lapply, multi.ace, output = "combined.matrix")
}, warning = function(w) {
  cat("WARNING in point_anc:", w$message, "\n")
})

trait_normal = list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))

cat("=== STARTING sample_anc ===\n"); flush.console()
withCallingHandlers({
  sample_anc <- lapply(fossil_anc, lapply, multi.ace, output = "combined.matrix", sample = 100, sample.fun = trait_normal)
}, warning = function(w) {
  cat("WARNING in sample_anc:", w$message, "\n")
})

cat("Ancestral states estimated...\n")

saveRDS(fossil_anc, paste0(file_path, sprintf("anc/fossil_anc_%03d.rds", replicate_id)))
saveRDS(point_anc, paste0(file_path, sprintf("anc/point_anc_%03d.rds", replicate_id)))
saveRDS(sample_anc, paste0(file_path, sprintf("anc/sample_anc_%03d.rds", replicate_id)))
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



cat("=== STARTING ord_sample ===\n"); flush.console()
withCallingHandlers({
  ord_sample <- lapply(sample_living, lapply, lapply, function(mat){
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })
}, warning = function(w) {
  cat("WARNING in ord_sample:", w$message, "\n")
})

cat("=== STARTING ord_point ===\n"); flush.console()
withCallingHandlers({
  ord_point <- lapply(point_living, lapply, function(mat){
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })
}, warning = function(w) {
  cat("WARNING in ord_point:", w$message, "\n")
})

cat("=== STARTING ord_no_ace ===\n"); flush.console()
withCallingHandlers({
  ord_no_ace <- lapply(no_ace_living, lapply, function(mat){
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })
}, warning = function(w) {
  cat("WARNING in ord_no_ace:", w$message, "\n")
})

cat("=== STARTING ord_true ===\n"); flush.console()
withCallingHandlers({
  ord_true <- lapply(true_living,  function(mat){
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })
}, warning = function(w) {
  cat("WARNING in ord_true:", w$message, "\n")
})

cat("Ordinations completed...\n")


cat("=== STARTING ord_fossil_tips ===\n"); flush.console()
withCallingHandlers({
  ord_fossil_tips <- lapply(fossil_matrices, lapply, function(x){
    mat <- x$matrix
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })
}, warning = function(w) {
  cat("WARNING in ord_fossil_tips:", w$message, "\n")
})



cat("=== STARTING post_ord_ace ===\n"); flush.console()
withCallingHandlers({
  post_ord_ace <- Map(function(rate_matrix, rate_tree){
    Map(function(fossil_matrix, fossil_tree){
      multi.ace(fossil_matrix, fossil_tree, models = "ML", output = "multi.ace")
    }, rate_matrix, rate_tree)
  }, ord_fossil_tips, fossil_trees)
}, warning = function(w) {
  cat("WARNING in post_ord_ace:", w$message, "\n")
})

cat("=== STARTING point_post_ord_ace ===\n"); flush.console()
withCallingHandlers({
  point_post_ord_ace <- lapply(post_ord_ace, lapply,  multi.ace, output = "combined.matrix")
}, warning = function(w) {
  cat("WARNING in point_post_ord_ace:", w$message, "\n")
})

point_post_ord_ace_living <- Map(function(rate_anc, rate_labels) {
    Map(function(fossil_anc, label_anc) {
      fossil_anc[label_anc, , drop = FALSE]
  }, rate_anc, rate_labels)
}, point_post_ord_ace, labels)

names(point_post_ord_ace_living) <- names(point_post_ord_ace)

trait_normal = list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))

cat("=== STARTING sample_post_ord_ace ===\n"); flush.console()
withCallingHandlers({
  sample_post_ord_ace <- lapply(post_ord_ace, lapply,  multi.ace, sample = 100, sample.fun = trait_normal, output = "combined.matrix")
}, warning = function(w) {
  cat("WARNING in sample_post_ord_ace:", w$message, "\n")
})

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

cat("=== CHECKING WARNINGS ===\n")
warning_list <- warnings()
if(!is.null(warning_list)) {
  cat("Number of warnings:", length(warning_list), "\n")
  print(warning_list)
} else {
  cat("No warnings detected\n")
}
cat("Script completed\n")