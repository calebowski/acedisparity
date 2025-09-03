args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
# library(parallel)
library(treats)

cat("Starting replicate", replicate_id, "\n")

set.seed(100 + replicate_id)

bd_params <- make.bd.params(speciation = 1, extinction = 0.7)

stop_rule <- list(max.living = 50)
# set.seed(123)
tree <- treats(stop.rule = stop_rule, bd.params = bd_params, null.error = 100)

## write the saved trees
write.tree(tree, sprintf("/users/bip24cns/acedisparity/discrete/out/discrete_tree_%03d.tre", replicate_id))

cat("Mapping traits...\n")


slow_binary_transitions <- matrix(c(
    0.99, 0.01,
    0.01, 0.99
), nrow = 2, byrow = TRUE)

slow_multi_transitions <- matrix(c(
    0.99, 0.025, 0.025,
    0.025, 0.99, 0.025,
    0.025, 0.025, 0.99
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

saveRDS(matrices, sprintf("/users/bip24cns/acedisparity/discrete/out/matrices_%03d.rds", replicate_id))

cat("Trait matrices saved for replicate", replicate_id, "\n")

source("/users/bip24cns/acedisparity/discrete/scripts/fossil.pres.R")

living <- lapply(matrices, remove.fossil, trees = tree, type = "discrete")
fossilised_high <- lapply(matrices, fossil.pres, trees = tree, preservation = 0.5, type = "discrete")
all_fossil <- lapply(matrices, fossil.pres, trees = tree, preservation = 1.0, type = "discrete")
fossilised_med <- lapply(matrices, fossil.pres, trees = tree, preservation = 0.15, type = "discrete")
fossilised_low <- lapply(matrices, fossil.pres, trees = tree, preservation = 0.05, type = "discrete")



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

cat("Fossil matrices created\n")
saveRDS(fossil_matrices, sprintf("/users/bip24cns/acedisparity/discrete/out/fossil_matrices_%03d.rds", replicate_id))



########################################################################################################################
## ANC STATES

# n_cores  <- (detectCores() - 1)
anc.states <- function(x) {
  # Run multi.ace for each tree
  anc_states <- multi.ace(data = x$matrix, 
                          tree = x$tree, 
                          models = "SYM", 
                          output = "multi.ace")
return(anc_states)}

fossil_anc <- lapply(fossil_matrices, lapply, anc.states)

cat("Ancestral states estimated\n")

sample_fossil_anc <- lapply(fossil_anc, lapply,  multi.ace, sample = 100)

strict_fossil_anc <- lapply(fossil_anc, lapply, multi.ace, threshold = FALSE, output = "combined.matrix", verbose = TRUE)

relative_fossil_anc <- lapply(fossil_anc, lapply,  multi.ace, output = "combined.matrix", verbose = TRUE)

saveRDS(fossil_anc, sprintf("/users/bip24cns/acedisparity/discrete/out/discrete_anc_%03d.rds", replicate_id))
saveRDS(sample_fossil_anc, sprintf("/users/bip24cns/acedisparity/discrete/out/sample_anc_%03d.rds", replicate_id))
saveRDS(relative_fossil_anc, sprintf("/users/bip24cns/acedisparity/discrete/out/rel_anc_%03d.rds", replicate_id))
saveRDS(strict_fossil_anc, sprintf("/users/bip24cns/acedisparity/discrete/out/strict_anc_%03d.rds", replicate_id))



########################################################################################################################


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


strict_living <- Map(function(rate_anc, rate_labels) {
    Map(function(fossil_anc, label_anc) {
      fossil_anc[label_anc, , drop = FALSE]
  }, rate_anc, rate_labels)
}, strict_fossil_anc, labels)

names(strict_living) <- names(strict_fossil_anc)

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

no_ace_living <- lapply(strict_living, lapply, function(matrix){
  no_node <- matrix[!grepl("^n", rownames(matrix)),]
})

true_living <-  lapply(seq_along(matrices), function(rate) {
  labs <- labels[[rate]]$all
  true_living <- matrices[[rate]][labs, , drop = FALSE]
})
names(true_living) <- names(matrices)


########################################################################################################################

ord_sample <- lapply(sample_living, lapply, lapply, function(rep){
  dist <- char.diff(rep, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})

saveRDS(ord_sample, sprintf("/users/bip24cns/acedisparity/discrete/out/ord_sample_%03d.rds", replicate_id))


ord_rel <- lapply(rel_living, lapply, function(rep){
  dist <- char.diff(rep, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})

saveRDS(ord_rel, sprintf("/users/bip24cns/acedisparity/discrete/out/ord_rel_%03d.rds", replicate_id))


ord_strict <- lapply(strict_living, lapply,  function(rep){
  dist <- char.diff(rep, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})

saveRDS(ord_strict, sprintf("/users/bip24cns/acedisparity/discrete/out/ord_strict_%03d.rds", replicate_id))


ord_true <- lapply(true_living, function(rep){
  dist <- char.diff(rep, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})

saveRDS(ord_true, sprintf("/users/bip24cns/acedisparity/discrete/out/ord_true_%03d.rds", replicate_id))


ord_no_ace <- lapply(no_ace_living, lapply, function(rep){
  dist <- char.diff(rep, method = "mord", by.col = FALSE)
  ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})

saveRDS(ord_no_ace, sprintf("/users/bip24cns/acedisparity/discrete/out/ord_no_ace_%03d.rds", replicate_id))

cat("ordinations calculated\n")

cat("Finished replicate", replicate_id, "\n")