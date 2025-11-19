
replicate_id <- 60
library(treats)
library(ape)

base_path <- "/mnt/parscratch/users/bip24cns/acedisparity/discrete/100t/"
job_id <- 8190277

write.path <- function(subfolder, filename) {
  paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id))
}

cat("Resuming post-ordination ACE for replicate", replicate_id, "\n")

# Load previously saved fossil matrices
fossil_matrices <- readRDS(write.path("matrices", "fossil_matrices_%03d.rds"))

# Recompute fossil_trees structure (same as original script)
fossil_trees <- lapply(fossil_matrices, lapply, function(level) level$tree)

# re-create the extract.living helper (same logic as original run)
extract.living <- function(fossils) {
  basal_node <- fossils$living$tree$node.label[1]
  living_nodes <- fossils$living$tree$node.label
  labels <- lapply(fossils, function(level){
    tree <- level$tree
    descendents <- (extract.clade(tree, node = basal_node))$tip.label
    c(descendents, living_nodes)
  })
  labels
}

# recompute labels used later
labels <- lapply(fossil_matrices, extract.living)

# post-ordination ACE block
cat("Running post ord ace", replicate_id, "\n")

ord_fossil_tips <- lapply(fossil_matrices, lapply, function(x){
  mat <- x$matrix 
  dist <- char.diff(mat, method = "mord", by.col = FALSE)
  (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})

post_ord_ace <- Map(function(rate_matrix, rate_tree){
  Map(function(fossil_matrix, fossil_tree){
    multi.ace(fossil_matrix, fossil_tree, models = "ML", output = "multi.ace")
  }, rate_matrix, rate_tree)
}, ord_fossil_tips, fossil_trees)

saveRDS(post_ord_ace, write.path("anc", "post_ord_ace_%03d.rds"))

trait_normal <- list(fun = rnorm, param = list(mean = mean, sd = function(x) return(diff(range(x))/4)))
point_post_ord_ace <- lapply(post_ord_ace, lapply, multi.ace, output = "combined.matrix")

point_post_ord_ace_living <- Map(function(rate_anc, rate_labels) {
  Map(function(fossil_anc, label_anc) {
    fossil_anc[label_anc, , drop = FALSE]
  }, rate_anc, rate_labels)
}, point_post_ord_ace, labels)
names(point_post_ord_ace_living) <- names(point_post_ord_ace)

sample_post_ord_ace <- lapply(post_ord_ace, lapply, multi.ace, sample = 100, sample.fun = trait_normal, output = "combined.matrix")

sample_post_ord_ace_living <- Map(function(rate_anc, rate_labels) {
  Map(function(fossil_anc, label_anc) {
    lapply(fossil_anc, function(rep) rep[label_anc, , drop = FALSE])
  }, rate_anc, rate_labels)
}, sample_post_ord_ace, labels)
names(sample_post_ord_ace_living) <- names(sample_post_ord_ace)

cat("Finished post ord ace", replicate_id, "\n")

saveRDS(point_post_ord_ace_living, write.path("ord", "post_ord_point_%03d.rds"))
saveRDS(sample_post_ord_ace_living, write.path("ord", "post_ord_sample_%03d.rds"))

cat("Finished replicate", replicate_id, "\n")
# ...existing code...