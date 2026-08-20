replicate_id <- i
seed <- 100 + replicate_id
tree_size  <- 50
set.seed(seed)

tree <- read.tree(paste0("../Data/trees_alternative/", "tree_", sprintf("%st_%03d.tre", tree_size, replicate_id)))

# ## Extract crown tree -----------------------------------------------------------------------------------
# ages <- tree.age(tree) # get tip ages
# extant <- ages$element[ages$ages == 0]
# living_tree <- keep.tip(tree, extant) ## root at mrca - alike to crown group analyses.
# crown_tree  <- extract.clade(tree, living_tree$node.label[1]) ## this is the crown tree
tree_height <- max(node.depth.edgelength(tree))

## Define trait model parameters -------------------------------------------------------------------------
bm  <-  make.traits(process = BM.process, n = 20, process.args = list(Sigma = diag(0.25, 20)))

bm_t <-  make.traits(process = BM.trend.process, n = 20, process.args = list(Sigma = diag(0.25, 20), trend = 0.3))

ou_st <-  make.traits(process = OU.process, n = 20, process.args = list(alpha = (log(2) / (tree_height * 0.25)), Sigma = diag(0.25, 20)))

traits <- list(bm = bm, bm_t = bm_t, ou_st = ou_st)

## Simulate traits across tree -------------------------------------------------------------------------
matrices <- lapply(traits, function(x){map.traits(x, tree)$data})
saveRDS(matrices, sprintf("../Data/continuous/matrices/matrices_%st_%03d.rds", tree_size, replicate_id))


## Simulate fossil sampling ----------------------------------
living <- lapply(matrices, remove.fossil, trees = tree, type = "continuous")
fossilised_high <- lapply(matrices, fossil.pres.alt, trees = tree, preservation = 0.5, type = "continuous", seed = seed)
all_fossil <- lapply(matrices, fossil.pres.alt, trees = tree, preservation = 1.0, type = "continuous", seed = seed)
fossilised_med <- lapply(matrices, fossil.pres, trees = tree, preservation = 0.15, type = "continuous", seed = seed)
fossilised_low <- lapply(matrices, fossil.pres.alt, trees = tree, preservation = 0.05, type = "continuous", seed = seed)

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

saveRDS(fossil_matrices, sprintf("../Data/continuous/matrices/fossil_matrices_%st_%03d.rds", tree_size,replicate_id))