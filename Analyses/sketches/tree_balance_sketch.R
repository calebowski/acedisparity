library(treats)
library(treestats)
source("../Functions/fossil.pres.R")
source("../Functions/utility.R")
tree_size <- 100
replicate_id <- 1
tree <- read.tree(paste0("../Data/trees/", "tree_", sprintf("%st_%03d.tre", tree_size, replicate_id)))

ages <- tree.age(tree) # get tip ages
extant <- ages$element[ages$ages == 0]
living_tree <- keep.tip(tree, extant) ## root at mrca - alike to crown group analyses.
crown_tree  <- extract.clade(tree, living_tree$node.label[1]) ## this is the crown tree
tree_height <- max(node.depth.edgelength(crown_tree))

## Define trait model parameters -------------------------------------------------------------------------
bm  <-  make.traits(process = BM.process, n = 20, process.args = list(Sigma = diag(0.25, 20)))

bm_t <-  make.traits(process = BM.trend.process, n = 20, process.args = list(Sigma = diag(0.25, 20), trend = 0.3))

ou_st <-  make.traits(process = OU.process, n = 20, process.args = list(alpha = (log(2) / (tree_height * 0.25)), Sigma = diag(0.25, 20)))

traits <- list(bm = bm, bm_t = bm_t, ou_st = ou_st)

## Simulate traits across tree -------------------------------------------------------------------------
matrices <- lapply(traits, function(x){map.traits(x, crown_tree)$data})


## Simulate fossil sampling ----------------------------------
living <- lapply(matrices, remove.fossil, trees = crown_tree, type = "continuous")
fossilised_high <- lapply(matrices, non.random.fossil.pres, trees = crown_tree, preservation = c(0.05, 0.5), type = "continuous", seed = seed)
all_fossil <- lapply(matrices, non.random.fossil.pres, trees = crown_tree, preservation = c(0.05, 1), type = "continuous", seed = seed)
fossilised_med <- lapply(matrices, non.random.fossil.pres, trees = crown_tree, preservation = c(0.05, 0.15), type = "continuous", seed = seed)
fossilised_low <- lapply(matrices, non.random.fossil.pres, trees = crown_tree, preservation = c(0.05, 0.05), type = "continuous", seed = seed)

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



trees <- lapply(fossil_matrices, lapply, function(x) x$tree)[[1]]


# library(ape)

# n_tips <- 10

# ladder_tree <- stree(n = n_tips, type = "right")


# ladder_tree <- compute.brlen(ladder_tree, method = "Grafen")

# plot(ladder_tree, main = "Completely Asymmetrical (Pectinate) Tree")





########################################################################################################################

bd_params <- make.bd.params(speciation = 1, extinction = 0.90)
living_tips <- 50

split.exponential <- function(n, r = 0.5) {
  weights <- r^(0:(n - 1))
  
  splits <- weights / sum(weights)
  
  return(splits)
}

select.last <- function(lineage) {
    return(as.integer(sample(1:lineage$n, 1, prob = split.exponential(lineage$n, r = 0.2)))) ## last one should be highest, not first
}


select.last <- function(lineage) {
    return(as.integer(sample(1:lineage$n, 1, prob = split.exponential(lineage$n, r = 0.2))))
}

ladderised_modifier <- make.modifiers(selection = select.last)

set.seed(0)
tree_size <- as.character(living_tips)
max_total_size <- living_tips * 10

# Keep generating trees until we get one within size limit

# set.seed(seed_val)
tree <- treats(stop.rule = list(max.living = living_tips), bd.params = bd_params, null.error = 1e6, modifiers = ladderised_modifier)
tree <- drop.singles(tree)
tree <- fix.zero.branches(tree)

total_tips <- length(tree$tip.label)



split.exponential(100, r = 0.5)






distance.modify <- function(x, trait.values, lineage) {
     ## Distance to the parent's trait
     parent_trait_val <- parent.traits(trait.values, lineage)[1]
     mean_trait_val <- mean(trait.values[, 1])
     distance <- abs(parent_trait_val - mean_trait_val)
     ## Scales x with the distance
     return(x + x * distance)
}

distance.speciation <- make.modifiers(speciation = speciation,
                                      modify = distance.modify)



#     # Check if tree is within acceptable size
#     if (total_tips <= max_total_size) {
#     break  # Accept this tree
#     } else {
#     cat("Tree too large (", total_tips, " > ", max_total_size, "), regenerating...\n")
#     # seed_val <- seed_val + 100  # Change seed for next attempt
#     }
# }


# ## Our function that only select taxa with positive trait values
# select.positive <- function(trait.values, lineage) {

#     ## Selecting the taxa names with positive values for the first trait
#     positives <- as.integer(rownames(trait.values)[which(trait.values[, 1] >= 0)])

#     ## Combine the descendants of the current lineages (lineage$parents)
#     ## with the species that have speciated (seq_along(lineages$split))
#     ## to have a table of pairs of parents/splits
#     parents_split_table <- cbind(lineage$parents, seq_along(lineage$split))
#     ## Select the current taxa that descend from a node with a positive value
#     positive_living <- parents_split_table[which(lineage$parents %in% positives), 2]

#     ## Select one tip randomly in the ones with descendants with positive values
#     return(sample(which(lineage$livings %in% positive_living), 1))
# }

## Creating the modifier
positive_skew <- make.modifiers(selection = select.positive)

## Creating a (default) trait object
BM_trait <- make.traits()

## Simulate a tree and trait with no modifier
set.seed(1)
default_treats <- treats(bd.params = bd_params,
                         stop.rule = stop_rule,
                         traits    = BM_trait)

## Simulate a tree and trait with the modifier
set.seed(1)
skewed_trait_treats <- treats(bd.params = bd_params,
                              stop.rule = stop_rule,
                              traits    = BM_trait,
                              modifiers = positive_skew)

## Plotting the differences in trees and traits
par(mfrow = c(1, 2))
plot(default_treats, main = "Default trait and tree")
plot(skewed_trait_treats, main = "Skewed trait and tree")