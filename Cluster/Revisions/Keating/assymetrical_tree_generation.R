library(treats)
library(treestats)
 ## three tree sizes
bd_params <- make.bd.params(speciation = 1, extinction = 0.90)

## Biased selection
bias.select <- function(lineage) {
    ## Sample one lineage among the existing lineages
    ## "lineage" is an internal treats object details in manual.
    ## The sample proportion is a decreasing log exponential distribution
    ## i.e. the last lineage has always more changes to be selected
    ## You can modify the rate parameter internally. Bigger = more ladderised
    probs <- rev(dexp(seq(1, lineage$n, by = 1), rate = 1.39))
    return(sample(1:lineage$n, 1, prob = probs))
}

n <- 10
lineages <- seq_len(n)
weights <- rev(dexp(lineages, rate = 1.39))
probabilities <- weights / sum(weights)

# data.frame(
#     lineage = lineages,
#     weight = weights,
#     probability = probabilities
# ) ## check the probabilities of each lineage being sampled is 75% at this

fix.zero.branches <- function(tree, min_length = 1e-8) {
  tree$edge.length[tree$edge.length <= 0] <- min_length
  return(tree)
}

## A modifier for selection
biased_modifier <- make.modifiers(selection = bias.select)

living_tips <- 50
n_replicates <- 100
metadata_list <- list()
for(i in 1:n_replicates) {
  replicate_id <- as.numeric(i)
  seed <- 100 + replicate_id
  set.seed(seed)

    tree_size <- as.character(living_tips)
    max_total_size <- living_tips * 10
    
    # Keep generating trees until we get one within size limit
    repeat {
      # set.seed(seed_val)
      tree <- treats(stop.rule = list(max.living = living_tips), bd.params = bd_params, null.error = 1e6, modifiers = biased_modifier)
      tree <- drop.singles(tree)
      tree <- fix.zero.branches(tree)
      
      total_tips <- length(tree$tip.label)
      
      # Check if tree is within acceptable size
      if (total_tips <= max_total_size) {
        break  # Accept this tree
      } else {
        cat("Tree too large (", total_tips, " > ", max_total_size, "), regenerating...\n")
        # seed_val <- seed_val + 100  # Change seed for next attempt
      }
    }
    tree_file <- sprintf("../../Data/revisions/keating/trees/assymetrical_tree_%st_%03d.tre", tree_size, replicate_id)
    write.tree(tree, tree_file)

    est <- crude.bd.est(tree, "estimate")
    colless <- colless_corr(tree)

    metadata_list[[i]] <- data.frame(
      replicate_id = replicate_id,
      total_size = total_tips,
      seed = seed,
      speciation = est$call$speciation, 
      extinction = est$call$extinction,
      colless = colless,
      stringsAsFactors = FALSE
    )
}
metadata_df <- do.call(rbind, metadata_list)
hist(metadata_df$colless)
metadata_file <- "../../Data/revisions/keating/trees/assymetrical_metadata_trees.csv"
write.csv(metadata_df, metadata_file, row.names = FALSE)
