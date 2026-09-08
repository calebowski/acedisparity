args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
tree_size <- args[2]
library(parallel)
library(treats)
library(readr)
library(phytools)

cat("Starting replicate", replicate_id, "\n")

set.seed(100 + replicate_id)

base_path <- "/mnt/parscratch/users/bip24cns/acedisparity/revisions/discrete_wagner2/"
job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")


write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}


tree <- read.tree(paste0("/mnt/parscratch/users/bip24cns/acedisparity/trees/overallDisparity/tree_", tree_size, sprintf("_%03d.tre", replicate_id)))

# tree <- set.root.time(tree)

ages <- tree.age(tree) # get tip ages
extant <- ages$element[ages$ages == 0]
living_tree <- keep.tip(tree, extant) ## root at mrca - alike to crown group analyses.

crown_tree  <- extract.clade(tree, living_tree$node.label[1]) ## this is the crown tree


simulate.punc.eq.traits <- function(tree,n.binary = 85,n.multi = 15, change.prob, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  tree <- ape::reorder.phylo(tree, order = "cladewise")

  n_characters <- n.binary + n.multi
  n_states <- c(
    rep(2, n.binary),
    rep(3, n.multi)
  )

  n_nodes <- ape::Ntip(tree) + tree$Nnode
  states <- matrix(
    0,
    nrow = n_nodes,
    ncol = n_characters
  )

  # Root states are all zero. Simulate each descendant from its parent.
  for (edge_index in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[edge_index, 1]
    child <- tree$edge[edge_index, 2]

    states[child, ] <- states[parent, ]

    n_changes <- 0
    while (n_changes == 0) {
      n_changes <- rbinom(
        n = 1,
        size = n_characters,
        prob = change.prob
      )
    }

    changed_characters <- sample(
      seq_len(n_characters),
      size = n_changes,
      replace = FALSE
    )

    for (character in changed_characters) {
      ancestral_state <- states[parent, character]

      possible_states <- setdiff(
        seq_len(n_states[character]) - 1,
        ancestral_state
      )

      states[child, character] <- sample(
        possible_states,
        size = 1
      )
    }
  }

  row_labels <- c(tree$tip.label, tree$node.label)

  rownames(states) <- row_labels
  return(states)
}

cat("Mapping traits...\n")

trait_change_probabilities <- c( ## assign same transition rates as before.
  slow = 0.01,
  med = 0.10,
  fast = 0.90
)

matrices <- lapply(trait_change_probabilities,function(changeprob) {
  simulate.punc.eq.traits(tree = crown_tree, n.binary = 85, n.multi = 15, change.prob = changeprob, seed = replicate_id + 100)
})

saveRDS(matrices, write.path("matrices", "matrices_%03d.rds"))

cat("Trait matrices saved for replicate", replicate_id, "\n")

source("/users/bip24cns/acedisparity/discrete/scripts/fossil.pres.R")
set_seed <- 100 + replicate_id
living <- lapply(matrices, remove.fossil, trees = crown_tree, type = "discrete")
fossilised_all <- lapply(matrices, fossil.pres.alt, trees = crown_tree, preservation = 0.5, type = "discrete", seed = set_seed)
fossilised_high <- lapply(matrices, fossil.pres.alt, trees = crown_tree, preservation = 0.5, type = "discrete", seed = set_seed)
fossilised_med <- lapply(matrices, fossil.pres.alt, trees = crown_tree, preservation = 0.15, type = "discrete", seed = set_seed)
fossilised_low <- lapply(matrices, fossil.pres.alt, trees = crown_tree, preservation = 0.05, type = "discrete", seed = set_seed)


fossil_matrices <- lapply(names(matrices), function(level) {
    list(
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

res_pre_ord_ace <- mclapply(seq_len(nrow(tasks)), function(i){ ## loop over each model combination by row
    task <- tasks[i,]
    level <- fossil_matrices[[task$rate]][[task$fossil_level]]
    tryCatch({
    multi.ace(level$matrix, level$tree, models = "ER", output = "multi.ace")
  }, error = function(e) {
    cat("ERROR:", task$fossil_level, e$message, "\n") ## error handeling
    NULL
  })
}, mc.cores = 15)

pre_ord_ace <- list()
for(i in seq_along(res_pre_ord_ace)) {
  r <- tasks$rate[i]  
  l <- tasks$fossil_level[i]
  pre_ord_ace[[r]][[l]] <- res_pre_ord_ace[[i]]
}

# fossil_anc <- lapply(fossil_matrices, lapply, anc.states)

cat("Ancestral states estimated\n")

sample_fossil_anc <- lapply(pre_ord_ace, lapply,  multi.ace, ml.collapse = list(type = "sample", sample = 100))

point_fossil_anc <- lapply(pre_ord_ace, lapply, multi.ace, ml.collapse = list(type = "majority", tie.breaker = TRUE), output = "combined.matrix", verbose = TRUE)


saveRDS(pre_ord_ace, write.path("anc", "pre_ord_anc_%03d.rds"))
saveRDS(sample_fossil_anc, write.path("anc", "pre_ord_sample_%03d.rds"))
saveRDS(point_fossil_anc, write.path("anc", "pre_ord_point_%03d.rds"))

########################################################################################################################

cat("Finished replicate", replicate_id, "\n")
