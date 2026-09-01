args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
tree_size <- args[2]
library(parallel)
library(treats)
library(readr)
library(phytools)

cat("Starting replicate", replicate_id, "\n")

set.seed(100 + replicate_id)

base_path <- "/mnt/parscratch/users/bip24cns/acedisparity/revisions/discrete/"
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


cat("Mapping traits...\n")


slow_binary_transitions <- matrix(c(
    0.99, 0.01,
    0.01, 0.99
), nrow = 2, byrow = TRUE)

slow_multi_transitions <- matrix(c(
    0.99, 0.005, 0.005,
    0.005, 0.99, 0.005,
    0.005, 0.005, 0.99
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
  mapped_binary <- (map.traits(traits$binary, crown_tree))$data
  mapped_multi <- (map.traits(traits$multi, crown_tree))$data
  cbind(mapped_binary, mapped_multi)
})

saveRDS(matrices, write.path("matrices", "matrices_%03d.rds"))

cat("Trait matrices saved for replicate", replicate_id, "\n")

source("/users/bip24cns/acedisparity/discrete/scripts/fossil.pres.R")
set_seed <- 100 + replicate_id
living <- lapply(matrices, remove.fossil, trees = crown_tree, type = "discrete")
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
