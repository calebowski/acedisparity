library(treats)
if(packageVersion("treats") != "1.1.7") {
    library(devtools)
    install_github("TGuillerme/treats", ref = "map.traits.events") ## requires this branch of treats
}

if(packageVersion("dispRity") != "1.9.8") {
    library(devtools)
    install_github("TGuillerme/dispRity")
}

library(MASS)
# library(here)
source("../Functions/fossil.pres.R")


dir.create("../Data/discrete", recursive = TRUE, showWarnings = FALSE)
dir.create("../Data/discrete/matrices", recursive = TRUE, showWarnings = FALSE)
dir.create("../Data/discrete/ord", recursive = TRUE, showWarnings = FALSE)
dir.create("../Data/discrete/ace", recursive = TRUE, showWarnings = FALSE)

run.sim.discrete.ace <- function(tree_size, replicate_id, samples = 100) {
    # replicate_id <- as.numeric(i)
    seed <- 100 + replicate_id
    set.seed(seed)
    tree <- read.tree(paste0("../Data/trees/", "tree_", sprintf("%st_%03d.tre", tree_size, replicate_id)))

    ## Extract crown tree -----------------------------------------------------------------------------------
    ages <- tree.age(tree) # get tip ages
    extant <- ages$element[ages$ages == 0]
    living_tree <- keep.tip(tree, extant) ## root at mrca - alike to crown group analyses.
    crown_tree  <- extract.clade(tree, living_tree$node.label[1]) ## this is the crown tree
    tree_height <- max(node.depth.edgelength(crown_tree))

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

    saveRDS(matrices, sprintf("../Data/discrete/matrices/matrices_%st_%03d.rds", tree_size, replicate_id))

    ## Simulate fossil sampling ----------------------------------
    living <- lapply(matrices, remove.fossil, trees = crown_tree, type = "discrete")
    fossilised_high <- lapply(matrices, fossil.pres, trees = crown_tree, preservation = 0.5, type = "discrete", seed = seed)
    all_fossil <- lapply(matrices, fossil.pres, trees = crown_tree, preservation = 1.0, type = "discrete", seed = seed)
    fossilised_med <- lapply(matrices, fossil.pres, trees = crown_tree, preservation = 0.15, type = "discrete", seed = seed)
    fossilised_low <- lapply(matrices, fossil.pres, trees = crown_tree, preservation = 0.05, type = "discrete", seed = seed)

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

    res_pre_ord_ace <- lapply(fossil_matrices, lapply, function(x) {
        multi.ace(x$matrix, x$tree, models = "ER", output = "multi.ace", verbose = TRUE)
    })

    ## Get point estimates ----------------------------------
    point_pre_ord_ace <- lapply(res_pre_ord_ace, lapply, multi.ace, ml.collapse = list(type = "majority", tie.breaker = TRUE), output = "combined.matrix")

    ## Get distributions of ace (100 samples) ----------------------------------
    sample_pre_ord_ace <- lapply(res_pre_ord_ace, lapply, multi.ace, output = "combined.matrix", ml.collapse = list(type = "sample", sample = samples))

    saveRDS(res_pre_ord_ace, sprintf("../Data/discrete/ace/pre_ord_ace_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(point_pre_ord_ace, sprintf("../Data/discrete/ace/pre_ord_point_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(sample_pre_ord_ace, sprintf("../Data/discrete/ace/pre_ord_sample_%st_%03d.rds", tree_size, replicate_id))

    ord_no_ace <- lapply(fossil_matrices, lapply, function(rep) {
    mat <- rep$matrix
    dist <- char.diff(mat, method = "mord", by.col = FALSE)
    ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
    })

    ord_true <- lapply(matrices, function(rep) {
    dist <- char.diff(rep, method = "mord", by.col = FALSE)
    ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
    })

    ord_pre_point <- lapply(point_pre_ord_ace, lapply,  function(rep) {
    dist <- char.diff(rep, method = "mord", by.col = FALSE)
    ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
    })

    ord_pre_sample <- lapply(sample_pre_ord_ace, lapply, lapply, function(rep) {
    dist <- char.diff(rep, method = "mord", by.col = FALSE)
    ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
    })

    ## Flatten tasks for post-ord ace ----------------------------------
    tasks_post_ord <- expand.grid(model = names(ord_no_ace), fossil_level = names(ord_no_ace[[1]]), stringsAsFactors = FALSE)

    res_post_ord_ace <- lapply(seq_len(nrow(tasks_post_ord)), function(i) {
        task <- tasks_post_ord[i, ]
        level_matrix <- ord_no_ace[[task$model]][[task$fossil_level]]
        level_tree <- fossil_trees[[task$model]][[task$fossil_level]]
        multi.ace(level_matrix, level_tree, models = "BM", output = "multi.ace", verbose = TRUE)
    })

    post_ord_ace  <- list()
    for(i in seq_along(res_post_ord_ace)){
        m <- tasks_post_ord$model[i]
        l <- tasks_post_ord$fossil_level[i]
        post_ord_ace[[m]][[l]] <- res_post_ord_ace[[i]]
    }

    ## Get point/distribution estimations for post-ord ace ----------------------------------
    point_post_ord_ace <- lapply(post_ord_ace, lapply, multi.ace, output = "combined.matrix")
    trait_normal  <-  list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4))) ## samples with normal distribution

    sample_post_ord_ace <- lapply(post_ord_ace, lapply, multi.ace, ml.collapse = list(type = "sample", sample = samples, sample.fun = trait_normal), output = "combined.matrix")

    saveRDS(ord_pre_sample, sprintf("../Data/discrete/ord/ord_pre_sample_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(ord_pre_point, sprintf("../Data/discrete/ord/ord_pre_point_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(ord_no_ace, sprintf("../Data/discrete/ord/ord_no_ace_%st_%03d.rds",tree_size, replicate_id))
    saveRDS(ord_true, sprintf("../Data/discrete/ord/ord_true_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(post_ord_ace, sprintf("../Data/discrete/ace/post_ord_ace_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(point_post_ord_ace, sprintf("../Data/discrete/ord/ord_post_point_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(sample_post_ord_ace, sprintf("../Data/discrete/ord/ord_post_sample_%st_%03d.rds", tree_size, replicate_id))

}

# Single replicate for testing reproducibility --------------------------------------
run.sim.discrete.ace(tree_size = 50, replicate_id = 1, samples = 10) ## lower sample count for computational speed

# Run all replicates
tree_sizes <- c(50, 100, 150)
for(tree_size in tree_sizes) {
    for(i in 1:100) {
        run.sim.discrete.ace(tree_size, i)
    }
}