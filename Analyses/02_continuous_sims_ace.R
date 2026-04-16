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
source("../Functions/fossil.pres.R")
source("../Functions/utility.R")


dir.create("../Data/continuous", recursive = TRUE, showWarnings = FALSE)
dir.create("../Data/continuous/matrices", recursive = TRUE, showWarnings = FALSE)
dir.create("../Data/continuous/ord", recursive = TRUE, showWarnings = FALSE)
dir.create("../Data/continuous/ace", recursive = TRUE, showWarnings = FALSE)



## iterate over each replicate
# tree_sizes <- c(50, 100, 150)
# for(tree_size in tree_sizes) {
#     for(i in 1:100) {

run.sim.cont.ace <- function(tree_size, replicate_id, samples = 100) {
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

    ## Define trait model parameters -------------------------------------------------------------------------
    bm  <-  make.traits(process = BM.process, n = 20, process.args = list(Sigma = diag(0.25, 20)))

    bm_t <-  make.traits(process = BM.trend.process, n = 20, process.args = list(Sigma = diag(0.25, 20), trend = 0.3))

    ou_st <-  make.traits(process = OU.process, n = 20, process.args = list(alpha = (log(2) / (tree_height * 0.25)), Sigma = diag(0.25, 20)))

    traits <- list(bm = bm, bm_t = bm_t, ou_st = ou_st)

    ## Simulate traits across tree -------------------------------------------------------------------------
    matrices <- lapply(traits, function(x){map.traits(x, crown_tree)$data})
    saveRDS(matrices, sprintf("../Data/continuous/matrices/matrices_%st_%03d.rds", tree_size, replicate_id))


    ## Simulate fossil sampling ----------------------------------
    living <- lapply(matrices, remove.fossil, trees = crown_tree, type = "continuous")
    fossilised_high <- lapply(matrices, fossil.pres, trees = crown_tree, preservation = 0.5, type = "continuous", seed = seed)
    all_fossil <- lapply(matrices, fossil.pres, trees = crown_tree, preservation = 1.0, type = "continuous", seed = seed)
    fossilised_med <- lapply(matrices, fossil.pres, trees = crown_tree, preservation = 0.15, type = "continuous", seed = seed)
    fossilised_low <- lapply(matrices, fossil.pres, trees = crown_tree, preservation = 0.05, type = "continuous", seed = seed)

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

    fossil_trees <- lapply(fossil_matrices, lapply, function(level) {
    tree <- level$tree
    })

    saveRDS(fossil_trees, sprintf("../Data/trees/fossil_trees_%st_%03d.rds", tree_size, replicate_id))

    ## Run ace ----------------------------------
    res_pre_ord_ace <- lapply(fossil_matrices, lapply,function(x) { 
        multi.ace(x$matrix, x$tree, models = "BM", output = "multi.ace", verbose = TRUE)
    })

    ## Get point estimates ----------------------------------

    point_pre_ord_ace <- lapply(res_pre_ord_ace, lapply, multi.ace, output = "combined.matrix")

    ## Get distributions of ace (100 samples) ----------------------------------
    trait_normal  <-  list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4))) ## samples with normal distribution

    sample_pre_ord_ace <- lapply(res_pre_ord_ace, lapply, multi.ace, output = "combined.matrix", ml.collapse = list(type = "sample", sample = samples, sample.fun = trait_normal)) 

    saveRDS(res_pre_ord_ace, sprintf("../Data/continuous/ace/pre_ord_ace_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(point_pre_ord_ace, sprintf("../Data/continuous/ace/pre_ord_point_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(sample_pre_ord_ace, sprintf("../Data/continuous/ace/pre_ord_sample_%st_%03d.rds", tree_size, replicate_id))

    ## Ordinations ----------------------------------
    ord_pre_sample <- lapply(sample_pre_ord_ace, lapply, lapply, function(mat){
        prcomp(mat, scale = FALSE, center = TRUE)$x
    })

    ord_pre_point <- lapply(point_pre_ord_ace, lapply, function(mat){
        prcomp(mat, scale = FALSE, center = TRUE)$x
    })

    ord_no_ace <- lapply(fossil_matrices, lapply, function(mat){
        prcomp(mat$matrix, scale = FALSE, center = TRUE)$x
    })

    ord_true <- lapply(matrices, function(mat) {
    prcomp(mat, scale = FALSE, center = TRUE)$x
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

    sample_post_ord_ace <- lapply(post_ord_ace, lapply, multi.ace, ml.collapse = list(type = "sample", sample = samples, sample.fun = trait_normal), output = "combined.matrix")

    saveRDS(ord_pre_sample, sprintf("../Data/continuous/ord/ord_pre_sample_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(ord_pre_point, sprintf("../Data/continuous/ord/ord_pre_point_%st_%03d.rds",tree_size, replicate_id))
    saveRDS(ord_no_ace, sprintf("../Data/continuous/ord/ord_no_ace_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(ord_true, sprintf("../Data/continuous/ord/ord_true_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(post_ord_ace, sprintf("../Data/continuous/ace/post_ord_ace_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(point_post_ord_ace, sprintf("../Data/continuous/ord/ord_post_point_%st_%03d.rds", tree_size, replicate_id))
    saveRDS(sample_post_ord_ace, sprintf("../Data/continuous/ord/ord_post_sample_%st_%03d.rds", tree_size, replicate_id))

}


# Single replicate for testing reproducibility --------------------------------------
run.sim.cont.ace(tree_size = 50, replicate_id = 1, samples = 10) ## lower sampling for computational speed

# Run all replicates
tree_sizes <- c(50, 100, 150)
for(tree_size in tree_sizes) {
    for(i in 1:100) {
        run.sim.cont.ace(tree_size, i, samples = 100)
    }
}