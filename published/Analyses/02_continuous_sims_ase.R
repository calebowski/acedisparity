library(treats)
library(MASS)
source("../Functions/fossil.pres.R")

dir.create("../Data/continuous", recursive = TRUE, showWarnings = FALSE)
dir.create("../Data/continuous", recursive = TRUE, showWarnings = FALSE)
dir.create("../Data/continuous/matrices", recursive = TRUE, showWarnings = FALSE)



## iterate over each replicate
tree_sizes <- c(50, 100, 150)
for(tree_size in tree_sizes) {
    for(i in 1:100) {
        base_path <- c("../Data/")
        replicate_id <- as.numeric(i)
        seed <- 100 + replicate_id
        set.seed(seed)

        tree <- read.tree(paste0("../Data/trees/", "tree_", tree_size, "t_", sprintf("%03d.tre", replicate_id)))

        ## Extract crown tree -----------------------------------------------------------------------------------
        ages <- tree.age(tree) # get tip ages
        extant <- ages$element[ages$ages == 0]
        living_tree <- keep.tip(tree, extant) ## root at mrca - alike to crown group analyses.
        crown_tree  <- extract.clade(tree, living_tree$node.label[1]) ## this is the crown tree
        tree_height <- max(node.depth.edgelength(crown_tree))

        ## Define trait model parameters -------------------------------------------------------------------------
        bm  <-  make.traits(process = BM.process, n = 20, 
                            process.args = list(Sigma = diag(0.25, 20)))

        bm_t <-  make.traits(process = BM.trend.process, n = 20, 
                            process.args = list(Sigma = diag(0.25, 20), trend = 0.3))

        ou_st <-  make.traits(process = OU.process, n = 20, 
                            process.args = list(alpha = (log(2) / (tree_height * 0.25)),
                                                Sigma = diag(0.25, 20)))

        traits <- list(bm = bm, bm_t = bm_t, ou_st = ou_st)

        ## Simulate traits across tree -------------------------------------------------------------------------
        matrices <- lapply(traits, function(x){map.traits(x, crown_tree)$data})
        saveRDS(matrices, sprintf("../Data/continuous/matrices/matrices_%03d.rds", replicate_id))


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

        saveRDS(fossil_matrices, sprintf("../Data/continuous/matrices/fossil_matrices_%03d.rds", replicate_id))


        # fossil_trees <- lapply(fossil_matrices, lapply, function(level){
        # tree <- level$tree
        # })

        ## Run ace ----------------------------------
        res_pre_ord_ace <- lapply(fossil_matrices, lapply,function(x){ 
            multi.ace(x$matrix, x$tree, models = "BM", output = "multi.ace", verbose = TRUE)
        })

        ## Get point estimates ----------------------------------

        point_ace <- lapply(res_pre_ord_ace, lapply, multi.ace, output = "combined.matrix") 

        ## Get distributions of ace (100 samples) ----------------------------------
        trait_normal  <-  list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4))) ## samples with normal distribution

        sample_ace <- lapply(res_pre_ord_ace, lapply, multi.ace, output = "combined.matrix", ml.collapse = list(type = "sample", sample = 100, sample.fun = trait_normal)) 

        saveRDS(fossil_anc, write.path("anc", "pre_ord_ace_%03d.rds"))
        saveRDS(point_anc, write.path("anc", "pre_ord_point_%03d.rds"))
        saveRDS(sample_anc, write.path("anc", "pre_ord_sample_%03d.rds"))






        






        # tree <- 
    }
}