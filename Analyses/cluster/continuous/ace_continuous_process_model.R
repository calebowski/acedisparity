args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
model_name <- args[2] 
tree_size <- args[3]
job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")

library(treats)
library(parallel)
library(MASS)

cat("Starting replicate", replicate_id, "\n")


base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/continuous/", tree_size, "/")
source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")

write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}

# Test 1: Read tree
# cat("\n=== STEP 1: Reading tree ===\n")
# tree_file <- write.path("trees", "tree_%03d.tre")
# cat("Tree file:", tree_file, "\n")

# tryCatch({
#   tree <- read.tree(tree_file)
#   cat("SUCCESS: Read tree with", length(tree$tip.label), "tips\n")
# }, error = function(e) {
#   cat("ERROR reading tree:", e$message, "\n")
#   quit(status = 1)
# })

cat("\n=== STEP 1: Generating tree ===\n")
set.seed(100 + replicate_id)  # Set seed for reproducible tree generation

tree_files <- list.files(paste0("/mnt/parscratch/users/bip24cns/acedisparity/discrete/", tree_size, "/trees"))
filename <- tree_files[1]
tree_job_id <- sub("_tree_.*$", "", filename)
tree <- read.tree(paste0("/mnt/parscratch/users/bip24cns/acedisparity/discrete/100t/trees/",tree_job_id, sprintf("_tree_%03d.tre", replicate_id)))

tree <- drop.singles(tree)
tree <- fix.zero.branches(tree)
tree <- set.root.time(tree)


tree_height <- max(node.depth.edgelength(tree))



cat("Tree generated for replicate", replicate_id, "\n")


traits <- switch(model_name,
  "bm" = make.traits(process = BM.process, n = 20, 
                     process.args = list(Sigma = diag(0.25, 20))),
  
  "bm_t" = make.traits(process = BM.trend.process, n = 20, 
                      process.args = list(Sigma = diag(0.25, 20), trend = 0.3)),
  
  "ou_w" = make.traits(process = OU.process, n = 20, 
                      process.args = list(alpha = (log(2) / (tree_height * 0.75)), 
                                        Sigma = diag(0.25, 20))),
  
  "ou_st" = make.traits(process = OU.process, n = 20, 
                       process.args = list(alpha = (log(2) / (tree_height * 0.25)), 
                                         Sigma = diag(0.25, 20))),
  
  "ou_sh" = make.traits(process = OU.process, n = 20, 
                       process.args = list(optimum = 2, 
                                         alpha = (log(2) / (tree_height * 0.75)), 
                                         Sigma = diag(0.25, 20))),
  
  stop("Unknown model: ", model_name)
)

cat("Starting model", model_name, "for replicate", replicate_id, "at", Sys.time(), "\n")

mat <- map.traits(traits, tree)$data

cat("Traits simulated for", model_name, "\n")

# matrices <- list(bm = bm_matrices, bm_t = bm_matrices_trend, ou_w = ou_w_matrices, ou_st = ou_s_matrices, ou_sh = ou_shift_matrices)

saveRDS(mat, write.path("matrices", paste0(model_name, "_matrix_%03d.rds")))

cat("Trait matrices saved for replicate", replicate_id, model_name, "\n")

source("/users/bip24cns/acedisparity/discrete/scripts/fossil.pres.R")

living <- remove.fossil(mat, trees = tree, type = "continuous")
fossilised_high <- fossil.pres(mat,  trees = tree, preservation = 0.5, type = "continuous")
all_fossil <- fossil.pres(mat,  trees = tree, preservation = 1.0, type = "continuous")
fossilised_med <- fossil.pres(mat,  trees = tree, preservation = 0.15, type = "continuous")
fossilised_low <- fossil.pres(mat,  trees = tree, preservation = 0.05, type = "continuous")


fossil_matrices <- list(
      all = all_fossil,
      fossil_high = fossilised_high,
      fossil_med = fossilised_med,
      fossil_low = fossilised_low,
      living = living
)

# names(fossil_matrices) <- names(mat)

fossil_trees <- lapply(fossil_matrices, function(level){
  tree <- level$tree
})
cat("Fossil matrices completed...\n")

saveRDS(fossil_matrices, write.path("matrices", paste0(model_name, "_fossil_matrices_%03d.rds")))
saveRDS(fossil_trees, write.path("trees", paste0(model_name, "_fossil_trees_%03d.rds")))


tasks  <- expand.grid(fossil_level = names(fossil_matrices), stringsAsFactors = FALSE) 

cl <- makeCluster(5) ## make 25 core cluster (one for each task)
clusterEvalQ(cl, library(treats))
clusterExport(cl, c("fossil_matrices", "tasks"))


res <- parLapply(cl, seq_len(nrow(tasks)), function(i){ ## loop over each model combination by row
    task <- tasks[i,]
    level <- fossil_matrices[[task]]
    tryCatch({
    multi.ace(level$matrix, level$tree, models = "BM", output = "multi.ace")
  }, error = function(e) {
    cat("ERROR:", task$fossil_level, e$message, "\n") ## error handeling
    NULL
  })
})

stopCluster(cl)

fossil_anc <- list()
for(i in seq_along(res)) {
  l <- tasks$fossil_level[i]
  fossil_anc[[l]] <- res[[i]]
}

cat("Starting point ancestral state estimation...\n")
point_anc <- lapply(fossil_anc, multi.ace, output = "combined.matrix")  ## takes median of confidence interval

trait_normal  <-  list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4))) ## samples with normal distribution

cat("Starting distribution ancestral state estimation...\n")
sample_anc <- lapply(fossil_anc, multi.ace, output = "combined.matrix", sample = 100, sample.fun = trait_normal) ## 100 matrices replicated, 


cat("Ancestral states estimated...\n")

saveRDS(fossil_anc, write.path("anc", paste0(model_name, "_fossil_anc_%03d.rds")))
saveRDS(point_anc, write.path("anc", paste0(model_name, "_point_anc_%03d.rds")))
saveRDS(sample_anc, write.path("anc", paste0(model_name, "_sample_anc_%03d.rds")))
cat("Ancestral states saved...\n")

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

labels <- extract.living(fossil_matrices)

point_living <- Map(function(fossil_anc, label_anc) {
      fossil_anc[label_anc, , drop = FALSE]
}, point_anc, labels)
names(point_living) <- names(point_anc)


sample_living <- Map(function(fossil_anc, label_anc) {
      lapply(fossil_anc, function(rep) {
        rep[label_anc, , drop = FALSE]
      })
}, sample_anc, labels)

names(sample_living) <- names(sample_anc)


no_ace_living <- lapply(point_living, function(matrix){
  no_node <- matrix[!grepl("^n", rownames(matrix)),]
})

true_living <- mat[labels$all, , drop = FALSE]


cat("Living matrices completed...\n")


cat("=== STARTING ord_sample ===\n")
  ord_sample <- lapply(sample_living, lapply, function(mat){
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })


cat("=== STARTING ord_point ===\n")
  ord_point <- lapply(point_living, function(mat){
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })


cat("=== STARTING ord_no_ace ===\n")
  ord_no_ace <- lapply(no_ace_living, function(mat){
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })


cat("=== STARTING ord_true ===\n")
  ord_true <- prcomp(true_living, scale = FALSE, center = TRUE)$x



cat("Ordinations completed...\n")


cat("=== STARTING ordination of fossil tips for post ordiation ace ===\n")
ord_fossil_tips <- lapply(fossil_matrices, function(x) {
    mat <- x$matrix
    prcomp(mat, scale = FALSE, center = TRUE)$x
  })

cat("=== STARTING ancestral statte estimation of post ordination ===\n")

tasks_post_ord <- expand.grid(fossil_level = names(ord_fossil_tips), stringsAsFactors = FALSE)

cl <- makeCluster(5)
clusterEvalQ(cl, library(treats))
clusterExport(cl, c("ord_fossil_tips", "fossil_trees", "tasks_post_ord"))

res_post_ord <- parLapply(cl, seq_len(nrow(tasks_post_ord)), function(i){
    task <- tasks_post_ord[i,]
    ord_matrix <- ord_fossil_tips[[task]]
    fossil_tree <- fossil_trees[[task]]
    tryCatch({
        multi.ace(ord_matrix, fossil_tree, models = "BM", output = "multi.ace")
    }, error = function(e) {
        cat("ERROR in level:", task, e$message, "\n")
    })
})
stopCluster(cl)

post_ord_ace  <- list()
for(i in seq_along(res_post_ord)){
    l <- tasks_post_ord$fossil_level[i]
    post_ord_ace[[l]] <- res_post_ord[[i]]
}

cat("=== STARTING point estimation post ord ace ===\n")

point_post_ord_ace <- lapply(post_ord_ace,  multi.ace, output = "combined.matrix")

point_post_ord_ace_living <- Map(function(fossil_anc, label_anc) {
      fossil_anc[label_anc, , drop = FALSE]
}, point_post_ord_ace, labels)

names(point_post_ord_ace_living) <- names(point_post_ord_ace)

trait_normal <-  list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))

sample_post_ord_ace <- lapply(post_ord_ace,  multi.ace, sample = 100, sample.fun = trait_normal, output = "combined.matrix")

sample_post_ord_ace_living <- Map(function(fossil_anc, label_anc) {
      lapply(fossil_anc, function(rep) {
        rep[label_anc, , drop = FALSE]
      })
}, sample_post_ord_ace, labels)
cat("Post ordination ace completed...\n")

saveRDS(ord_sample, write.path("ord/temp", paste0(model_name, "_ord_sample_%03d.rds")))
saveRDS(ord_point, write.path("ord/temp", paste0(model_name, "_ord_point_%03d.rds")))
saveRDS(ord_no_ace, write.path("ord/temp", paste0(model_name, "_ord_no_ace_%03d.rds")))
saveRDS(ord_true, write.path("ord/temp", paste0(model_name, "_ord_true_%03d.rds")))
saveRDS(point_post_ord_ace_living, write.path("ord/temp", paste0(model_name, "_post_ord_point_%03d.rds")))
saveRDS(sample_post_ord_ace_living, write.path("ord/temp", paste0(model_name, "_post_ord_sample_%03d.rds")))


cat("Finished model", model_name, "for replicate", replicate_id, "at", Sys.time(), "\n")
