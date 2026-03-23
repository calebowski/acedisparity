args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
tree_size <- args[2]
job_id <- args[3]
library(treats)

cat("Starting replicate", replicate_id, "\n")
set.seed(100 + replicate_id)
# job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")
base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/discrete/crown/", tree_size, "/")
write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}


# fossil_matrices <- readRDS(write.path("matrices", "fossil_matrices_%03d.rds"))
# relative_fossil_anc <- readRDS(write.path("anc", "pre_ord_rel_%03d.rds"))
# point_fossil_anc <- readRDS(write.path("anc", "pre_ord_point_%03d.rds"))
# matrices <- readRDS(write.path("matrices", "matrices_%03d.rds"))

# no_ace_matrix <- lapply(point_fossil_anc, lapply, function(matrix){matrix[!grepl("^n", rownames(matrix)),]}) # remove nodes

cat("Writing ord_no_ace...\n")

fossil_matrices <- readRDS(write.path("matrices", "fossil_matrices_%03d.rds"))
ord_no_ace <- lapply(fossil_matrices, lapply, function(rep){
mat <- rep$matrix
dist <- char.diff(mat, method = "mord", by.col = FALSE)
ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})
saveRDS(ord_no_ace, write.path("ord", "ord_no_ace_%03d.rds"))



cat("Writing ord_true...\n")
matrices <- readRDS(write.path("matrices", "matrices_%03d.rds"))
ord_true <- lapply(matrices, function(rep){
dist <- char.diff(rep, method = "mord", by.col = FALSE)
ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})
saveRDS(ord_true, write.path("ord", "ord_true_%03d.rds"))



cat("Writing ord_rel...\n")
relative_fossil_anc <- readRDS(write.path("anc", "pre_ord_rel_%03d.rds"))
ord_rel <- lapply(relative_fossil_anc, lapply, function(rep){
dist <- char.diff(rep, method = "mord", by.col = FALSE)
ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
}) ## ordinate full matrix, prune out certain datapoints later.
saveRDS(ord_rel, write.path("ord", "ord_rel_%03d.rds"))


cat("Writing ord_point...\n")
point_fossil_anc <- readRDS(write.path("anc", "pre_ord_point_%03d.rds"))
ord_point <- lapply(point_fossil_anc, lapply,  function(rep){
dist <- char.diff(rep, method = "mord", by.col = FALSE)
ord <- (cmdscale(dist, k = ncol(dist) - 2, add = TRUE))$points
})
saveRDS(ord_point, write.path("ord", "ord_point_%03d.rds"))

cat("Finished replicate", replicate_id, "\n")
