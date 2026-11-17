args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
tree_size <- args[2]
job_id <- args[3]

library(dispRity)


base_path <- "/mnt/parscratch/users/bip24cns/acedisparity/revisions/discrete_wagner/"
write.path <- function(subfolder, filename) {
  paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id))
}


matrices_true <- readRDS(write.path("matrices", "matrices_%03d.rds"))



fossil_matrices <- readRDS(write.path("matrices", "fossil_matrices_%03d.rds"))



point_anc <- readRDS(write.path("anc", "pre_ord_point_%03d.rds"))



distances_true <- lapply(matrices_true,  char.diff, method = "hamming", by.col = FALSE)

distances_sampled_only <- lapply(fossil_matrices, lapply,function(x) char.diff(x$matrix, method = "hamming", by.col = FALSE))

distances_ace <- lapply(point_anc, lapply, char.diff, method = "hamming", by.col = FALSE)

mean.pairwise <- function(distance_matrix) {
  mean(distance_matrix[upper.tri(distance_matrix)], na.rm = TRUE)
}

mean_pairwise_true <- lapply(distances_true,mean.pairwise)

mean_pairwise_sampled <- lapply(distances_sampled_only, lapply,  mean.pairwise)

mean_pairwise_ace <- lapply(distances_ace, lapply,  mean.pairwise)


ace_errors <- Map(function(rate_true, rate_ace){
        lapply(rate_ace, function(fossil_level){
            fossil_level - rate_true
        })
}, mean_pairwise_true, mean_pairwise_ace)


sampled_errors <- Map(function(rate_true, rate_sampled){
    lapply(rate_sampled, function(fossil_level){
        fossil_level - rate_true
    })
}, mean_pairwise_true, mean_pairwise_sampled)

saveRDS(ace_errors, write.path("dist_disparity", "ace_dist_disparity_errors_%03d.rds"))
saveRDS(sampled_errors, write.path("dist_disparity", "sampled_only_dist_disparity_errors_%03d.rds"))

