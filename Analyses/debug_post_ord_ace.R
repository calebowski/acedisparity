library(dispRity)

reps  <- c(1, 47, 86)

job_id <- "8558401"
ord_true <- list()
post_ord_ace_sample <- list()
post_ord_ace_point <- list()
fossil_matrices <- list()
ord_no_ace <- list()
post_ord_ace_sample_all <- list()


for(i in seq_along(reps)) {
    rep <- reps[i]
    fossil_matrices[[i]] <- readRDS(paste0("../Data/cluster/discrete_crown/50t/matrices/", job_id, sprintf("_fossil_matrices_%03d.rds", rep)))
    ord_true[[i]] <- readRDS(paste0("../Data/cluster/discrete_crown/50t/ord/", job_id, sprintf("_ord_true_%03d.rds", rep)))
    post_ord_ace_sample[[i]] <- readRDS(paste0("../Data/cluster/discrete_crown/50t/checkpoints_ace/", sprintf("rep%03d", rep), "/fast/living/", job_id, "_sample_post_ord_ace.rds"))
    post_ord_ace_sample_all[[i]] <- readRDS(paste0("../Data/cluster/discrete_crown/50t/checkpoints_ace/", sprintf("rep%03d", rep), "/fast/all/", job_id, "_sample_post_ord_ace.rds"))
    post_ord_ace_point[[i]] <- readRDS(paste0("../Data/cluster/discrete_crown/50t/checkpoints_ace/", sprintf("rep%03d", rep), "/fast/living/", job_id, "_point_post_ord_ace.rds"))
    ord_no_ace[[i]] <- readRDS(paste0("../Data/cluster/discrete_crown/50t/ord/", job_id, sprintf("_ord_no_ace_%03d.rds", rep)))
}


scale_01 <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}



post_ord_ace_sample_rm_axes <- lapply(post_ord_ace_sample, lapply, function(x){
    apply(x[, 1:48], 2, scale_01)
})

ord_true_rm_axes <- lapply(ord_true, lapply, function(x){
    apply(x[, 1:48], 2, scale_01)
})

post_ord_ace_sample_rm_axes_all <- lapply(post_ord_ace_sample_all, lapply, function(x){
    apply(x[, 1:48], 2, scale_01)
})


disp_true <- lapply(ord_true_rm_axes, lapply, function(x){ get.disparity(dispRity(x, metric = c(sum, variances)))})


disp_post_ord_ace_sample <- lapply(post_ord_ace_sample_rm_axes, lapply,  function(x){ get.disparity(dispRity(x, metric = c(sum, variances)))})

disp_post_ord_ace_sample_scaled <- lapply(post_ord_ace_sample_scaled, lapply,  function(x){ get.disparity(dispRity(x, metric = c(sum, variances)))})


disp_post_ord_ace_sample_all <- lapply(post_ord_ace_sample_rm_axes_all, lapply,  function(x){ get.disparity(dispRity(x, metric = c(sum, variances)))})

disp_post_ord_ace_point <- lapply(post_ord_ace_point,  function(x){ get.disparity(dispRity(x, metric = c(sum, variances)))})

disp_no_ace <- lapply(ord_no_ace, lapply, lapply, function(x){ get.disparity(dispRity(x, metric = c(sum, variances)))})



rep_001_sample <- disp_post_ord_ace_sample[[1]]
rep_001_sample_diffs <- lapply(disp_post_ord_ace_sample[[1]], function(x){
    est_disp <- x[[1]][1]
    (est_disp - disp_true[[1]]$fast[[1]])/disp_true[[1]]$fast[[1]] 
})

rep_001_sample_diffs_all <- lapply(disp_post_ord_ace_sample_all[[1]], function(x){
    est_disp <- x[[1]][1]
    (est_disp - disp_true[[1]]$fast[[1]])/disp_true[[1]]$fast[[1]] 
})



scale_01 <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}