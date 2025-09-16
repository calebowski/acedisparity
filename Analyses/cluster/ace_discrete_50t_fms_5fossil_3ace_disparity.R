args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])

source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")

library(dispRity)

ord_rel <- readRDS(sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/ord/ord_rel_%03d.rds", replicate_id))
ord_sample <- readRDS(sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/ord/ord_sample_%03d.rds", replicate_id))
ord_strict <- readRDS(sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/ord/ord_strict_%03d.rds", replicate_id))
ord_no_ace <- readRDS(sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/ord/ord_no_ace_%03d.rds", replicate_id))
ord_true <- readRDS(sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/ord/ord_true_%03d.rds", replicate_id))
point_post_ord_ace <- readRDS(sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/ord/post_ord_ace_point_%03d.rds", replicate_id))
sample_post_ord_ace <- readRDS(sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/ord/post_ord_ace_sample_%03d.rds", replicate_id))

# ord_rel <- readRDS("../Data/cluster/discrete//ord/ord_rel_001.rds")
# ord_sample <- readRDS("../Data/cluster/discrete/ord/ord_sample_001.rds")
# ord_strict <- readRDS("../Data/cluster/discrete/ord/ord_strict_001.rds")
# ord_no_ace <- readRDS("../Data/cluster/discrete/ord/ord_no_ace_001.rds")
# ord_true <- readRDS("../Data/cluster/discrete/ord/ord_true_001.rds")

scale.pc <- function(ord){
    pc1 <- ord[,1]
    min <- abs(min(pc1))
    max <- abs(max(pc1)) + min
    scal <- (ord + min) / max
    return(scal)
}


sum_var_rel <- lapply(ord_rel, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_sample <- lapply(ord_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_strict <- lapply(ord_strict, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_no_ace <- lapply(ord_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_true <- lapply(ord_true,  function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_point_post_ord <- lapply(point_post_ord_ace, lapply,  function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_sample_post_ord <- lapply(sample_post_ord_ace, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

saveRDS(sum_var_rel, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/sumvar_rel_%03d.rds", replicate_id))
saveRDS(sum_var_sample, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/sumvar_sample_%03d.rds", replicate_id))
saveRDS(sum_var_strict, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/sumvar_strict_%03d.rds", replicate_id))
saveRDS(sum_var_no_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/sumvar_no_ace_%03d.rds", replicate_id))
saveRDS(sum_var_true, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/sumvar_true_%03d.rds", replicate_id))
saveRDS(sum_var_point_post_ord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/sumvar_point_postord_%03d.rds", replicate_id))
saveRDS(sum_var_sample_post_ord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/sumvar_sample_postord_%03d.rds", replicate_id))




sum_var_rel_diff <- Map(function(rate_rel, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_rel)
}, sum_var_rel, sum_var_true)


sum_var_no_ace_diff <- Map(function(rate_no_ace, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_no_ace)
}, sum_var_no_ace, sum_var_true)


sum_var_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
     (disp.diff(sample, rate_true))
      }, fossil_data)
      return(unlist(diffs))
  }, rate_sample)
}, sum_var_sample, sum_var_true)

sum_var_strict_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, sum_var_strict, sum_var_true)


sum_var_post_ord_point_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, sum_var_point_post_ord, sum_var_true)

sum_var_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
     (disp.diff(sample, rate_true))
      }, fossil_data)
      return(unlist(diffs))
  }, rate_sample)
}, sum_var_sample_post_ord, sum_var_true)



saveRDS(sum_var_rel_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_sumvar_rel_%03d.rds", replicate_id))
saveRDS(sum_var_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_sumvar_sample_%03d.rds", replicate_id))
saveRDS(sum_var_strict_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_sumvar_strict_%03d.rds", replicate_id))
saveRDS(sum_var_no_ace_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_sumvar_no_ace_%03d.rds", replicate_id))
saveRDS(sum_var_post_ord_point_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_sumvar_point_postord_%03d.rds", replicate_id))
saveRDS(sum_var_post_ord_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_sumvar_sample_postord_%03d.rds", replicate_id))


neighbours_rel <- lapply(ord_rel, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

neighbours_sample <- lapply(ord_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

neighbours_strict <- lapply(ord_strict, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

neighbours_no_ace <- lapply(ord_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

neighbours_true <- lapply(ord_true,  function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

neighbours_point_post_ord <- lapply(point_post_ord_ace, lapply,  function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

neighbours_sample_post_ord <- lapply(sample_post_ord_ace, lapply, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})


saveRDS(neighbours_rel, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/neighbours_rel_%03d.rds", replicate_id))
saveRDS(neighbours_sample, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/neighbours_sample_%03d.rds", replicate_id))
saveRDS(neighbours_strict, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/neighbours_strict_%03d.rds", replicate_id))
saveRDS(neighbours_no_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/neighbours_no_ace_%03d.rds", replicate_id))
saveRDS(neighbours_true, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/neighbours_true_%03d.rds", replicate_id))
saveRDS(neighbours_point_post_ord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/neighbours_point_postord_%03d.rds", replicate_id))
saveRDS(neighbours_sample_post_ord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/neighbours_sample_postord_%03d.rds", replicate_id))



neighbours_rel_diff <- Map(function(rate_rel, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_rel)
}, neighbours_rel, neighbours_true)


neighbours_no_ace_diff <- Map(function(rate_no_ace, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_no_ace)
}, neighbours_no_ace, neighbours_true)


neighbours_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
     (disp.diff(sample, rate_true))
      }, fossil_data)
      return(unlist(diffs))
  }, rate_sample)
}, neighbours_sample, neighbours_true)

neighbours_strict_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, neighbours_strict, neighbours_true)

neighbours_point_postord_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, neighbours_point_post_ord, neighbours_true)


neighbours_sample_postord_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
     (disp.diff(sample, rate_true))
      }, fossil_data)
      return(unlist(diffs))
  }, rate_sample)
}, neighbours_sample_post_ord, neighbours_true)

saveRDS(neighbours_rel_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_neighbours_rel_%03d.rds", replicate_id))
saveRDS(neighbours_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_neighbours_sample_%03d.rds", replicate_id))
saveRDS(neighbours_strict_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_neighbours_strict_%03d.rds", replicate_id))
saveRDS(neighbours_no_ace_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_neighbours_no_ace_%03d.rds", replicate_id))
saveRDS(neighbours_point_postord_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_neighbours_point_postord_%03d.rds", replicate_id))
saveRDS(neighbours_sample_postord_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_neighbours_sample_postord_%03d.rds", replicate_id))