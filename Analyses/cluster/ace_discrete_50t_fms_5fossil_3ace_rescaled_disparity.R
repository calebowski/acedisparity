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

scale.pc.remove.axes <- function(ord){
    pc1 <- ord[,1]
    min <- abs(min(pc1))
    max <- abs(max(pc1)) + min
    scal <- (ord + min) / max
    return(scal)
}

rescaled_rel <- lapply(ord_rel, lapply, scale.pc)
rescaled_sample <- lapply(ord_sample, lapply, lapply, scale.pc)
rescaled_strict <- lapply(ord_strict, lapply, scale.pc)
rescaled_no_ace <- lapply(ord_no_ace, lapply, scale.pc)
rescaled_true <- lapply(ord_true,  scale.pc)
rescaled_point_postord <- lapply(point_post_ord_ace, lapply,  scale.pc)
rescaled_sample_postord <- lapply(sample_post_ord_ace, lapply, lapply,  scale.pc)



sum_var_rel <- lapply(rescaled_rel, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_sample <- lapply(rescaled_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_strict <- lapply(rescaled_strict, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_no_ace <- lapply(rescaled_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_true <- lapply(rescaled_true,  function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_point_post_ord <- lapply(rescaled_point_postord, lapply,  function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_sample_post_ord <- lapply(rescaled_sample_postord, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

saveRDS(sum_var_rel, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/sumvar_rel_%03d.rds", replicate_id))
saveRDS(sum_var_sample, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/sumvar_sample_%03d.rds", replicate_id))
saveRDS(sum_var_strict, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/sumvar_strict_%03d.rds", replicate_id))
saveRDS(sum_var_no_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/sumvar_no_ace_%03d.rds", replicate_id))
saveRDS(sum_var_true, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/sumvar_true_%03d.rds", replicate_id))
saveRDS(sum_var_point_post_ord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/sumvar_point_postord_%03d.rds", replicate_id))
saveRDS(sum_var_sample_post_ord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/sumvar_sample_postord_%03d.rds", replicate_id))




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



saveRDS(sum_var_rel_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_sumvar_rel_%03d.rds", replicate_id))
saveRDS(sum_var_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_sumvar_sample_%03d.rds", replicate_id))
saveRDS(sum_var_strict_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_sumvar_strict_%03d.rds", replicate_id))
saveRDS(sum_var_no_ace_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_sumvar_no_ace_%03d.rds", replicate_id))
saveRDS(sum_var_post_ord_point_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_sumvar_point_postord_%03d.rds", replicate_id))
saveRDS(sum_var_post_ord_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_sumvar_sample_postord_%03d.rds", replicate_id))


neighbours_rel <- lapply(rescaled_rel, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

neighbours_sample <- lapply(rescaled_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

neighbours_strict <- lapply(rescaled_strict, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

neighbours_no_ace <- lapply(rescaled_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

neighbours_true <- lapply(rescaled_true,  function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

neighbours_point_post_ord <- lapply(rescaled_point_postord, lapply,  function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

neighbours_sample_post_ord <- lapply(rescaled_sample_postord, lapply, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})


saveRDS(neighbours_rel, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/neighbours_rel_%03d.rds", replicate_id))
saveRDS(neighbours_sample, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/neighbours_sample_%03d.rds", replicate_id))
saveRDS(neighbours_strict, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/neighbours_strict_%03d.rds", replicate_id))
saveRDS(neighbours_no_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/neighbours_no_ace_%03d.rds", replicate_id))
saveRDS(neighbours_true, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/neighbours_true_%03d.rds", replicate_id))
saveRDS(neighbours_point_post_ord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/neighbours_point_postord_%03d.rds", replicate_id))
saveRDS(neighbours_sample_post_ord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/neighbours_sample_postord_%03d.rds", replicate_id))



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

saveRDS(neighbours_rel_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_neighbours_rel_%03d.rds", replicate_id))
saveRDS(neighbours_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_neighbours_sample_%03d.rds", replicate_id))
saveRDS(neighbours_strict_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_neighbours_strict_%03d.rds", replicate_id))
saveRDS(neighbours_no_ace_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_neighbours_no_ace_%03d.rds", replicate_id))
saveRDS(neighbours_point_postord_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_neighbours_point_postord_%03d.rds", replicate_id))
saveRDS(neighbours_sample_postord_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_neighbours_sample_postord_%03d.rds", replicate_id))