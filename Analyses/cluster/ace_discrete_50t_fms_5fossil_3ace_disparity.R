args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])

source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")

library(dispRity)

ord_rel <- readRDS(sprintf("/users/bip24cns/acedisparity/discrete/out/ord/ord_rel_%03d.rds", replicate_id))
ord_sample <- readRDS(sprintf("/users/bip24cns/acedisparity/discrete/out/ord/ord_sample_%03d.rds", replicate_id))
ord_strict <- readRDS(sprintf("/users/bip24cns/acedisparity/discrete/out/ord/ord_strict_%03d.rds", replicate_id))
ord_no_ace <- readRDS(sprintf("/users/bip24cns/acedisparity/discrete/out/ord/ord_no_ace_%03d.rds", replicate_id))
ord_true <- readRDS(sprintf("/users/bip24cns/acedisparity/discrete/out/ord/ord_true_%03d.rds", replicate_id))

# ord_rel <- readRDS("../Data/cluster/discrete//ord/ord_rel_001.rds")
# ord_sample <- readRDS("../Data/cluster/discrete/ord/ord_sample_001.rds")
# ord_strict <- readRDS("../Data/cluster/discrete/ord/ord_strict_001.rds")
# ord_no_ace <- readRDS("../Data/cluster/discrete/ord/ord_no_ace_001.rds")
# ord_true <- readRDS("../Data/cluster/discrete/ord/ord_true_001.rds")




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

saveRDS(sum_var_rel, sprintf("/users/bip24cns/acedisparity/discrete/out/disparity/sumvar_rel_%03d.rds", replicate_id))
saveRDS(sum_var_sample, sprintf("/users/bip24cns/acedisparity/discrete/out/disparity/sumvar_sample_%03d.rds", replicate_id))
saveRDS(sum_var_strict, sprintf("/users/bip24cns/acedisparity/discrete/out/disparity/sumvar_strict_%03d.rds", replicate_id))
saveRDS(sum_var_no_ace, sprintf("/users/bip24cns/acedisparity/discrete/out/disparity/sumvar_no_ace_%03d.rds", replicate_id))
saveRDS(sum_var_true, sprintf("/users/bip24cns/acedisparity/discrete/out/disparity/sumvar_true_%03d.rds", replicate_id))



disp_rel_diff <- Map(function(rate_rel, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_rel)
}, disp_rel, disp_true)


disp_no_ace_diff <- Map(function(rate_no_ace, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_no_ace)
}, disp_no_ace, disp_true)


disp_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
     (disp.diff(sample, rate_true))
      }, fossil_data)
      return(unlist(diffs))
  }, rate_sample)
}, disp_sample, disp_true)

disp_strict_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, disp_strict, disp_true)


saveRDS(disp_rel_diff, sprintf("/users/bip24cns/acedisparity/discrete/out/diff_sumvar_rel_%03d.rds", replicate_id))
saveRDS(disp_sample_diff, sprintf("/users/bip24cns/acedisparity/discrete/out/diff_sumvar_sample_%03d.rds", replicate_id))
saveRDS(disp_strict_diff, sprintf("/users/bip24cns/acedisparity/discrete/out/diff_sumvar_strict_%03d.rds", replicate_id))
saveRDS(disp_no_ace_diff, sprintf("/users/bip24cns/acedisparity/discrete/out/diff_sumvar_no_ace_%03d.rds", replicate_id))


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

saveRDS(neighbours_rel, sprintf("/users/bip24cns/acedisparity/discrete/out/disparity/neighbours_rel_%03d.rds", replicate_id))
saveRDS(neighbours_sample, sprintf("/users/bip24cns/acedisparity/discrete/out/disparity/neighbours_sample_%03d.rds", replicate_id))
saveRDS(neighbours_strict, sprintf("/users/bip24cns/acedisparity/discrete/out/disparity/neighbours_strict_%03d.rds", replicate_id))
saveRDS(neighbours_no_ace, sprintf("/users/bip24cns/acedisparity/discrete/out/disparity/neighbours_no_ace_%03d.rds", replicate_id))
saveRDS(neighbours_true, sprintf("/users/bip24cns/acedisparity/discrete/out/disparity/neighbours_true_%03d.rds", replicate_id))


neighbours_rel_diff <- Map(function(rate_rel, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_rel)
}, disp_rel, disp_true)


neighbours_no_ace_diff <- Map(function(rate_no_ace, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_no_ace)
}, disp_no_ace, disp_true)


neighbours_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
     (disp.diff(sample, rate_true))
      }, fossil_data)
      return(unlist(diffs))
  }, rate_sample)
}, disp_sample, disp_true)

neighbours_strict_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, disp_strict, disp_true)

saveRDS(disp_rel_diff, sprintf("/users/bip24cns/acedisparity/discrete/out/diff_neighbours_rel_%03d.rds", replicate_id))
saveRDS(disp_sample_diff, sprintf("/users/bip24cns/acedisparity/discrete/out/diff_neighbours_sample_%03d.rds", replicate_id))
saveRDS(disp_strict_diff, sprintf("/users/bip24cns/acedisparity/discrete/out/diff_neighbours_strict_%03d.rds", replicate_id))
saveRDS(disp_no_ace_diff, sprintf("/users/bip24cns/acedisparity/discrete/out/diff_neighbours_no_ace_%03d.rds", replicate_id))