args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])

source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")

library(dispRity)

ord_rel <- readRDS(sprintf("/users/bip24cns/acedisparity/discrete/out/ord_rel_%03d.rds", replicate_id))
ord_sample <- readRDS(sprintf("/users/bip24cns/acedisparity/discrete/out/ord_sample_%03d.rds", replicate_id))
ord_strict <- readRDS(sprintf("/users/bip24cns/acedisparity/discrete/out/ord_strict_%03d.rds", replicate_id))
ord_no_ace <- readRDS(sprintf("/users/bip24cns/acedisparity/discrete/out/ord_no_ace_%03d.rds", replicate_id))
ord_true <- readRDS(sprintf("/users/bip24cns/acedisparity/discrete/out/ord_true_%03d.rds", replicate_id))

# ord_rel <- readRDS("../Data/cluster/ord_rel_001.rds")
# ord_sample <- readRDS("../Data/cluster/ord_sample_001.rds")
# ord_strict <- readRDS("../Data/cluster/ord_strict_001.rds")
# ord_no_ace <- readRDS("../Data/cluster/ord_no_ace_001.rds")
# ord_true <- readRDS("../Data/cluster/ord_true_001.rds")





disp_rel <- lapply(ord_rel, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

disp_sample <- lapply(ord_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

disp_strict <- lapply(ord_strict, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

disp_no_ace <- lapply(ord_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

disp_true <- lapply(ord_true,  function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

saveRDS(disp_rel, sprintf("/users/bip24cns/acedisparity/discrete/out/disp_rel_%03d.rds", replicate_id))
saveRDS(disp_sample, sprintf("/users/bip24cns/acedisparity/discrete/out/disp_sample_%03d.rds", replicate_id))
saveRDS(disp_strict, sprintf("/users/bip24cns/acedisparity/discrete/out/disp_strict_%03d.rds", replicate_id))
saveRDS(disp_no_ace, sprintf("/users/bip24cns/acedisparity/discrete/out/disp_no_ace_%03d.rds", replicate_id))
saveRDS(disp_true, sprintf("/users/bip24cns/acedisparity/discrete/out/disp_true_%03d.rds", replicate_id))





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


saveRDS(disp_rel_diff, sprintf("/users/bip24cns/acedisparity/discrete/out/diff_disp_rel_%03d.rds", replicate_id))
saveRDS(disp_sample_diff, sprintf("/users/bip24cns/acedisparity/discrete/out/diff_disp_sample_%03d.rds", replicate_id))
saveRDS(disp_strict_diff, sprintf("/users/bip24cns/acedisparity/discrete/out/diff_disp_strict_%03d.rds", replicate_id))
saveRDS(disp_no_ace_diff, sprintf("/users/bip24cns/acedisparity/discrete/out/diff_disp_no_ace_%03d.rds", replicate_id))
