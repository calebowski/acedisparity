args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])

source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")

library(dispRity)

base_path <- "/mnt/parscratch/users/bip24cns/acedisparity/discrete/50t/"

ord_rel <- readRDS(paste0(base_path, sprintf("ord/ord_rel_%03d.rds", replicate_id)))
ord_sample <- readRDS(paste0(base_path, sprintf("ord/ord_sample_%03d.rds", replicate_id)))
ord_strict <- readRDS(paste0(base_path, sprintf("ord/ord_strict_%03d.rds", replicate_id)))
ord_no_ace <- readRDS(paste0(base_path, sprintf("ord/ord_no_ace_%03d.rds", replicate_id)))
ord_true <- readRDS(paste0(base_path, sprintf("ord/ord_true_%03d.rds", replicate_id)))
point_post_ord_ace <- readRDS(paste0(base_path, sprintf("ord/post_ord_ace_point_%03d.rds", replicate_id)))
sample_post_ord_ace <- readRDS(paste0(base_path, sprintf("ord/post_ord_ace_sample_%03d.rds", replicate_id)))

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



sum_quant_rel <- lapply(ord_rel, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

sum_quant_sample <- lapply(ord_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

sum_quant_strict <- lapply(ord_strict, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

sum_quant_no_ace <- lapply(ord_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

sum_quant_true <- lapply(ord_true,  function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

sum_quant_point_post_ord <- lapply(point_post_ord_ace, lapply,  function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

sum_quant_sample_post_ord <- lapply(sample_post_ord_ace, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

saveRDS(sum_quant_rel, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/quant_rel_%03d.rds", replicate_id))
saveRDS(sum_quant_sample, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/quant_sample_%03d.rds", replicate_id))
saveRDS(sum_quant_strict, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/quant_strict_%03d.rds", replicate_id))
saveRDS(sum_quant_no_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/quant_no_ace_%03d.rds", replicate_id))
saveRDS(sum_quant_true, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/quant_true_%03d.rds", replicate_id))
saveRDS(sum_quant_point_post_ord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/quant_point_postord_%03d.rds", replicate_id))
saveRDS(sum_quant_sample_post_ord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/quant_sample_postord_%03d.rds", replicate_id))




sum_quant_rel_diff <- Map(function(rate_rel, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_rel)
}, sum_quant_rel, sum_quant_true)


sum_quant_no_ace_diff <- Map(function(rate_no_ace, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_no_ace)
}, sum_quant_no_ace, sum_quant_true)


sum_quant_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
     (disp.diff(sample, rate_true))
      }, fossil_data)
      return(unlist(diffs))
  }, rate_sample)
}, sum_quant_sample, sum_quant_true)

sum_quant_strict_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, sum_var_strict, sum_quant_true)


sum_quant_post_ord_point_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, sum_quant_point_post_ord, sum_quant_true)

sum_quant_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
     (disp.diff(sample, rate_true))
      }, fossil_data)
      return(unlist(diffs))
  }, rate_sample)
}, sum_quant_sample_post_ord, sum_quant_true)

saveRDS(sum_quant_rel_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_quant_rel_%03d.rds", replicate_id))
saveRDS(sum_quant_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_quant_sample_%03d.rds", replicate_id))
saveRDS(sum_quant_strict_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_quant_strict_%03d.rds", replicate_id))
saveRDS(sum_quant_no_ace_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_quant_no_ace_%03d.rds", replicate_id))
saveRDS(sum_quant_post_ord_point_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_quant_point_postord_%03d.rds", replicate_id))
saveRDS(sum_quant_post_ord_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/diff_quant_sample_postord_%03d.rds", replicate_id))


#########################################################################################################################                                                                                                                      #
#                                                                                                                      #
#                                                                                                                      #
#                                            RESCALED DISPARITY                                                        #
#                                                                                                                      #
#                                                                                                                      #
########################################################################################################################                                                                                                                     

scale.pc <- function(ord){
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



# Sum of variances
rescaled_sum_var_rel <- lapply(rescaled_rel, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rescaled_sum_var_sample <- lapply(rescaled_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rescaled_sum_var_strict <- lapply(rescaled_strict, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rescaled_sum_var_no_ace <- lapply(rescaled_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rescaled_sum_var_true <- lapply(rescaled_true,  function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rescaled_sum_var_point_post_ord <- lapply(rescaled_point_postord, lapply,  function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rescaled_sum_var_sample_post_ord <- lapply(rescaled_sample_postord, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

saveRDS(rescaled_sum_var_rel, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/sumvar_rel_%03d.rds", replicate_id))
saveRDS(rescaled_sum_var_sample, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/sumvar_sample_%03d.rds", replicate_id))
saveRDS(rescaled_sum_var_strict, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/sumvar_strict_%03d.rds", replicate_id))
saveRDS(rescaled_sum_var_no_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/sumvar_no_ace_%03d.rds", replicate_id))
saveRDS(rescaled_sum_var_true, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/sumvar_true_%03d.rds", replicate_id))
saveRDS(rescaled_sum_var_point_post_ord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/sumvar_point_postord_%03d.rds", replicate_id))
saveRDS(rescaled_sum_var_sample_post_ord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/sumvar_sample_postord_%03d.rds", replicate_id))

# Diffs
rescaled_sum_var_rel_diff <- Map(function(rate_rel, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_rel)
}, rescaled_sum_var_rel, rescaled_sum_var_true)

rescaled_sum_var_no_ace_diff <- Map(function(rate_no_ace, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_no_ace)
}, rescaled_sum_var_no_ace, rescaled_sum_var_true)

rescaled_sum_var_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rescaled_sum_var_sample, rescaled_sum_var_true)

rescaled_sum_var_strict_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, rescaled_sum_var_strict, rescaled_sum_var_true)

rescaled_sum_var_post_ord_point_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, rescaled_sum_var_point_post_ord, rescaled_sum_var_true)

rescaled_sum_var_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rescaled_sum_var_sample_post_ord, rescaled_sum_var_true)

saveRDS(rescaled_sum_var_rel_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_sumvar_rel_%03d.rds", replicate_id))
saveRDS(rescaled_sum_var_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_sumvar_sample_%03d.rds", replicate_id))
saveRDS(rescaled_sum_var_strict_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_sumvar_strict_%03d.rds", replicate_id))
saveRDS(rescaled_sum_var_no_ace_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_sumvar_no_ace_%03d.rds", replicate_id))
saveRDS(rescaled_sum_var_post_ord_point_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_sumvar_point_postord_%03d.rds", replicate_id))
saveRDS(rescaled_sum_var_post_ord_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_sumvar_sample_postord_%03d.rds", replicate_id))

# Neighbours
rescaled_neighbours_rel <- lapply(rescaled_rel, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

rescaled_neighbours_sample <- lapply(rescaled_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

rescaled_neighbours_strict <- lapply(rescaled_strict, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

rescaled_neighbours_no_ace <- lapply(rescaled_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

rescaled_neighbours_true <- lapply(rescaled_true,  function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

rescaled_neighbours_point_post_ord <- lapply(rescaled_point_postord, lapply,  function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

rescaled_neighbours_sample_post_ord <- lapply(rescaled_sample_postord, lapply, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

saveRDS(rescaled_neighbours_rel, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/neighbours_rel_%03d.rds", replicate_id))
saveRDS(rescaled_neighbours_sample, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/neighbours_sample_%03d.rds", replicate_id))
saveRDS(rescaled_neighbours_strict, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/neighbours_strict_%03d.rds", replicate_id))
saveRDS(rescaled_neighbours_no_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/neighbours_no_ace_%03d.rds", replicate_id))
saveRDS(rescaled_neighbours_true, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/neighbours_true_%03d.rds", replicate_id))
saveRDS(rescaled_neighbours_point_post_ord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/neighbours_point_postord_%03d.rds", replicate_id))
saveRDS(rescaled_neighbours_sample_post_ord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/neighbours_sample_postord_%03d.rds", replicate_id))

# Neighbours diffs
rescaled_neighbours_rel_diff <- Map(function(rate_rel, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_rel)
}, rescaled_neighbours_rel, rescaled_neighbours_true)

rescaled_neighbours_no_ace_diff <- Map(function(rate_no_ace, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_no_ace)
}, rescaled_neighbours_no_ace, rescaled_neighbours_true)

rescaled_neighbours_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rescaled_neighbours_sample, rescaled_neighbours_true)

rescaled_neighbours_strict_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, rescaled_neighbours_strict, rescaled_neighbours_true)

rescaled_neighbours_point_postord_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, rescaled_neighbours_point_post_ord, rescaled_neighbours_true)

rescaled_neighbours_sample_postord_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rescaled_neighbours_sample_post_ord, rescaled_neighbours_true)

saveRDS(rescaled_neighbours_rel_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_neighbours_rel_%03d.rds", replicate_id))
saveRDS(rescaled_neighbours_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_neighbours_sample_%03d.rds", replicate_id))
saveRDS(rescaled_neighbours_strict_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_neighbours_strict_%03d.rds", replicate_id))
saveRDS(rescaled_neighbours_no_ace_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_neighbours_no_ace_%03d.rds", replicate_id))
saveRDS(rescaled_neighbours_point_postord_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_neighbours_point_postord_%03d.rds", replicate_id))
saveRDS(rescaled_neighbours_sample_postord_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled/diff_neighbours_sample_postord_%03d.rds", replicate_id))

