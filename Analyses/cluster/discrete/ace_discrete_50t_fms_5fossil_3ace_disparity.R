args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])

source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")

library(dispRity)

base_path <- "/mnt/parscratch/users/bip24cns/acedisparity/discrete/50t/"

job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")



write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}


ord_rel <- readRDS(paste0(base_path, sprintf("ord/8187111_ord_rel_%03d.rds", replicate_id)))
ord_sample <- readRDS(paste0(base_path, sprintf("ord/8187111_ord_sample_%03d.rds", replicate_id)))
ord_strict <- readRDS(paste0(base_path, sprintf("ord/8187111_ord_strict_%03d.rds", replicate_id)))
ord_no_ace <- readRDS(paste0(base_path, sprintf("ord/8187111_ord_no_ace_%03d.rds", replicate_id)))
ord_true <- readRDS(paste0(base_path, sprintf("ord/8187111_ord_true_%03d.rds", replicate_id)))
point_post_ord_ace <- readRDS(paste0(base_path, sprintf("ord/8187111_post_ord_point_%03d.rds", replicate_id)))
sample_post_ord_ace <- readRDS(paste0(base_path, sprintf("ord/8187111_post_ord_sample_%03d.rds", replicate_id)))

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

saveRDS(sum_var_rel, write.path("disparity/raw", "sumvar_rel_%03d.rds"))
saveRDS(sum_var_sample, write.path("disparity/raw", "sumvar_sample_%03d.rds"))
saveRDS(sum_var_strict, write.path("disparity/raw", "sumvar_strict_%03d.rds"))
saveRDS(sum_var_no_ace, write.path("disparity/raw", "sumvar_no_ace_%03d.rds"))
saveRDS(sum_var_true, write.path("disparity/raw", "sumvar_true_%03d.rds"))
saveRDS(sum_var_point_post_ord, write.path("disparity/raw", "sumvar_point_postord_%03d.rds"))
saveRDS(sum_var_sample_post_ord, write.path("disparity/raw", "sumvar_sample_postord_%03d.rds"))





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
     disp.diff(sample, rate_true)
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
     disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
  }, rate_sample)
}, sum_var_sample_post_ord, sum_var_true)



saveRDS(sum_var_rel_diff, write.path("disparity/raw", "diff_sumvar_rel_%03d.rds"))
saveRDS(sum_var_sample_diff, write.path("disparity/raw", "diff_sumvar_sample_%03d.rds"))
saveRDS(sum_var_strict_diff, write.path("disparity/raw", "diff_sumvar_strict_%03d.rds"))
saveRDS(sum_var_no_ace_diff, write.path("disparity/raw", "diff_sumvar_no_ace_%03d.rds"))
saveRDS(sum_var_post_ord_point_diff, write.path("disparity/raw", "diff_sumvar_point_postord_%03d.rds"))
saveRDS(sum_var_post_ord_sample_diff, write.path("disparity/raw", "diff_sumvar_sample_postord_%03d.rds"))


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


saveRDS(neighbours_rel, write.path("disparity/raw", "neighbours_rel_%03d.rds"))
saveRDS(neighbours_sample, write.path("disparity/raw", "neighbours_sample_%03d.rds"))
saveRDS(neighbours_strict, write.path("disparity/raw", "neighbours_strict_%03d.rds"))
saveRDS(neighbours_no_ace, write.path("disparity/raw", "neighbours_no_ace_%03d.rds"))
saveRDS(neighbours_true, write.path("disparity/raw", "neighbours_true_%03d.rds"))
saveRDS(neighbours_point_post_ord, write.path("disparity/raw", "neighbours_point_postord_%03d.rds"))
saveRDS(neighbours_sample_post_ord, write.path("disparity/raw", "neighbours_sample_postord_%03d.rds"))




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
     disp.diff(sample, rate_true)
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
     disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
  }, rate_sample)
}, neighbours_sample_post_ord, neighbours_true)

saveRDS(neighbours_rel_diff, write.path("disparity/raw", "diff_neighbours_rel_%03d.rds"))
saveRDS(neighbours_sample_diff, write.path("disparity/raw", "diff_neighbours_sample_%03d.rds"))
saveRDS(neighbours_strict_diff, write.path("disparity/raw", "diff_neighbours_strict_%03d.rds"))
saveRDS(neighbours_no_ace_diff, write.path("disparity/raw", "diff_neighbours_no_ace_%03d.rds"))
saveRDS(neighbours_point_postord_diff, write.path("disparity/raw", "diff_neighbours_point_postord_%03d.rds"))
saveRDS(neighbours_sample_postord_diff, write.path("disparity/raw", "diff_neighbours_sample_postord_%03d.rds"))



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

saveRDS(sum_quant_rel, write.path("disparity/raw", "quant_rel_%03d.rds"))
saveRDS(sum_quant_sample, write.path("disparity/raw", "quant_sample_%03d.rds"))
saveRDS(sum_quant_strict, write.path("disparity/raw", "quant_strict_%03d.rds"))
saveRDS(sum_quant_no_ace, write.path("disparity/raw", "quant_no_ace_%03d.rds"))
saveRDS(sum_quant_true, write.path("disparity/raw", "quant_true_%03d.rds"))
saveRDS(sum_quant_point_post_ord, write.path("disparity/raw", "quant_point_postord_%03d.rds"))
saveRDS(sum_quant_sample_post_ord, write.path("disparity/raw", "quant_sample_postord_%03d.rds"))





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
     disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
  }, rate_sample)
}, sum_quant_sample, sum_quant_true)

sum_quant_strict_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, sum_quant_strict, sum_quant_true)


sum_quant_post_ord_point_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, sum_quant_point_post_ord, sum_quant_true)

sum_quant_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
     disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
  }, rate_sample)
}, sum_quant_sample_post_ord, sum_quant_true)

saveRDS(sum_quant_rel_diff, write.path("disparity/raw", "diff_quant_rel_%03d.rds"))
saveRDS(sum_quant_sample_diff, write.path("disparity/raw", "diff_quant_sample_%03d.rds"))
saveRDS(sum_quant_strict_diff, write.path("disparity/raw", "diff_quant_strict_%03d.rds"))
saveRDS(sum_quant_no_ace_diff, write.path("disparity/raw", "diff_quant_no_ace_%03d.rds"))
saveRDS(sum_quant_post_ord_point_diff, write.path("disparity/raw", "diff_quant_point_postord_%03d.rds"))
saveRDS(sum_quant_post_ord_sample_diff, write.path("disparity/raw", "diff_quant_sample_postord_%03d.rds"))



#########################################################################################################################                                                                                                                      #
#                                                                                                                      #
#                                                                                                                      #
#                                            RESCALED DISPARITY                                                        #
#                                                                                                                      #
#                                                                                                                      #
########################################################################################################################                                                                                                                     

scale.pc <- function(ord){
    pc1_min <- min(ord[,1])
    pc1_range <- diff(range(ord[,1]))
    scal <- (ord - pc1_min) / pc1_range  # Shift then scale
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

saveRDS(rescaled_sum_var_rel, write.path("disparity/rescaled", "sumvar_rel_%03d.rds"))
saveRDS(rescaled_sum_var_sample, write.path("disparity/rescaled", "sumvar_sample_%03d.rds"))
saveRDS(rescaled_sum_var_strict, write.path("disparity/rescaled", "sumvar_strict_%03d.rds"))
saveRDS(rescaled_sum_var_no_ace, write.path("disparity/rescaled", "sumvar_no_ace_%03d.rds"))
saveRDS(rescaled_sum_var_true, write.path("disparity/rescaled", "sumvar_true_%03d.rds"))
saveRDS(rescaled_sum_var_point_post_ord, write.path("disparity/rescaled", "sumvar_point_postord_%03d.rds"))
saveRDS(rescaled_sum_var_sample_post_ord, write.path("disparity/rescaled", "sumvar_sample_postord_%03d.rds"))


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

saveRDS(rescaled_sum_var_rel_diff, write.path("disparity/rescaled", "diff_sumvar_rel_%03d.rds"))
saveRDS(rescaled_sum_var_sample_diff, write.path("disparity/rescaled", "diff_sumvar_sample_%03d.rds"))
saveRDS(rescaled_sum_var_strict_diff, write.path("disparity/rescaled", "diff_sumvar_strict_%03d.rds"))
saveRDS(rescaled_sum_var_no_ace_diff, write.path("disparity/rescaled", "diff_sumvar_no_ace_%03d.rds"))
saveRDS(rescaled_sum_var_post_ord_point_diff, write.path("disparity/rescaled", "diff_sumvar_point_postord_%03d.rds"))
saveRDS(rescaled_sum_var_post_ord_sample_diff, write.path("disparity/rescaled", "diff_sumvar_sample_postord_%03d.rds"))


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

saveRDS(rescaled_neighbours_rel, write.path("disparity/rescaled", "neighbours_rel_%03d.rds"))
saveRDS(rescaled_neighbours_sample, write.path("disparity/rescaled", "neighbours_sample_%03d.rds"))
saveRDS(rescaled_neighbours_strict, write.path("disparity/rescaled", "neighbours_strict_%03d.rds"))
saveRDS(rescaled_neighbours_no_ace, write.path("disparity/rescaled", "neighbours_no_ace_%03d.rds"))
saveRDS(rescaled_neighbours_true, write.path("disparity/rescaled", "neighbours_true_%03d.rds"))
saveRDS(rescaled_neighbours_point_post_ord, write.path("disparity/rescaled", "neighbours_point_postord_%03d.rds"))
saveRDS(rescaled_neighbours_sample_post_ord, write.path("disparity/rescaled", "neighbours_sample_postord_%03d.rds"))


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

saveRDS(rescaled_neighbours_rel_diff, write.path("disparity/rescaled", "diff_neighbours_rel_%03d.rds"))
saveRDS(rescaled_neighbours_sample_diff, write.path("disparity/rescaled", "diff_neighbours_sample_%03d.rds"))
saveRDS(rescaled_neighbours_strict_diff, write.path("disparity/rescaled", "diff_neighbours_strict_%03d.rds"))
saveRDS(rescaled_neighbours_no_ace_diff, write.path("disparity/rescaled", "diff_neighbours_no_ace_%03d.rds"))
saveRDS(rescaled_neighbours_point_postord_diff, write.path("disparity/rescaled", "diff_neighbours_point_postord_%03d.rds"))
saveRDS(rescaled_neighbours_sample_postord_diff, write.path("disparity/rescaled", "diff_neighbours_sample_postord_%03d.rds"))


rescaled_sum_quant_rel <- lapply(rescaled_rel, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rescaled_sum_quant_sample <- lapply(rescaled_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rescaled_sum_quant_strict <- lapply(rescaled_strict, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rescaled_sum_quant_no_ace <- lapply(rescaled_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rescaled_sum_quant_true <- lapply(rescaled_true, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rescaled_sum_quant_point_post_ord <- lapply(rescaled_point_postord, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rescaled_sum_quant_sample_post_ord <- lapply(rescaled_sample_postord, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

saveRDS(rescaled_sum_quant_rel, write.path("disparity/rescaled", "quant_rel_%03d.rds"))
saveRDS(rescaled_sum_quant_sample, write.path("disparity/rescaled", "quant_sample_%03d.rds"))
saveRDS(rescaled_sum_quant_strict, write.path("disparity/rescaled", "quant_strict_%03d.rds"))
saveRDS(rescaled_sum_quant_no_ace, write.path("disparity/rescaled", "quant_no_ace_%03d.rds"))
saveRDS(rescaled_sum_quant_true, write.path("disparity/rescaled", "quant_true_%03d.rds"))
saveRDS(rescaled_sum_quant_point_post_ord, write.path("disparity/rescaled", "quant_point_postord_%03d.rds"))
saveRDS(rescaled_sum_quant_sample_post_ord, write.path("disparity/rescaled", "quant_sample_postord_%03d.rds"))

rescaled_sum_quant_rel_diff <- Map(function(rate_rel, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_rel)
}, rescaled_sum_quant_rel, rescaled_sum_quant_true)

rescaled_sum_quant_no_ace_diff <- Map(function(rate_no_ace, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_no_ace)
}, rescaled_sum_quant_no_ace, rescaled_sum_quant_true)

rescaled_sum_quant_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rescaled_sum_quant_sample, rescaled_sum_quant_true)

rescaled_sum_quant_strict_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, rescaled_sum_quant_strict, rescaled_sum_quant_true)

rescaled_sum_quant_post_ord_point_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, rescaled_sum_quant_point_post_ord, rescaled_sum_quant_true)

rescaled_sum_quant_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rescaled_sum_quant_sample_post_ord, rescaled_sum_quant_true)

saveRDS(rescaled_sum_quant_rel_diff, write.path("disparity/rescaled", "diff_quant_rel_%03d.rds"))
saveRDS(rescaled_sum_quant_sample_diff, write.path("disparity/rescaled", "diff_quant_sample_%03d.rds"))
saveRDS(rescaled_sum_quant_strict_diff, write.path("disparity/rescaled", "diff_quant_strict_%03d.rds"))
saveRDS(rescaled_sum_quant_no_ace_diff, write.path("disparity/rescaled", "diff_quant_no_ace_%03d.rds"))
saveRDS(rescaled_sum_quant_post_ord_point_diff, write.path("disparity/rescaled", "diff_quant_point_postord_%03d.rds"))
saveRDS(rescaled_sum_quant_post_ord_sample_diff, write.path("disparity/rescaled", "diff_quant_sample_postord_%03d.rds"))

#########################################################################################################################                                                                                                                      #
#                                                                                                                      #
#                                                                                                                      #
#                                             REMOVE AXES DISPARITY                                                    #
#                                                                                                                      #
#                                                                                                                      #
########################################################################################################################

remove.axes <- function(ord){
    select_axes <- ord[,1:48]
    return(select_axes)
}

rmaxes_rel <- lapply(ord_rel, lapply, remove.axes)
rmaxes_sample <- lapply(ord_sample, lapply, lapply, remove.axes)
rmaxes_strict <- lapply(ord_strict, lapply, remove.axes)
rmaxes_no_ace <- lapply(ord_no_ace, lapply, remove.axes)
rmaxes_true <- lapply(ord_true, remove.axes)
rmaxes_point_postord <- lapply(point_post_ord_ace, lapply, remove.axes)
rmaxes_sample_postord <- lapply(sample_post_ord_ace, lapply, lapply, remove.axes)

# Sum of variances
rmaxes_sum_var_rel <- lapply(rmaxes_rel, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rmaxes_sum_var_sample <- lapply(rmaxes_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rmaxes_sum_var_strict <- lapply(rmaxes_strict, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rmaxes_sum_var_no_ace <- lapply(rmaxes_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rmaxes_sum_var_true <- lapply(rmaxes_true, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rmaxes_sum_var_point_post_ord <- lapply(rmaxes_point_postord, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rmaxes_sum_var_sample_post_ord <- lapply(rmaxes_sample_postord, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

saveRDS(rmaxes_sum_var_rel, write.path("disparity/rmaxes", "sumvar_rel_%03d.rds"))
saveRDS(rmaxes_sum_var_sample, write.path("disparity/rmaxes", "sumvar_sample_%03d.rds"))
saveRDS(rmaxes_sum_var_strict, write.path("disparity/rmaxes", "sumvar_strict_%03d.rds"))
saveRDS(rmaxes_sum_var_no_ace, write.path("disparity/rmaxes", "sumvar_no_ace_%03d.rds"))
saveRDS(rmaxes_sum_var_true, write.path("disparity/rmaxes", "sumvar_true_%03d.rds"))
saveRDS(rmaxes_sum_var_point_post_ord, write.path("disparity/rmaxes", "sumvar_point_postord_%03d.rds"))
saveRDS(rmaxes_sum_var_sample_post_ord, write.path("disparity/rmaxes", "sumvar_sample_postord_%03d.rds"))

# Diffs
rmaxes_sum_var_rel_diff <- Map(function(rate_rel, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_rel)
}, rmaxes_sum_var_rel, rmaxes_sum_var_true)

rmaxes_sum_var_no_ace_diff <- Map(function(rate_no_ace, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_no_ace)
}, rmaxes_sum_var_no_ace, rmaxes_sum_var_true)

rmaxes_sum_var_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rmaxes_sum_var_sample, rmaxes_sum_var_true)

rmaxes_sum_var_strict_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, rmaxes_sum_var_strict, rmaxes_sum_var_true)

rmaxes_sum_var_post_ord_point_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, rmaxes_sum_var_point_post_ord, rmaxes_sum_var_true)

rmaxes_sum_var_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rmaxes_sum_var_sample_post_ord, rmaxes_sum_var_true)

saveRDS(rmaxes_sum_var_rel_diff, write.path("disparity/rmaxes", "diff_sumvar_rel_%03d.rds"))
saveRDS(rmaxes_sum_var_sample_diff, write.path("disparity/rmaxes", "diff_sumvar_sample_%03d.rds"))
saveRDS(rmaxes_sum_var_strict_diff, write.path("disparity/rmaxes", "diff_sumvar_strict_%03d.rds"))
saveRDS(rmaxes_sum_var_no_ace_diff, write.path("disparity/rmaxes", "diff_sumvar_no_ace_%03d.rds"))
saveRDS(rmaxes_sum_var_post_ord_point_diff, write.path("disparity/rmaxes", "diff_sumvar_point_postord_%03d.rds"))
saveRDS(rmaxes_sum_var_post_ord_sample_diff, write.path("disparity/rmaxes", "diff_sumvar_sample_postord_%03d.rds"))

# Neighbours
rmaxes_neighbours_rel <- lapply(rmaxes_rel, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

rmaxes_neighbours_sample <- lapply(rmaxes_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

rmaxes_neighbours_strict <- lapply(rmaxes_strict, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

rmaxes_neighbours_no_ace <- lapply(rmaxes_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

rmaxes_neighbours_true <- lapply(rmaxes_true, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

rmaxes_neighbours_point_post_ord <- lapply(rmaxes_point_postord, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

rmaxes_neighbours_sample_post_ord <- lapply(rmaxes_sample_postord, lapply, lapply, function(rep){
    dispRity(rep, metric = c(mean, neighbours))$disparity
})

saveRDS(rmaxes_neighbours_rel, write.path("disparity/rmaxes", "neighbours_rel_%03d.rds"))
saveRDS(rmaxes_neighbours_sample, write.path("disparity/rmaxes", "neighbours_sample_%03d.rds"))
saveRDS(rmaxes_neighbours_strict, write.path("disparity/rmaxes", "neighbours_strict_%03d.rds"))
saveRDS(rmaxes_neighbours_no_ace, write.path("disparity/rmaxes", "neighbours_no_ace_%03d.rds"))
saveRDS(rmaxes_neighbours_true, write.path("disparity/rmaxes", "neighbours_true_%03d.rds"))
saveRDS(rmaxes_neighbours_point_post_ord, write.path("disparity/rmaxes", "neighbours_point_postord_%03d.rds"))
saveRDS(rmaxes_neighbours_sample_post_ord, write.path("disparity/rmaxes", "neighbours_sample_postord_%03d.rds"))

# Neighbours diffs
rmaxes_neighbours_rel_diff <- Map(function(rate_rel, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_rel)
}, rmaxes_neighbours_rel, rmaxes_neighbours_true)

rmaxes_neighbours_no_ace_diff <- Map(function(rate_no_ace, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_no_ace)
}, rmaxes_neighbours_no_ace, rmaxes_neighbours_true)

rmaxes_neighbours_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rmaxes_neighbours_sample, rmaxes_neighbours_true)

rmaxes_neighbours_strict_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, rmaxes_neighbours_strict, rmaxes_neighbours_true)

rmaxes_neighbours_point_postord_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, rmaxes_neighbours_point_post_ord, rmaxes_neighbours_true)

rmaxes_neighbours_sample_postord_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rmaxes_neighbours_sample_post_ord, rmaxes_neighbours_true)

saveRDS(rmaxes_neighbours_rel_diff, write.path("disparity/rmaxes", "diff_neighbours_rel_%03d.rds"))
saveRDS(rmaxes_neighbours_sample_diff, write.path("disparity/rmaxes", "diff_neighbours_sample_%03d.rds"))
saveRDS(rmaxes_neighbours_strict_diff, write.path("disparity/rmaxes", "diff_neighbours_strict_%03d.rds"))
saveRDS(rmaxes_neighbours_no_ace_diff, write.path("disparity/rmaxes", "diff_neighbours_no_ace_%03d.rds"))
saveRDS(rmaxes_neighbours_point_postord_diff, write.path("disparity/rmaxes", "diff_neighbours_point_postord_%03d.rds"))
saveRDS(rmaxes_neighbours_sample_postord_diff, write.path("disparity/rmaxes", "diff_neighbours_sample_postord_%03d.rds"))

# Add this after the rmaxes neighbours section:

# Sum of quantiles
rmaxes_sum_quant_rel <- lapply(rmaxes_rel, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rmaxes_sum_quant_sample <- lapply(rmaxes_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rmaxes_sum_quant_strict <- lapply(rmaxes_strict, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rmaxes_sum_quant_no_ace <- lapply(rmaxes_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rmaxes_sum_quant_true <- lapply(rmaxes_true, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rmaxes_sum_quant_point_post_ord <- lapply(rmaxes_point_postord, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rmaxes_sum_quant_sample_post_ord <- lapply(rmaxes_sample_postord, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

saveRDS(rmaxes_sum_quant_rel, write.path("disparity/rmaxes", "quant_rel_%03d.rds"))
saveRDS(rmaxes_sum_quant_sample, write.path("disparity/rmaxes", "quant_sample_%03d.rds"))
saveRDS(rmaxes_sum_quant_strict, write.path("disparity/rmaxes", "quant_strict_%03d.rds"))
saveRDS(rmaxes_sum_quant_no_ace, write.path("disparity/rmaxes", "quant_no_ace_%03d.rds"))
saveRDS(rmaxes_sum_quant_true, write.path("disparity/rmaxes", "quant_true_%03d.rds"))
saveRDS(rmaxes_sum_quant_point_post_ord, write.path("disparity/rmaxes", "quant_point_postord_%03d.rds"))
saveRDS(rmaxes_sum_quant_sample_post_ord, write.path("disparity/rmaxes", "quant_sample_postord_%03d.rds"))

# Sum of quantiles diffs
rmaxes_sum_quant_rel_diff <- Map(function(rate_rel, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_rel)
}, rmaxes_sum_quant_rel, rmaxes_sum_quant_true)

rmaxes_sum_quant_no_ace_diff <- Map(function(rate_no_ace, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_no_ace)
}, rmaxes_sum_quant_no_ace, rmaxes_sum_quant_true)

rmaxes_sum_quant_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rmaxes_sum_quant_sample, rmaxes_sum_quant_true)

rmaxes_sum_quant_strict_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, rmaxes_sum_quant_strict, rmaxes_sum_quant_true)

rmaxes_sum_quant_post_ord_point_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, rmaxes_sum_quant_point_post_ord, rmaxes_sum_quant_true)

rmaxes_sum_quant_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rmaxes_sum_quant_sample_post_ord, rmaxes_sum_quant_true)

saveRDS(rmaxes_sum_quant_rel_diff, write.path("disparity/rmaxes", "diff_quant_rel_%03d.rds"))
saveRDS(rmaxes_sum_quant_sample_diff, write.path("disparity/rmaxes", "diff_quant_sample_%03d.rds"))
saveRDS(rmaxes_sum_quant_strict_diff, write.path("disparity/rmaxes", "diff_quant_strict_%03d.rds"))
saveRDS(rmaxes_sum_quant_no_ace_diff, write.path("disparity/rmaxes", "diff_quant_no_ace_%03d.rds"))
saveRDS(rmaxes_sum_quant_post_ord_point_diff, write.path("disparity/rmaxes", "diff_quant_point_postord_%03d.rds"))
saveRDS(rmaxes_sum_quant_post_ord_sample_diff, write.path("disparity/rmaxes", "diff_quant_sample_postord_%03d.rds"))