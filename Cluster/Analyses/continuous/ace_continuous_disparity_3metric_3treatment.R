args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
tree_size <- args[2]

source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")
library(dispRity)

base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/continuous/", tree_size, "/")
job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")

write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}

# Load data

point_pre_ord_ace <- readRDS(paste0(base_path, sprintf("ord/8368504_ord_point_%03d.rds", replicate_id)))
sample_pre_ord_ace <- readRDS(paste0(base_path, sprintf("ord/8368504_ord_sample_%03d.rds", replicate_id)))
ord_no_ace <- readRDS(paste0(base_path, sprintf("ord/8368504_ord_no_ace_%03d.rds", replicate_id)))
ord_true <- readRDS(paste0(base_path, sprintf("ord/8368504_ord_true_%03d.rds", replicate_id)))
point_post_ord_ace <- readRDS(paste0(base_path, sprintf("ord/8368504_post_ord_point_%03d.rds", replicate_id)))
sample_post_ord_ace <- readRDS(paste0(base_path, sprintf("ord/8368504_post_ord_sample_%03d.rds", replicate_id)))

# Sum of variances
sum_var_sample <- lapply(sample_pre_ord_ace, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_point <- lapply(point_pre_ord_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_no_ace <- lapply(ord_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_true <- lapply(ord_true, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_point_post_ord <- lapply(point_post_ord_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

sum_var_sample_post_ord <- lapply(sample_post_ord_ace, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

saveRDS(sum_var_sample, write.path("disparity/raw", "sumvar_sample_preord_%03d.rds"))
saveRDS(sum_var_point, write.path("disparity/raw", "sumvar_point_preord_%03d.rds"))
saveRDS(sum_var_no_ace, write.path("disparity/raw", "sumvar_no_ace_%03d.rds"))
saveRDS(sum_var_true, write.path("disparity/raw", "sumvar_true_%03d.rds"))
saveRDS(sum_var_point_post_ord, write.path("disparity/raw", "sumvar_point_postord_%03d.rds"))
saveRDS(sum_var_sample_post_ord, write.path("disparity/raw", "sumvar_sample_postord_%03d.rds"))

# Sum of variances diffs
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

sum_var_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, sum_var_point, sum_var_true)

sum_var_post_ord_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, sum_var_point_post_ord, sum_var_true)

sum_var_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, sum_var_sample_post_ord, sum_var_true)

saveRDS(sum_var_sample_diff, write.path("disparity/raw", "diff_sumvar_sample_%03d.rds"))
saveRDS(sum_var_point_diff, write.path("disparity/raw", "diff_sumvar_point_%03d.rds"))
saveRDS(sum_var_no_ace_diff, write.path("disparity/raw", "diff_sumvar_no_ace_%03d.rds"))
saveRDS(sum_var_post_ord_point_diff, write.path("disparity/raw", "diff_sumvar_point_postord_%03d.rds"))
saveRDS(sum_var_post_ord_sample_diff, write.path("disparity/raw", "diff_sumvar_sample_postord_%03d.rds"))

# Neighbours
neighbours_sample <- lapply(sample_pre_ord_ace, lapply, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in neighbours_sample:", e$message, "\n")
        return(NA)
    })
})

neighbours_point <- lapply(point_pre_ord_ace, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in neighbours_point:", e$message, "\n")
        return(NA)
    })
})

neighbours_no_ace <- lapply(ord_no_ace, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in neighbours_no_ace:", e$message, "\n")
        return(NA)
    })
})

neighbours_true <- lapply(ord_true, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in neighbours_true:", e$message, "\n")
        return(NA)
    })
})

neighbours_point_post_ord <- lapply(point_post_ord_ace, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in neighbours_point_post_ord:", e$message, "\n")
        return(NA)
    })
})

neighbours_sample_post_ord <- lapply(sample_post_ord_ace, lapply, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in neighbours_sample_post_ord:", e$message, "\n")
        return(NA)
    })
})

saveRDS(neighbours_sample, write.path("disparity/raw", "neighbours_sample_%03d.rds"))
saveRDS(neighbours_point, write.path("disparity/raw", "neighbours_point_%03d.rds"))
saveRDS(neighbours_no_ace, write.path("disparity/raw", "neighbours_no_ace_%03d.rds"))
saveRDS(neighbours_true, write.path("disparity/raw", "neighbours_true_%03d.rds"))
saveRDS(neighbours_point_post_ord, write.path("disparity/raw", "neighbours_point_postord_%03d.rds"))
saveRDS(neighbours_sample_post_ord, write.path("disparity/raw", "neighbours_sample_postord_%03d.rds"))

# Neighbours diffs
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

neighbours_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, neighbours_point, neighbours_true)

neighbours_point_postord_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, neighbours_point_post_ord, neighbours_true)

neighbours_sample_postord_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, neighbours_sample_post_ord, neighbours_true)

saveRDS(neighbours_sample_diff, write.path("disparity/raw", "diff_neighbours_sample_%03d.rds"))
saveRDS(neighbours_point_diff, write.path("disparity/raw", "diff_neighbours_point_%03d.rds"))
saveRDS(neighbours_no_ace_diff, write.path("disparity/raw", "diff_neighbours_no_ace_%03d.rds"))
saveRDS(neighbours_point_postord_diff, write.path("disparity/raw", "diff_neighbours_point_postord_%03d.rds"))
saveRDS(neighbours_sample_postord_diff, write.path("disparity/raw", "diff_neighbours_sample_postord_%03d.rds"))

# Sum of quantiles
sum_quant_sample <- lapply(sample_pre_ord_ace, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

sum_quant_point <- lapply(point_pre_ord_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

sum_quant_no_ace <- lapply(ord_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

sum_quant_true <- lapply(ord_true, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

sum_quant_point_post_ord <- lapply(point_post_ord_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

sum_quant_sample_post_ord <- lapply(sample_post_ord_ace, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

saveRDS(sum_quant_sample, write.path("disparity/raw", "quant_sample_%03d.rds"))
saveRDS(sum_quant_point, write.path("disparity/raw", "quant_point_%03d.rds"))
saveRDS(sum_quant_no_ace, write.path("disparity/raw", "quant_no_ace_%03d.rds"))
saveRDS(sum_quant_true, write.path("disparity/raw", "quant_true_%03d.rds"))
saveRDS(sum_quant_point_post_ord, write.path("disparity/raw", "quant_point_postord_%03d.rds"))
saveRDS(sum_quant_sample_post_ord, write.path("disparity/raw", "quant_sample_postord_%03d.rds"))

# Sum of quantiles diffs
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

sum_quant_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, sum_quant_point, sum_quant_true)

sum_quant_post_ord_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, sum_quant_point_post_ord, sum_quant_true)

sum_quant_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, sum_quant_sample_post_ord, sum_quant_true)

saveRDS(sum_quant_sample_diff, write.path("disparity/raw", "diff_quant_sample_%03d.rds"))
saveRDS(sum_quant_point_diff, write.path("disparity/raw", "diff_quant_point_%03d.rds"))
saveRDS(sum_quant_no_ace_diff, write.path("disparity/raw", "diff_quant_no_ace_%03d.rds"))
saveRDS(sum_quant_post_ord_point_diff, write.path("disparity/raw", "diff_quant_point_postord_%03d.rds"))
saveRDS(sum_quant_post_ord_sample_diff, write.path("disparity/raw", "diff_quant_sample_postord_%03d.rds"))

#########################################################################################################################
#                                                                                                                      #
#                                            RESCALED DISPARITY                                                        #
#                                                                                                                      #
#########################################################################################################################

scale.pc <- function(ord){
    pc1 <- ord[,1]
    min <- abs(min(pc1))
    max <- abs(max(pc1)) + min
    scal <- (ord + min) / max
    return(scal)
}

rescaled_sample <- lapply(sample_pre_ord_ace, lapply, lapply, scale.pc)
rescaled_point <- lapply(point_pre_ord_ace, lapply, scale.pc)
rescaled_no_ace <- lapply(ord_no_ace, lapply, scale.pc)
rescaled_true <- lapply(ord_true, scale.pc)
rescaled_point_postord <- lapply(point_post_ord_ace, lapply, scale.pc)
rescaled_sample_postord <- lapply(sample_post_ord_ace, lapply, lapply, scale.pc)

# Sum of variances
rescaled_sum_var_sample <- lapply(rescaled_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rescaled_sum_var_point <- lapply(rescaled_point, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rescaled_sum_var_no_ace <- lapply(rescaled_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rescaled_sum_var_true <- lapply(rescaled_true, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rescaled_sum_var_point_post_ord <- lapply(rescaled_point_postord, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rescaled_sum_var_sample_post_ord <- lapply(rescaled_sample_postord, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

saveRDS(rescaled_sum_var_sample, write.path("disparity/rescaled", "sumvar_sample_%03d.rds"))
saveRDS(rescaled_sum_var_point, write.path("disparity/rescaled", "sumvar_point_%03d.rds"))
saveRDS(rescaled_sum_var_no_ace, write.path("disparity/rescaled", "sumvar_no_ace_%03d.rds"))
saveRDS(rescaled_sum_var_true, write.path("disparity/rescaled", "sumvar_true_%03d.rds"))
saveRDS(rescaled_sum_var_point_post_ord, write.path("disparity/rescaled", "sumvar_point_postord_%03d.rds"))
saveRDS(rescaled_sum_var_sample_post_ord, write.path("disparity/rescaled", "sumvar_sample_postord_%03d.rds"))

# Rescaled sum of variances diffs
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

rescaled_sum_var_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, rescaled_sum_var_point, rescaled_sum_var_true)

rescaled_sum_var_post_ord_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, rescaled_sum_var_point_post_ord, rescaled_sum_var_true)

rescaled_sum_var_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rescaled_sum_var_sample_post_ord, rescaled_sum_var_true)

saveRDS(rescaled_sum_var_sample_diff, write.path("disparity/rescaled", "diff_sumvar_sample_%03d.rds"))
saveRDS(rescaled_sum_var_point_diff, write.path("disparity/rescaled", "diff_sumvar_point_%03d.rds"))
saveRDS(rescaled_sum_var_no_ace_diff, write.path("disparity/rescaled", "diff_sumvar_no_ace_%03d.rds"))
saveRDS(rescaled_sum_var_post_ord_point_diff, write.path("disparity/rescaled", "diff_sumvar_point_postord_%03d.rds"))
saveRDS(rescaled_sum_var_post_ord_sample_diff, write.path("disparity/rescaled", "diff_sumvar_sample_postord_%03d.rds"))

# Neighbours
rescaled_neighbours_sample <- lapply(rescaled_sample, lapply, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in rescaled_neighbours_sample:", e$message, "\n")
        return(NA)
    })
})

rescaled_neighbours_point <- lapply(rescaled_point, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in rescaled_neighbours_point:", e$message, "\n")
        return(NA)
    })
})

rescaled_neighbours_no_ace <- lapply(rescaled_no_ace, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in rescaled_neighbours_no_ace:", e$message, "\n")
        return(NA)
    })
})

rescaled_neighbours_true <- lapply(rescaled_true, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in rescaled_neighbours_true:", e$message, "\n")
        return(NA)
    })
})

rescaled_neighbours_point_post_ord <- lapply(rescaled_point_postord, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in rescaled_neighbours_point_post_ord:", e$message, "\n")
        return(NA)
    })
})

rescaled_neighbours_sample_post_ord <- lapply(rescaled_sample_postord, lapply, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in rescaled_neighbours_sample_post_ord:", e$message, "\n")
        return(NA)
    })
})

saveRDS(rescaled_neighbours_sample, write.path("disparity/rescaled", "neighbours_sample_%03d.rds"))
saveRDS(rescaled_neighbours_point, write.path("disparity/rescaled", "neighbours_point_%03d.rds"))
saveRDS(rescaled_neighbours_no_ace, write.path("disparity/rescaled", "neighbours_no_ace_%03d.rds"))
saveRDS(rescaled_neighbours_true, write.path("disparity/rescaled", "neighbours_true_%03d.rds"))
saveRDS(rescaled_neighbours_point_post_ord, write.path("disparity/rescaled", "neighbours_point_postord_%03d.rds"))
saveRDS(rescaled_neighbours_sample_post_ord, write.path("disparity/rescaled", "neighbours_sample_postord_%03d.rds"))

# Neighbours diffs
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

rescaled_neighbours_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, rescaled_neighbours_point, rescaled_neighbours_true)

rescaled_neighbours_point_postord_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, rescaled_neighbours_point_post_ord, rescaled_neighbours_true)

rescaled_neighbours_sample_postord_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rescaled_neighbours_sample_post_ord, rescaled_neighbours_true)

saveRDS(rescaled_neighbours_sample_diff, write.path("disparity/rescaled", "diff_neighbours_sample_%03d.rds"))
saveRDS(rescaled_neighbours_point_diff, write.path("disparity/rescaled", "diff_neighbours_point_%03d.rds"))
saveRDS(rescaled_neighbours_no_ace_diff, write.path("disparity/rescaled", "diff_neighbours_no_ace_%03d.rds"))
saveRDS(rescaled_neighbours_point_postord_diff, write.path("disparity/rescaled", "diff_neighbours_point_postord_%03d.rds"))
saveRDS(rescaled_neighbours_sample_postord_diff, write.path("disparity/rescaled", "diff_neighbours_sample_postord_%03d.rds"))

# Sum of quantiles
rescaled_sum_quant_sample <- lapply(rescaled_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rescaled_sum_quant_point <- lapply(rescaled_point, lapply,  function(rep){
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

saveRDS(rescaled_sum_quant_sample, write.path("disparity/rescaled", "quant_sample_%03d.rds"))
saveRDS(rescaled_sum_quant_point, write.path("disparity/rescaled", "quant_point_%03d.rds"))
saveRDS(rescaled_sum_quant_no_ace, write.path("disparity/rescaled", "quant_no_ace_%03d.rds"))
saveRDS(rescaled_sum_quant_true, write.path("disparity/rescaled", "quant_true_%03d.rds"))
saveRDS(rescaled_sum_quant_point_post_ord, write.path("disparity/rescaled", "quant_point_postord_%03d.rds"))
saveRDS(rescaled_sum_quant_sample_post_ord, write.path("disparity/rescaled", "quant_sample_postord_%03d.rds"))

# Sum of quantiles diffs
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

rescaled_sum_quant_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, rescaled_sum_quant_point, rescaled_sum_quant_true)

rescaled_sum_quant_post_ord_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, rescaled_sum_quant_point_post_ord, rescaled_sum_quant_true)

rescaled_sum_quant_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rescaled_sum_quant_sample_post_ord, rescaled_sum_quant_true)

saveRDS(rescaled_sum_quant_sample_diff, write.path("disparity/rescaled", "diff_quant_sample_%03d.rds"))
saveRDS(rescaled_sum_quant_point_diff, write.path("disparity/rescaled", "diff_quant_point_%03d.rds"))
saveRDS(rescaled_sum_quant_no_ace_diff, write.path("disparity/rescaled", "diff_quant_no_ace_%03d.rds"))
saveRDS(rescaled_sum_quant_post_ord_point_diff, write.path("disparity/rescaled", "diff_quant_point_postord_%03d.rds"))
saveRDS(rescaled_sum_quant_post_ord_sample_diff, write.path("disparity/rescaled", "diff_quant_sample_postord_%03d.rds"))

#########################################################################################################################
#                                                                                                                      #
#                                             REMOVE AXES DISPARITY                                                    #
#                                                                                                                      #
#########################################################################################################################

remove.axes <- function(ord){
    select_axes <- ord[,1:48]
    return(select_axes)
}

rmaxes_sample <- lapply(sample_pre_ord_ace, lapply, lapply, remove.axes)
rmaxes_point <- lapply(point_pre_ord_ace, lapply, remove.axes)
rmaxes_no_ace <- lapply(ord_no_ace, lapply, remove.axes)
rmaxes_true <- lapply(ord_true, remove.axes)
rmaxes_point_postord <- lapply(point_post_ord_ace, lapply, remove.axes)
rmaxes_sample_postord <- lapply(sample_post_ord_ace, lapply, lapply, remove.axes)

# Sum of variances
rmaxes_sum_var_sample <- lapply(rmaxes_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, variances))$disparity
})

rmaxes_sum_var_point <- lapply(rmaxes_point, lapply, function(rep){
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

saveRDS(rmaxes_sum_var_sample, write.path("disparity/rmaxes", "sumvar_sample_%03d.rds"))
saveRDS(rmaxes_sum_var_point, write.path("disparity/rmaxes", "sumvar_point_%03d.rds"))
saveRDS(rmaxes_sum_var_no_ace, write.path("disparity/rmaxes", "sumvar_no_ace_%03d.rds"))
saveRDS(rmaxes_sum_var_true, write.path("disparity/rmaxes", "sumvar_true_%03d.rds"))
saveRDS(rmaxes_sum_var_point_post_ord, write.path("disparity/rmaxes", "sumvar_point_postord_%03d.rds"))
saveRDS(rmaxes_sum_var_sample_post_ord, write.path("disparity/rmaxes", "sumvar_sample_postord_%03d.rds"))

# Diffs
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

rmaxes_sum_var_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, rmaxes_sum_var_point, rmaxes_sum_var_true)

rmaxes_sum_var_post_ord_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, rmaxes_sum_var_point_post_ord, rmaxes_sum_var_true)

rmaxes_sum_var_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rmaxes_sum_var_sample_post_ord, rmaxes_sum_var_true)

saveRDS(rmaxes_sum_var_sample_diff, write.path("disparity/rmaxes", "diff_sumvar_sample_%03d.rds"))
saveRDS(rmaxes_sum_var_point_diff, write.path("disparity/rmaxes", "diff_sumvar_point_%03d.rds"))
saveRDS(rmaxes_sum_var_no_ace_diff, write.path("disparity/rmaxes", "diff_sumvar_no_ace_%03d.rds"))
saveRDS(rmaxes_sum_var_post_ord_point_diff, write.path("disparity/rmaxes", "diff_sumvar_point_postord_%03d.rds"))
saveRDS(rmaxes_sum_var_post_ord_sample_diff, write.path("disparity/rmaxes", "diff_sumvar_sample_postord_%03d.rds"))

# Neighbours
rmaxes_neighbours_sample <- lapply(rmaxes_sample, lapply, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in rmaxes_neighbours_sample:", e$message, "\n")
        return(NA)
    })
})

rmaxes_neighbours_point <- lapply(rmaxes_point, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in rmaxes_neighbours_point:", e$message, "\n")
        return(NA)
    })
})

rmaxes_neighbours_no_ace <- lapply(rmaxes_no_ace, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in rmaxes_neighbours_no_ace:", e$message, "\n")
        return(NA)
    })
})

rmaxes_neighbours_true <- lapply(rmaxes_true, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in rmaxes_neighbours_true:", e$message, "\n")
        return(NA)
    })
})

rmaxes_neighbours_point_post_ord <- lapply(rmaxes_point_postord, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in rmaxes_neighbours_point_post_ord:", e$message, "\n")
        return(NA)
    })
})

rmaxes_neighbours_sample_post_ord <- lapply(rmaxes_sample_postord, lapply, lapply, function(rep){
    tryCatch({
        dispRity(rep, metric = c(mean, neighbours))$disparity
    }, error = function(e) {
        cat("Error in rmaxes_neighbours_sample_post_ord:", e$message, "\n")
        return(NA)
    })
})

saveRDS(rmaxes_neighbours_sample, write.path("disparity/rmaxes", "neighbours_sample_%03d.rds"))
saveRDS(rmaxes_neighbours_point, write.path("disparity/rmaxes", "neighbours_point_%03d.rds"))
saveRDS(rmaxes_neighbours_no_ace, write.path("disparity/rmaxes", "neighbours_no_ace_%03d.rds"))
saveRDS(rmaxes_neighbours_true, write.path("disparity/rmaxes", "neighbours_true_%03d.rds"))
saveRDS(rmaxes_neighbours_point_post_ord, write.path("disparity/rmaxes", "neighbours_point_postord_%03d.rds"))
saveRDS(rmaxes_neighbours_sample_post_ord, write.path("disparity/rmaxes", "neighbours_sample_postord_%03d.rds"))

# Neighbours diffs
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

rmaxes_neighbours_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, rmaxes_neighbours_point, rmaxes_neighbours_true)

rmaxes_neighbours_point_postord_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, rmaxes_neighbours_point_post_ord, rmaxes_neighbours_true)

rmaxes_neighbours_sample_postord_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rmaxes_neighbours_sample_post_ord, rmaxes_neighbours_true)

saveRDS(rmaxes_neighbours_sample_diff, write.path("disparity/rmaxes", "diff_neighbours_sample_%03d.rds"))
saveRDS(rmaxes_neighbours_point_diff, write.path("disparity/rmaxes", "diff_neighbours_point_%03d.rds"))
saveRDS(rmaxes_neighbours_no_ace_diff, write.path("disparity/rmaxes", "diff_neighbours_no_ace_%03d.rds"))
saveRDS(rmaxes_neighbours_point_postord_diff, write.path("disparity/rmaxes", "diff_neighbours_point_postord_%03d.rds"))
saveRDS(rmaxes_neighbours_sample_postord_diff, write.path("disparity/rmaxes", "diff_neighbours_sample_postord_%03d.rds"))

# Sum of quantiles
rmaxes_sum_quant_sample <- lapply(rmaxes_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rmaxes_sum_quant_point <- lapply(rmaxes_point, lapply,  function(rep){
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

saveRDS(rmaxes_sum_quant_sample, write.path("disparity/rmaxes", "quant_sample_%03d.rds"))
saveRDS(rmaxes_sum_quant_point, write.path("disparity/rmaxes", "quant_point_%03d.rds"))
saveRDS(rmaxes_sum_quant_no_ace, write.path("disparity/rmaxes", "quant_no_ace_%03d.rds"))
saveRDS(rmaxes_sum_quant_true, write.path("disparity/rmaxes", "quant_true_%03d.rds"))
saveRDS(rmaxes_sum_quant_point_post_ord, write.path("disparity/rmaxes", "quant_point_postord_%03d.rds"))
saveRDS(rmaxes_sum_quant_sample_post_ord, write.path("disparity/rmaxes", "quant_sample_postord_%03d.rds"))

# Sum of quantiles diffs
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

rmaxes_sum_quant_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, rmaxes_sum_quant_point, rmaxes_sum_quant_true)

rmaxes_sum_quant_post_ord_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, rmaxes_sum_quant_point_post_ord, rmaxes_sum_quant_true)

rmaxes_sum_quant_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        disp.diff(sample, rate_true)
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rmaxes_sum_quant_sample_post_ord, rmaxes_sum_quant_true)

saveRDS(rmaxes_sum_quant_sample_diff, write.path("disparity/rmaxes", "diff_quant_sample_%03d.rds"))
saveRDS(rmaxes_sum_quant_point_diff, write.path("disparity/rmaxes", "diff_quant_point_%03d.rds"))
saveRDS(rmaxes_sum_quant_no_ace_diff, write.path("disparity/rmaxes", "diff_quant_no_ace_%03d.rds"))
saveRDS(rmaxes_sum_quant_post_ord_point_diff, write.path("disparity/rmaxes", "diff_quant_point_postord_%03d.rds"))
saveRDS(rmaxes_sum_quant_post_ord_sample_diff, write.path("disparity/rmaxes", "diff_quant_sample_postord_%03d.rds"))

cat("Disparity calculations completed for replicate", replicate_id, "\n")