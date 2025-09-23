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

###########################################################################################################


# ord_rel <- readRDS("../Data/cluster/discrete/ord/ord_rel_080.rds")
# ord_strict <- readRDS("../Data/cluster/discrete/ord/ord_strict_080.rds")
# ord_no_ace <- readRDS("../Data/cluster/discrete/ord/ord_no_ace_080.rds")
# ord_sample <- readRDS("../Data/cluster/discrete/ord/ord_sample_080.rds")
# point_post_ord_ace <- readRDS("../Data/cluster/discrete/ord/post_ord_ace_point_080.rds")
# sample_post_ord_ace <- readRDS("../Data/cluster/discrete/ord/post_ord_ace_sample_080.rds")
# ord_true <- readRDS("../Data/cluster/discrete/ord/ord_true_080.rds")

scale.pc <- function(ord){
    pc1 <- ord[,1]
    min <- abs(min(pc1))
    max <- abs(max(pc1)) + min
    scal <- (ord + min) / max
    return(scal)
}

remove.axes <- function(ord){
    select_axes <- ord[,1:48]
    return(select_axes)
}

rescaled_rel <- lapply(ord_rel, lapply, scale.pc)
rescaled_sample <- lapply(ord_sample, lapply, lapply, scale.pc)
rescaled_strict <- lapply(ord_strict, lapply, scale.pc)
rescaled_no_ace <- lapply(ord_no_ace, lapply, scale.pc)
rescaled_true <- lapply(ord_true,  scale.pc)
rescaled_point_postord <- lapply(point_post_ord_ace, lapply,  scale.pc)
rescaled_sample_postord <- lapply(sample_post_ord_ace, lapply, lapply,  scale.pc)


rm_rel <- lapply(ord_rel, lapply, remove.axes)
rm_sample <- lapply(ord_sample, lapply, lapply, remove.axes)
rm_strict <- lapply(ord_strict, lapply, remove.axes)
rm_no_ace <- lapply(ord_no_ace, lapply, remove.axes)
rm_true <- lapply(ord_true,  remove.axes)
rm_point_postord <- lapply(point_post_ord_ace, lapply,  remove.axes)
rm_sample_postord <- lapply(sample_post_ord_ace, lapply, lapply,  remove.axes)



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
}, sum_quant_strict, sum_quant_true)


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


# Save original difference results (add after line 115)
saveRDS(sum_quant_rel_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/quant_rel_diff_%03d.rds", replicate_id))
saveRDS(sum_quant_no_ace_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/quant_no_ace_diff_%03d.rds", replicate_id))
saveRDS(sum_quant_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/quant_sample_diff_%03d.rds", replicate_id))
saveRDS(sum_quant_strict_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/quant_strict_diff_%03d.rds", replicate_id))
saveRDS(sum_quant_post_ord_point_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/quant_point_postord_diff_%03d.rds", replicate_id))
saveRDS(sum_quant_post_ord_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/quant_sample_postord_diff_%03d.rds", replicate_id))


# Apply sum, quantiles to rescaled variables
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

rescaled_sum_quant_point_postord <- lapply(rescaled_point_postord, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rescaled_sum_quant_sample_postord <- lapply(rescaled_sample_postord, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

# Apply sum, quantiles to remove.axes variables
rm_sum_quant_rel <- lapply(rm_rel, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rm_sum_quant_sample <- lapply(rm_sample, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rm_sum_quant_strict <- lapply(rm_strict, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rm_sum_quant_no_ace <- lapply(rm_no_ace, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rm_sum_quant_true <- lapply(rm_true, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rm_sum_quant_point_postord <- lapply(rm_point_postord, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

rm_sum_quant_sample_postord <- lapply(rm_sample_postord, lapply, lapply, function(rep){
    dispRity(rep, metric = c(sum, quantiles))$disparity
})

# Save rescaled results
saveRDS(rescaled_sum_quant_rel, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled_quant_rel_%03d.rds", replicate_id))
saveRDS(rescaled_sum_quant_sample, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled_quant_sample_%03d.rds", replicate_id))
saveRDS(rescaled_sum_quant_strict, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled_quant_strict_%03d.rds", replicate_id))
saveRDS(rescaled_sum_quant_no_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled_quant_no_ace_%03d.rds", replicate_id))
saveRDS(rescaled_sum_quant_true, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled_quant_true_%03d.rds", replicate_id))
saveRDS(rescaled_sum_quant_point_postord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled_quant_point_postord_%03d.rds", replicate_id))
saveRDS(rescaled_sum_quant_sample_postord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled_quant_sample_postord_%03d.rds", replicate_id))

# Save remove.axes results
saveRDS(rm_sum_quant_rel, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rm_quant_rel_%03d.rds", replicate_id))
saveRDS(rm_sum_quant_sample, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rm_quant_sample_%03d.rds", replicate_id))
saveRDS(rm_sum_quant_strict, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rm_quant_strict_%03d.rds", replicate_id))
saveRDS(rm_sum_quant_no_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rm_quant_no_ace_%03d.rds", replicate_id))
saveRDS(rm_sum_quant_true, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rm_quant_true_%03d.rds", replicate_id))
saveRDS(rm_sum_quant_point_postord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rm_quant_point_postord_%03d.rds", replicate_id))
saveRDS(rm_sum_quant_sample_postord, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rm_quant_sample_postord_%03d.rds", replicate_id))



# Disparity differences for rescaled variables
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
        (disp.diff(sample, rate_true))
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rescaled_sum_quant_sample, rescaled_sum_quant_true)

rescaled_sum_quant_strict_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, rescaled_sum_quant_strict, rescaled_sum_quant_true)

rescaled_sum_quant_post_ord_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, rescaled_sum_quant_point_postord, rescaled_sum_quant_true)

rescaled_sum_quant_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        (disp.diff(sample, rate_true))
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rescaled_sum_quant_sample_postord, rescaled_sum_quant_true)

# Disparity differences for remove.axes variables
rm_sum_quant_rel_diff <- Map(function(rate_rel, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_rel)
}, rm_sum_quant_rel, rm_sum_quant_true)

rm_sum_quant_no_ace_diff <- Map(function(rate_no_ace, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_no_ace)
}, rm_sum_quant_no_ace, rm_sum_quant_true)

rm_sum_quant_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        (disp.diff(sample, rate_true))
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rm_sum_quant_sample, rm_sum_quant_true)

rm_sum_quant_strict_diff <- Map(function(rate_strict, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_strict)
}, rm_sum_quant_strict, rm_sum_quant_true)

rm_sum_quant_post_ord_point_diff <- Map(function(rate_point, rate_true) {
    Map(function(fossil_data) {
      disp.diff(fossil_data, rate_true)
    }, rate_point)
}, rm_sum_quant_point_postord, rm_sum_quant_true)

rm_sum_quant_post_ord_sample_diff <- Map(function(rate_sample, rate_true) {
    Map(function(fossil_data) {
      diffs <- Map(function(sample) {
        (disp.diff(sample, rate_true))
      }, fossil_data)
      return(unlist(diffs))
    }, rate_sample)
}, rm_sum_quant_sample_postord, rm_sum_quant_true)

# Save rescaled difference results
saveRDS(rescaled_sum_quant_rel_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled_quant_rel_diff_%03d.rds", replicate_id))
saveRDS(rescaled_sum_quant_no_ace_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled_quant_no_ace_diff_%03d.rds", replicate_id))
saveRDS(rescaled_sum_quant_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled_quant_sample_diff_%03d.rds", replicate_id))
saveRDS(rescaled_sum_quant_strict_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled_quant_strict_diff_%03d.rds", replicate_id))
saveRDS(rescaled_sum_quant_post_ord_point_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled_quant_point_postord_diff_%03d.rds", replicate_id))
saveRDS(rescaled_sum_quant_post_ord_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rescaled_quant_sample_postord_diff_%03d.rds", replicate_id))

# Save remove.axes difference results
saveRDS(rm_sum_quant_rel_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rm_quant_rel_diff_%03d.rds", replicate_id))
saveRDS(rm_sum_quant_no_ace_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rm_quant_no_ace_diff_%03d.rds", replicate_id))
saveRDS(rm_sum_quant_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rm_quant_sample_diff_%03d.rds", replicate_id))
saveRDS(rm_sum_quant_strict_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rm_quant_strict_diff_%03d.rds", replicate_id))
saveRDS(rm_sum_quant_post_ord_point_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rm_quant_point_postord_diff_%03d.rds", replicate_id))
saveRDS(rm_sum_quant_post_ord_sample_diff, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/discrete/out/disparity/rm_quant_sample_postord_diff_%03d.rds", replicate_id))