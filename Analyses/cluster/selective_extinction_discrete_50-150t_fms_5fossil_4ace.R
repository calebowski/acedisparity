args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])

################################################################################
#                                                                              #
#                           Random Extinction SIMULATION                       #
#                                                                              #
################################################################################

library(treats)
source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")
source("/users/bip24cns/acedisparity/randomExtinction/scripts/find.extinction.time.R")


cat("Starting replicate", replicate_id, "\n")

set.seed(100 + replicate_id) # set seed to change for each replicate

bd_params <- make.bd.params(speciation = 0.1, extinction = 0.07)

stop_rule <- list(max.living = sample(x = c(50, 100, 150), size = 1)) # different tree sizes

slow_multi_transitions <- matrix(c(
    0.99, 0.025, 0.025,
    0.025, 0.99, 0.025,
    0.025, 0.025, 0.99
), nrow = 3, byrow = TRUE)


med_multi_transitions <- matrix(c(
  0.9, 0.05, 0.05,
  0.05, 0.9, 0.05,
  0.05, 0.05, 0.9
), nrow = 3, byrow = TRUE)


fast_multi_transitions <- matrix(c(
    0.1, 0.45, 0.45, 
    0.45, 0.1, 0.45,
    0.45, 0.45, 0.1 
), nrow = 3, byrow = TRUE)


slow_cor <- correlate.traits(correlation = 0.3, n.traits = 100, transition.matrix = slow_multi_transitions, n.correlated.traits = 50) ## high phylogenetic signal in correlation

med_cor <- correlate.traits(correlation = 0.25, n.traits = 100, transition.matrix = med_multi_transitions, n.correlated.traits = 50) 

fast_cor <- correlate.traits(correlation = 0.25, n.traits = 100, transition.matrix = fast_multi_transitions, n.correlated.traits = 50) ## low phylogenetic signal in the correlation

multi_traits  <- list(slow = slow_cor, med = med_cor, fast = fast_cor)

extinction_intensity <- runif(n = 1, min = 0.85, max = 0.95)
selective_extinction <- make.events(
                      target       = "taxa",
                      condition    = taxa.condition(stop_rule[[1]][[1]] - 10),
                      modification = trait.extinction(x = 1, trait = c(seq(from =  1, to = 50, by = 1)),
                                                      condition = `<`, severity = extinction_intensity, threshold = 0.3))



selective <- lapply(multi_traits, function(rate) {
  treats(traits = rate, stop.rule = stop_rule, bd.params = bd_params, null.error = 10000, events = selective_extinction)
})

selective <- lapply(selective, drop.singles)

matrices <- lapply(selective, function(x){
  x$data
})
saveRDS(matrices, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/matrices/selec_ext_matrices_%03d.rds", replicate_id))

trees <- lapply(selective, function(x){
  x$tree
})


trees <- lapply(trees, fix.zero.branches)
trees <- lapply(trees,  set.root.time)

saveRDS(matrices, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/trees/selec_ext_trees_%03d.rds", replicate_id))

source("/users/bip24cns/acedisparity/discrete/scripts/fossil.pres.R")



living <- Map(function(tree, matrix){ 
  remove.fossil(tree, matrix, type = "discrete")}, trees, matrices)


fossilised_high <- Map(function(tree, matrix){ 
  fossil.pres(tree, matrix, preservation = 0.5, type = "discrete")}, trees, matrices)

fossilised_med <- Map(function(tree, matrix){ 
  fossil.pres(tree, matrix, preservation = 0.15, type = "discrete")}, trees, matrices)

fossilised_low <- Map(function(tree, matrix){ 
  fossil.pres(tree, matrix, preservation = 0.05, type = "discrete")}, trees, matrices)

all_fossil <- Map(function(tree, matrix){ 
  fossil.pres(tree, matrix, preservation = 1.0, type = "discrete")}, trees, matrices)



fossil_matrices <- lapply(names(matrices), function(level) {
    list(
      all = all_fossil[[level]],
      fossil_high = fossilised_high[[level]],
      fossil_med = fossilised_med[[level]],
      fossil_low = fossilised_low[[level]],
      living = living[[level]]
    )
})


# Assign names to the outer list
names(fossil_matrices) <- names(matrices)

fossil_trees <- lapply(fossil_matrices, lapply,  function(level){
  tree <- level$tree
})

saveRDS(fossil_trees, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/trees/selec_ext_fossil_tree_%03d.rds", replicate_id))


cat("Fossil matrices created\n")
saveRDS(fossil_matrices, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/matrices/selec_ext_fossil_matrices_%03d.rds", replicate_id))

anc.states <- function(x) {
  # Run multi.ace for each tree
  anc_states <- multi.ace(data = x$matrix, 
                          tree = x$tree, 
                          models = "SYM", 
                          output = "multi.ace"
                          )
return(anc_states)}

fossil_anc <- lapply(fossil_matrices, lapply, anc.states)

sample_fossil_anc <- lapply(fossil_anc, lapply,  multi.ace, sample = 5)

strict_fossil_anc <- lapply(fossil_anc, lapply, multi.ace, threshold = FALSE, output = "combined.matrix", verbose = TRUE)

relative_fossil_anc <- lapply(fossil_anc, lapply,  multi.ace, output = "combined.matrix", verbose = TRUE)



cat("Ancestral state estimation completed...\n")

saveRDS(fossil_anc, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/anc/selec_ext_discrete_anc_%03d.rds", replicate_id))
saveRDS(sample_fossil_anc, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/anc/selec_ext_sample_anc_%03d.rds", replicate_id))
saveRDS(relative_fossil_anc, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/anc/selec_ext_rel_anc_%03d.rds", replicate_id))
saveRDS(strict_fossil_anc, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/anc/selec_ext_strict_anc_%03d.rds", replicate_id))


distances_true <- lapply(matrices, function(x){
  matrix <-  apply(x, c(1,2), as.character)
  dist <- char.diff(matrix, method = "mord", by.col = FALSE)
  return(dist)
})

distances_no_ace <- lapply(fossil_matrices, lapply, function(x){
  matrix <- x$matrix
    dist <- char.diff(matrix, method = "mord", by.col = FALSE)
    return(dist)
})

distances_strict <- lapply(strict_fossil_anc, lapply,  function(matrices) {
  dist <- char.diff(matrices, method = "mord", by.col = FALSE)
  return(dist)
})

distances_rel <- lapply(relative_fossil_anc, lapply, function(matrices) {
  dist <- char.diff(matrices, method = "mord", by.col = FALSE)
  return(dist)
})

distances_sample <- lapply(sample_fossil_anc, lapply, lapply, function(matrices) {
  dist <- char.diff(matrices, method = "mord", by.col = FALSE)
  return(dist)
})


cat("Distances completed...\n")

ord_no_ace <- lapply(distances_no_ace, lapply,  function(matrix){
    ord <-  (cmdscale(matrix, k = ncol(matrix) - 2, add = TRUE))$points
    return(ord)
})

ord_true <- lapply(distances_true, function(matrix){
    ord <-  (cmdscale(matrix, k = ncol(matrix) - 2, add = TRUE))$points
    return(ord)
})


ord_rel <- lapply(distances_rel, lapply, function(matrix) {
 cmdscale(matrix, k = ncol(matrix) - 2, add = TRUE)$points
})


ord_strict <- lapply(distances_strict, lapply, function(matrix) {
 cmdscale(matrix, k = ncol(matrix) - 2, add = TRUE)$points
})

ord_sample <- lapply(distances_sample, lapply, lapply, function(matrix) {
 cmdscale(matrix, k = ncol(matrix) - 2, add = TRUE)$points
})

cat("Ordinations completed...\n")


saveRDS(ord_no_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/ord/selec_ext_ord_no_ace_%03d.rds", replicate_id))
saveRDS(ord_true, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/ord/selec_ext_ord_true_%03d.rds", replicate_id))
saveRDS(ord_rel, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/ord/selec_ext_ord_rel_%03d.rds", replicate_id))
saveRDS(ord_strict, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/ord/selec_ext_ord_strict_%03d.rds", replicate_id))
saveRDS(ord_sample, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/ord/selec_ext_ord_sample_%03d.rds", replicate_id))



tip_ages <- lapply(trees,  tip.ages)
extinction_times <- lapply(tip_ages,  find.extinction.time)

time_slices <- lapply(extinction_times, function(time){
  intervals <- c(time + 20, time - 0.01,  (time - 0.01) - 20)
  return(intervals)
})


chrono_subsets_strict <- lapply(names(ord_strict), function(rate) {
  fossil_result <- lapply(names(ord_strict[[rate]]), function(fossil){ 
    tryCatch({
      chrono <- chrono.subsets(ord_strict[[rate]][[fossil]],
      tree = fossil_trees[[rate]][[fossil]],
      method = "d",
      time = time_slices[[rate]],
      inc.nodes = TRUE)
      return(chrono)
    }, warning = function(w) {
      rm(w)
    })
  })
  names(fossil_result) <- names(ord_strict[[rate]])
  return(fossil_result)
})
names(chrono_subsets_strict) <- names(ord_strict)


chrono_subsets_rel <- lapply(names(ord_rel), function(rate) {
  fossil_result <- lapply(names(ord_rel[[rate]]), function(fossil){ 
    tryCatch({
      chrono <- chrono.subsets(ord_rel[[rate]][[fossil]],
      tree = fossil_trees[[rate]][[fossil]],
      method = "d",
      time = time_slices[[rate]],
      inc.nodes = TRUE)
      return(chrono)
    }, warning = function(w) {
      rm(w)
    })
  })
  names(fossil_result) <- names(ord_rel[[rate]])
  return(fossil_result)
})
names(chrono_subsets_rel) <- names(ord_rel)

chrono_subsets_no_ace <- lapply(names(ord_no_ace), function(rate) {
  fossil_result <- lapply(names(ord_no_ace[[rate]]), function(fossil){ 
    tryCatch({
      chrono <- chrono.subsets(ord_no_ace[[rate]][[fossil]],
      tree = fossil_trees[[rate]][[fossil]],
      method = "d",
      time = time_slices[[rate]]) # do not include nodes
      return(chrono)
    }, warning = function(w) {
      rm(w)
    })
  })
  names(fossil_result) <- names(ord_no_ace[[rate]])
  return(fossil_result)
})
names(chrono_subsets_no_ace) <- names(ord_no_ace)


chrono_subsets_true <- lapply(names(ord_true), function(rate) {
    tryCatch({
      chrono <- chrono.subsets(ord_true[[rate]],
      tree = trees[[rate]], # use tree rather than fossil trees
      method = "d",
      time = time_slices[[rate]],
      inc.nodes = TRUE)
      return(chrono)
    }, warning = function(w) {
      rm(w)
    })
})
names(chrono_subsets_true) <- names(ord_true)


chrono_subsets_sample <- lapply(names(ord_sample), function(rate) {
  fossil_result <- lapply(names(ord_sample[[rate]]), function(fossil) {
    rep_result <- lapply(seq_along(ord_sample[[rate]][[fossil]]), function(i) {
      tryCatch({
        chrono <- chrono.subsets(ord_sample[[rate]][[fossil]][[i]],
        tree = fossil_trees[[rate]][[fossil]],
        method = "d",
        time = time_slices[[rate]],
        inc.nodes = TRUE)
        return(chrono)
      }, warning = function(w) {
        rm(w)
      })
    })
    return(rep_result)
  })
  names(fossil_result) <- names(ord_sample[[rate]])
  return(fossil_result)
})
names(chrono_subsets_sample) <- names(ord_sample)

################################################################################
#                                                                              #
#                           SUM OF VARIANCES                                   #
#                                                                              #
################################################################################


true_sum_var_change <- lapply(chrono_subsets_true, function(mat) {
  disp <- dispRity(mat, metric = c(sum, variances))
  values <- get.disparity(disp)
  change <- (values[[2]] - values[[1]]) / values[[1]]
  return(change)
})

strict_sum_var_change <- lapply(chrono_subsets_strict, lapply,  function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(sum, variances))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)
    }, warning = function(w) {
      NULL
    })
})


rel_sum_var_change <- lapply(chrono_subsets_rel, lapply, function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(sum, variances))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)
    }, warning = function(w) {
      NULL
    })
})

sample_sum_var_change <- lapply(chrono_subsets_sample, lapply, lapply,  function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(sum, variances))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)
    }, warning = function(w) {
      NULL
    })
})

no_ace_sum_var_change <- lapply(chrono_subsets_no_ace, lapply,  function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(sum, variances))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)}, warning = function(w) {
      NULL
    })
})


saveRDS(strict_sum_var_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_sum_var_strict_%03d.rds", replicate_id))
saveRDS(rel_sum_var_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_sum_var_rel_%03d.rds", replicate_id))
saveRDS(sample_sum_var_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_sum_var_sample_%03d.rds", replicate_id))
saveRDS(no_ace_sum_var_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_sum_var_no_ace_%03d.rds", replicate_id))
saveRDS(true_sum_var_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_sum_var_true_%03d.rds", replicate_id))


################################################################################
#                                                                              #
#                           AV NEAREST NEIGHBOUR                               #
#                                                                              #
################################################################################

true_neighbours_change <- lapply(chrono_subsets_true, function(mat) {
  disp <- dispRity(mat, metric = c(mean, neighbours))
  values <- get.disparity(disp)
  change <- (values[[2]] - values[[1]]) / values[[1]]
  return(change)
})

strict_neighbours_change <- lapply(chrono_subsets_strict, lapply,  function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(mean, neighbours))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)
    }, warning = function(w) {
      NULL
    })
})


rel_neighbours_change <- lapply(chrono_subsets_rel, lapply, function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(mean, neighbours))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)}, warning = function(w) {
      NULL
    })
})

sample_neighbours_change <- lapply(chrono_subsets_sample, lapply, lapply,  function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(mean, neighbours))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)}, warning = function(w) {
      NULL
    })
})

no_ace_neighbours_change <- lapply(chrono_subsets_no_ace, lapply,  function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(mean, neighbours))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)}, warning = function(w) {
      NULL
    })
})


saveRDS(strict_neighbours_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_neighbours_strict_%03d.rds", replicate_id))
saveRDS(rel_neighbours_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_neighbours_rel_%03d.rds", replicate_id))
saveRDS(sample_neighbours_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_neighbours_sample_%03d.rds", replicate_id))
saveRDS(no_ace_neighbours_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_neighbours_no_ace_%03d.rds", replicate_id))
saveRDS(true_neighbours_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_neighbours_true_%03d.rds", replicate_id))


################################################################################
#                                                                              #
#                           AV DISPLACEMENT                                    #
#                                                                              #
################################################################################

true_displacement_change <- lapply(chrono_subsets_true, function(mat) {
  disp <- dispRity(mat, metric = c(mean, displacements))
  values <- get.disparity(disp)
  change <- (values[[2]] - values[[1]]) / values[[1]]
  return(change)
})

strict_displacement_change <- lapply(chrono_subsets_strict, lapply,  function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(mean, displacements))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)
    }, warning = function(w) {
      NULL
    })
})


rel_displacement_change <- lapply(chrono_subsets_rel, lapply, function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(mean, displacements))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)}, warning = function(w) {
      NULL
    })
})

sample_displacement_change <- lapply(chrono_subsets_sample, lapply, lapply,  function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(mean, displacements))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)}, warning = function(w) {
      NULL
    })
})

no_ace_displacement_change <- lapply(chrono_subsets_no_ace, lapply,  function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(mean, displacements))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)}, warning = function(w) {
      NULL
    })
})


saveRDS(strict_displacement_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_displacement_strict_%03d.rds", replicate_id))
saveRDS(rel_displacement_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_displacement_rel_%03d.rds", replicate_id))
saveRDS(sample_displacement_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_displacement_sample_%03d.rds", replicate_id))
saveRDS(no_ace_displacement_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_displacement_no_ace_%03d.rds", replicate_id))
saveRDS(true_displacement_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/disp/selec_ext_change_displacement_true_%03d.rds", replicate_id))



metadata <- list(
  replicate_id = replicate_id,
  tree_size = stop_rule$max.living,
  extinction_intensity = extinction_intensity,  # If you want to track this too
  seed = 100 + replicate_id
)

saveRDS(metadata, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/selectiveExtinction/out/metadata/metadata_%03d.rds", replicate_id))