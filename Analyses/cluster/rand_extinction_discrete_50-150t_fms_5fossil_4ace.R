args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])

################################################################################
#                                                                              #
#                           Random Extinction SIMULATION                       #
#                                                                              #
################################################################################


# library(remotes)

# install_github("TGuillerme/treats", ref = "selective.extinction.alteration")

library(treats)
# library(parallel)
source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")
source("/users/bip24cns/acedisparity/randomExtinction/scripts/find.extinction.time.R")


cat("Starting replicate", replicate_id, "\n")

set.seed(100 + replicate_id) # set seed to change for each replicate

bd_params <- make.bd.params(speciation = 0.1, extinction = 0.07)

stop_rule <- list(max.living = sample(x = c(50, 100, 150), size = 1)) # different tree sizes

extinction_intensity <- runif(n = 1, min = 0.75, max = 0.9)


random_extinction <- make.events(
                      target       = "taxa",
                      condition    = taxa.condition(stop_rule[[1]][[1]] - 10),
                      modification = random.extinction(extinction_intensity)
)

tree <- treats(stop.rule = stop_rule, bd.params = bd_params, null.error = 100, events = random_extinction)


tree <- drop.singles(tree)
tree <- fix.zero.branches(tree)
tree <- set.root.time(tree)


## fix zero branches
## set root time

## fix single nodes.

write.tree(tree, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/trees/rand_ext_tree_%03d.tre", replicate_id))

cat("Mapping traits...\n")

slow_binary_transitions <- matrix(c(
    0.99, 0.01,
    0.01, 0.99
), nrow = 2, byrow = TRUE)

slow_multi_transitions <- matrix(c(
    0.99, 0.025, 0.025,
    0.025, 0.99, 0.025,
    0.025, 0.025, 0.99
), nrow = 3, byrow = TRUE)

slow_binary <- treats::make.traits(process = discrete.process, n = 85, process.args = list(transitions = slow_binary_transitions))
slow_multi <- treats::make.traits(process = discrete.process, n = 15, process.args = list(transitions = slow_multi_transitions))


med_binary_transitions <- matrix(c(
  0.90, 0.1,
  0.1, 0.9
), nrow = 2, byrow = TRUE)


med_multi_transitions <- matrix(c(
  0.9, 0.05, 0.05,
  0.05, 0.9, 0.05,
  0.05, 0.05, 0.9
), nrow = 3, byrow = TRUE)

med_binary <- treats::make.traits(process = discrete.process, n = 85, process.args = list(transitions = med_binary_transitions))
med_multi <- treats::make.traits(process = discrete.process, n = 15, process.args = list(transitions = med_multi_transitions))

fast_binary_transitions <- matrix(c(
    0.1, 0.9,
    0.9, 0.1
), nrow = 2, byrow = TRUE)

fast_multi_transitions <- matrix(c(
    0.1, 0.45, 0.45, 
    0.45, 0.1, 0.45,
    0.45, 0.45, 0.1 
), nrow = 3, byrow = TRUE)

fast_binary <- treats::make.traits(process = discrete.process, n = 85, process.args = list(transitions = fast_binary_transitions))
fast_multi <- treats::make.traits(process = discrete.process, n = 15, process.args = list(transitions = fast_multi_transitions))



## map the traits

# Create list of trait objects
trait_sets <- list(
  slow = list(binary = slow_binary, multi = slow_multi),
  med = list(binary = med_binary, multi = med_multi),
  fast = list(binary = fast_binary, multi = fast_multi)
)

# Map all traits and combine
matrices <- lapply(trait_sets, function(traits) {
  mapped_binary <- (map.traits(traits$binary, tree))$data
  mapped_multi <- (map.traits(traits$multi, tree))$data
  cbind(mapped_binary, mapped_multi)
})

saveRDS(matrices, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/matrices/rand_ext_matrices_%03d.rds", replicate_id))

cat("Trait matrices saved for replicate", replicate_id, "\n")

source("/users/bip24cns/acedisparity/discrete/scripts/fossil.pres.R")

living <- lapply(matrices, remove.fossil, trees = tree, type = "discrete")
fossilised_high <- lapply(matrices, fossil.pres, trees = tree, preservation = 0.5, type = "discrete")
all_fossil <- lapply(matrices, fossil.pres, trees = tree, preservation = 1.0, type = "discrete")
fossilised_med <- lapply(matrices, fossil.pres, trees = tree, preservation = 0.15, type = "discrete")
fossilised_low <- lapply(matrices, fossil.pres, trees = tree, preservation = 0.05, type = "discrete")



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

saveRDS(fossil_trees, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/trees/rand_ext_fossil_tree_%03d.rds", replicate_id))


cat("Fossil matrices created\n")
saveRDS(fossil_matrices, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/matrices/rand_ext_fossil_matrices_%03d.rds", replicate_id))


# n_cores  <- 5

anc.states <- function(x) {
  # Run multi.ace for each tree
  anc_states <- multi.ace(data = x$matrix, 
                          tree = x$tree, 
                          models = "SYM", 
                          output = "multi.ace"
                          )
return(anc_states)}

fossil_anc <- lapply(fossil_matrices, lapply, anc.states)

sample_fossil_anc <- lapply(fossil_anc, lapply,  multi.ace, sample = 100)

strict_fossil_anc <- lapply(fossil_anc, lapply, multi.ace, threshold = FALSE, output = "combined.matrix", verbose = TRUE)


relative_fossil_anc <- lapply(fossil_anc, lapply,  multi.ace, output = "combined.matrix", verbose = TRUE)

cat("Ancestral state estimation completed...\n")


saveRDS(fossil_anc, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/anc/rand_ext_discrete_anc_%03d.rds", replicate_id))
saveRDS(sample_fossil_anc, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/anc/rand_ext_sample_anc_%03d.rds", replicate_id))
saveRDS(relative_fossil_anc, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/anc/rand_ext_rel_anc_%03d.rds", replicate_id))
saveRDS(strict_fossil_anc, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/anc/rand_ext_strict_anc_%03d.rds", replicate_id))


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


post_ord_ace <- Map(function(rate_matrix, rate_tree) {
  Map(function(fossil_matrix, fossil_tree) {
    # clean <- clean.data(fossil_matrix, fossil_tree)
    # tree <- clean$tree
    # matrix <- clean$data
    # return(list(matrix = matrix, tree = tree))
    multi.ace(fossil_matrix, fossil_tree, models = "ML", output = "multi.ace")
  }, rate_matrix, rate_tree)
}, ord_no_ace, fossil_trees)

trait_normal = list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))
point_post_ord_ace <- lapply(post_ord_ace, lapply,  multi.ace, output = "combined.matrix")
sample_post_ord_ace <- lapply(post_ord_ace, lapply,  multi.ace, sample = 100, sample.fun = trait_normal, output = "combined.matrix")

saveRDS(point_post_ord_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/ord/post_ord_ace_point_%03d.rds", replicate_id))
saveRDS(sample_post_ord_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/ord/post_ord_ace_sample_%03d.rds", replicate_id))


cat("Ordinations completed...\n")

saveRDS(ord_no_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/ord/rand_ext_ord_no_ace_%03d.rds", replicate_id))
saveRDS(ord_true, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/ord/rand_ext_ord_true_%03d.rds", replicate_id))
saveRDS(ord_rel, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/ord/rand_ext_ord_rel_%03d.rds", replicate_id))
saveRDS(ord_strict, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/ord/rand_ext_ord_strict_%03d.rds", replicate_id))
saveRDS(ord_sample, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/ord/rand_ext_ord_sample_%03d.rds", replicate_id))

tip_ages <- tip.ages(tree)
extinction_time <- find.extinction.time(tip_ages)

time_slices <- c(extinction_time + 20, extinction_time - 0.01,  (extinction_time - 0.01) - 20)

chrono_subsets_strict <- lapply(names(ord_strict), function(rate) {
  fossil_result <- lapply(names(ord_strict[[rate]]), function(fossil){ 
    tryCatch({
      chrono <- chrono.subsets(ord_strict[[rate]][[fossil]],
      tree = fossil_trees[[rate]][[fossil]],
      method = "d",
      time = time_slices,
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
      time = time_slices,
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
      time = time_slices) # do not include nodes
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
      tree = tree, # use tree rather than fossil trees
      method = "d",
      time = time_slices,
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
        time = time_slices,
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


chrono_subsets_point_post_ord <- lapply(names(point_post_ord_ace), function(rate) {
  fossil_result <- lapply(names(point_post_ord_ace[[rate]]), function(fossil){ 
    tryCatch({
      chrono <- chrono.subsets(point_post_ord_ace[[rate]][[fossil]],
      tree = fossil_trees[[rate]][[fossil]],
      method = "d",
      time = time_slices,
      inc.nodes = TRUE)
      return(chrono)
    }, warning = function(w) {
      rm(w)
    })
  })
  names(fossil_result) <- names(point_post_ord_ace[[rate]])
  return(fossil_result)
})
names(chrono_subsets_point_post_ord) <- names(point_post_ord_ace)

chrono_subsets_sample_post_ord <- lapply(names(sample_post_ord_ace), function(rate) {
  fossil_result <- lapply(names(sample_post_ord_ace[[rate]]), function(fossil) {
    rep_result <- lapply(seq_along(sample_post_ord_ace[[rate]][[fossil]]), function(i) {
      tryCatch({
        chrono <- chrono.subsets(sample_post_ord_ace[[rate]][[fossil]][[i]],
        tree = fossil_trees[[rate]][[fossil]],
        method = "d",
        time = time_slices,
        inc.nodes = TRUE)
        return(chrono)
      }, warning = function(w) {
        rm(w)
      })
    })
    return(rep_result)
  })
  names(fossil_result) <- names(sample_post_ord_ace[[rate]])
  return(fossil_result)
})
names(chrono_subsets_sample_post_ord) <- names(sample_post_ord_ace)

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

point_post_ord_sum_var_change <- lapply(chrono_subsets_point_post_ord, lapply,  function(chrono) {
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

sample_post_ord_sum_var_change <- lapply(chrono_subsets_sample_post_ord, lapply, lapply,  function(chrono) {
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
# saveRDS(sum_var_true, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_sum_var_true_%03d.rds", replicate_id))
# saveRDS(sum_var_strict, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_sum_var_strict_%03d.rds", replicate_id))
# saveRDS(sum_var_rel, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_sum_var_rel_%03d.rds", replicate_id))
# saveRDS(sum_var_sample, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_sum_var_sample_%03d.rds", replicate_id))
# saveRDS(sum_var_no_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_sum_var_no_ace_%03d.rds", replicate_id))


# sum_var_list <- list(strict = sum_var_strict, rel = sum_var_rel, sample = sum_var_sample, no_ace = sum_var_no_ace)
### get disparity

# strict_sum_var_change <- lapply(sum_var_strict, lapply,  function(disp) {
#   if(is.null(disp)) {
#     return(NULL)
#   } else{
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
#   }
# })

# rel_sum_var_change <- lapply(sum_var_rel, lapply,  function(disp) {
#   if(is.null(disp)) {
#     return(NULL)
#   } else{
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
#   }
# })

# sample_sum_var_change <- lapply(sum_var_sample, lapply, lapply, function(disp) {
#   if(is.null(disp)) {
#     return(NULL)
#   } else{
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
#   }
# })

# no_ace_sum_var_change <- lapply(sum_var_no_ace, lapply,  function(disp) {
#   if(is.null(disp)) {
#     return(NULL)
#   } else{
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
#   }
# })

# true_sum_var_change <- lapply(sum_var_true,  function(disp){
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
# })

saveRDS(strict_sum_var_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_sum_var_strict_%03d.rds", replicate_id))
saveRDS(rel_sum_var_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_sum_var_rel_%03d.rds", replicate_id))
saveRDS(sample_sum_var_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_sum_var_sample_%03d.rds", replicate_id))
saveRDS(no_ace_sum_var_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_sum_var_no_ace_%03d.rds", replicate_id))
saveRDS(true_sum_var_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_sum_var_true_%03d.rds", replicate_id))
saveRDS(sample_post_ord_sum_var_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_sum_var_sample_postord_%03d.rds", replicate_id))
saveRDS(point_post_ord_sum_var_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_sum_var_point_postord_%03d.rds", replicate_id))



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

point_post_ord_neighbours_change <- lapply(chrono_subsets_point_post_ord, lapply,  function(chrono) {
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

sample_post_ord_neighbours_change <- lapply(chrono_subsets_sample_post_ord, lapply, lapply,  function(chrono) {
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
# saveRDS(neighbours_true, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_neighbours_true_%03d.rds", replicate_id))
# saveRDS(neighbours_strict, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_neighbours_strict_%03d.rds", replicate_id))
# saveRDS(neighbours_rel, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_neighbours_rel_%03d.rds", replicate_id))
# saveRDS(neighbours_sample, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_neighbours_sample_%03d.rds", replicate_id))
# saveRDS(neighbours_no_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_neighbours_no_ace_%03d.rds", replicate_id))


# sum_var_list <- list(strict = sum_var_strict, rel = sum_var_rel, sample = sum_var_sample, no_ace = sum_var_no_ace)
### get disparity

# strict_neighbours_change <- lapply(neighbours_strict, lapply,  function(disp) {
#   if(is.null(disp)) {
#     return(NULL)
#   } else{
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
#   }
# })

# rel_neighbours_change <- lapply(neighbours_rel, lapply,  function(disp) {
#   if(is.null(disp)) {
#     return(NULL)
#   } else{
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
#   }
# })

# sample_neighbours_change <- lapply(neighbours_sample, lapply, lapply, function(disp) {
#   if(is.null(disp)) {
#     return(NULL)
#   } else{
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
#   }
# })

# no_ace_neighbours_change <- lapply(neighbours_no_ace, lapply,  function(disp) {
#   if(is.null(disp)) {
#     return(NULL)
#   } else{
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
#   }
# })

# true_neighbours_change <- lapply(neighbours_true,  function(disp){
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
# })

saveRDS(strict_neighbours_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_neighbours_strict_%03d.rds", replicate_id))
saveRDS(rel_neighbours_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_neighbours_rel_%03d.rds", replicate_id))
saveRDS(sample_neighbours_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_neighbours_sample_%03d.rds", replicate_id))
saveRDS(no_ace_neighbours_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_neighbours_no_ace_%03d.rds", replicate_id))
saveRDS(true_neighbours_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_neighbours_true_%03d.rds", replicate_id))
saveRDS(point_post_ord_neighbours_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_neighbours_point_postord_%03d.rds", replicate_id))
saveRDS(sample_post_ord_neighbours_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_neighbours_sample_postord_%03d.rds", replicate_id))

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

point_post_ord_displacement_change <- lapply(chrono_subsets_point_post_ord, lapply,  function(chrono) {
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

sample_post_ord_displacement_change <- lapply(chrono_subsets_sample_post_ord, lapply, lapply,  function(chrono) {
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

# saveRDS(displacement_true, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_displacement_true_%03d.rds", replicate_id))
# saveRDS(displacement_strict, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_displacement_strict_%03d.rds", replicate_id))
# saveRDS(displacement_rel, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_displacement_rel_%03d.rds", replicate_id))
# saveRDS(displacement_sample, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_displacement_sample_%03d.rds", replicate_id))
# saveRDS(displacement_no_ace, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_displacement_no_ace_%03d.rds", replicate_id))


# sum_var_list <- list(strict = sum_var_strict, rel = sum_var_rel, sample = sum_var_sample, no_ace = sum_var_no_ace)
### get disparity

# strict_displacement_change <- lapply(displacement_strict, lapply,  function(disp) {
#   if(is.null(disp)) {
#     return(NULL)
#   } else{
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
#   }
# })

# rel_displacement_change <- lapply(displacement_rel, lapply,  function(disp) {
#   if(is.null(disp)) {
#     return(NULL)
#   } else{
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
#   }
# })

# sample_displacement_change <- lapply(displacement_sample, lapply, lapply, function(disp) {
#   if(is.null(disp)) {
#     return(NULL)
#   } else{
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
#   }
# })

# no_ace_displacement_change <- lapply(displacement_no_ace, lapply,  function(disp) {
#   if(is.null(disp)) {
#     return(NULL)
#   } else{
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
#   }
# })

# true_displacement_change <- lapply(displacement_true,  function(disp){
#   values <- get.disparity(disp)
#   change <- (values[[2]] - values[[1]]) / values[[1]]
# })

saveRDS(strict_displacement_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_displacement_strict_%03d.rds", replicate_id))
saveRDS(rel_displacement_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_displacement_rel_%03d.rds", replicate_id))
saveRDS(sample_displacement_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_displacement_sample_%03d.rds", replicate_id))
saveRDS(no_ace_displacement_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_displacement_no_ace_%03d.rds", replicate_id))
saveRDS(true_displacement_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_displacement_true_%03d.rds", replicate_id))
saveRDS(point_post_ord_displacement_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_displacement_point_postord_%03d.rds", replicate_id))
saveRDS(sample_post_ord_displacement_change, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/disp/rand_ext_change_displacement_sample_postord_%03d.rds", replicate_id))

metadata <- list(
  replicate_id = replicate_id,
  tree_size = stop_rule$max.living,
  extinction_intensity = extinction_intensity,  # If you want to track this too
  seed = 100 + replicate_id
)

saveRDS(metadata, sprintf("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/out/metadata/metadata_%03d.rds", replicate_id))