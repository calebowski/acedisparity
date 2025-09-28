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
library(dispRity)
# library(parallel)
source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")
source("/users/bip24cns/acedisparity/randomExtinction/scripts/find.extinction.time.R")


# Build the base path with job-specific directory
base_path <- "/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/100t/"
job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")



write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}

cat("Starting replicate", replicate_id, "\n")

set.seed(100 + replicate_id) # set seed to change for each replicate

bd_params <- make.bd.params(speciation = 1.0, extinction = 0.5)

stop_rule <- list(max.living = 100) # different tree sizes

extinction_intensity <- runif(n = 1, min = 0.75, max = 0.95)


random_extinction <- make.events(
                      target       = "taxa",
                      condition    = taxa.condition(stop_rule[[1]][[1]] - 10),
                      modification = random.extinction(extinction_intensity)
)

# tree <- treats(stop.rule = stop_rule, bd.params = bd_params, null.error = 100, events = random_extinction)

# Generate tree with feedback loop to prevent >300 tips
max_attempts <- 200  # Prevent infinite loops
attempt <- 1


repeat {
  tree <- treats(stop.rule = stop_rule, bd.params = bd_params, null.error = 100, events = random_extinction)
  n_tips <- length(tree$tip.label)
  
  if(n_tips <= 300) {
    cat("Tree found with", n_tips, "tips\n")
    break
  }
  
  if(attempt >= max_attempts) {
    cat("Warning: Reached maximum attempts, accepting tree with", n_tips, "tips\n")
    break
  }
  
  attempt <- attempt + 1
}
# gc()

# est <- crude.bd.est(tree, method = "estimate")
# ape_est <- birthdeath(tree)

tree <- drop.singles(tree)
tree <- fix.zero.branches(tree)
tree <- set.root.time(tree)

write.tree(tree, write.path("trees", "rand_ext_tree_%03d.tre"))

# fossils <- sim.fossils.poisson(rate = 0.25, tree)

metadata_df <- data.frame(
  replicate_id = replicate_id,
  tree_size = stop_rule$max.living,
  extinction_intensity = extinction_intensity,
  seed = 100 + replicate_id
)

write.csv(metadata_df, 
          write.path("metadata","metadata_%03d.csv"),
          row.names = FALSE)


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

saveRDS(matrices, write.path("matrices", "matrices_%03d.rds"))

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

saveRDS(fossil_trees, write.path("trees", "fossil_trees_%03d.rds"))


cat("Fossil matrices created\n")
saveRDS(fossil_matrices, write.path("matrices", "fossil_matrices_%03d.rds"))


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


saveRDS(fossil_anc, write.path("anc", "fossil_anc_%03d.rds"))
saveRDS(sample_fossil_anc, write.path("anc", "sample_anc_%03d.rds"))
saveRDS(relative_fossil_anc, write.path("anc", "rel_anc_%03d.rds"))
saveRDS(strict_fossil_anc, write.path("anc", "strict_anc_%03d.rds"))


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
ord_no_ace <- lapply(distances_no_ace, lapply, function(matrix) {
  # Check for NA values in the distance matrix
  if(any(is.na(matrix))) {
    cat("Warning: NA values found in distance matrix, skipping ordination\n")
    return(NULL)
  }
  tryCatch({
    ord <- (cmdscale(matrix, k = ncol(matrix) - 2, add = TRUE))$points
    
    # Check if ordination result has NAs
    if(any(is.na(ord))) {
      cat("Warning: NA values in ordination result\n")
      # Still return the result even with NAs
    }
    
    return(ord)
  }, error = function(e) {
    cat("Error in cmdscale:", e$message, "\n")
    return(NULL)
  })
})

ord_true <- lapply(distances_true, function(matrix){
  if(any(is.na(matrix))) {
    cat("Warning: NA values in true distance matrix\n")
    return(NULL)
  }
  
  tryCatch({
    ord <- (cmdscale(matrix, k = ncol(matrix) - 2, add = TRUE))$points
    if(any(is.na(ord))) {
      cat("Warning: NA values in true ordination result\n")
      # Still return the result even with NAs
    }
    return(ord)
  }, error = function(e) {
    cat("Error in ordination:", e$message, "\n")
    return(NULL)
  })
})

ord_rel <- lapply(distances_rel, lapply, function(matrix) {
  if(any(is.na(matrix))) {
    cat("Warning: NA values in distance matrix\n")
    return(NULL)
  }
  tryCatch({
    ord <- cmdscale(matrix, k = ncol(matrix) - 2, add = TRUE)$points
    if(any(is.na(ord))) {
      cat("Warning: NA values in ordination result\n")
      # Still return the result even with NAs
    }
    return(ord)
  }, error = function(e) {
    cat("Error in ordination:", e$message, "\n")
    return(NULL)
  })
})

ord_strict <- lapply(distances_strict, lapply, function(matrix) {
  if(any(is.na(matrix))) {
    cat("Warning: NA values in distance matrix\n")
    return(NULL)
  }
  tryCatch({
    ord <- cmdscale(matrix, k = ncol(matrix) - 2, add = TRUE)$points
    if(any(is.na(ord))) {
      cat("Warning: NA values in ordination result\n")
      # Still return the result even with NAs
    }
    return(ord)
  }, error = function(e) {
    cat("Error in ordination:", e$message, "\n")
    return(NULL)
  })
})

ord_sample <- lapply(distances_sample, lapply, lapply, function(matrix) {
  if(any(is.na(matrix))) {
    cat("Warning: NA values in distance matrix\n")
    return(NULL)
  }
  tryCatch({
    ord <- cmdscale(matrix, k = ncol(matrix) - 2, add = TRUE)$points
    if(any(is.na(ord))) {
      cat("Warning: NA values in ordination result\n")
      # Still return the result even with NAs
    }
    return(ord)
  }, error = function(e) {
    cat("Error in ordination:", e$message, "\n")
    return(NULL)
  })
})

cat("Ordinations completed...\n")

saveRDS(ord_no_ace, write.path("ord", "ord_no_ace_%03d.rds"))
saveRDS(ord_true, write.path("ord", "ord_true_%03d.rds"))
saveRDS(ord_rel, write.path("ord", "ord_rel_%03d.rds"))
saveRDS(ord_strict, write.path("ord", "ord_strict_%03d.rds"))
saveRDS(ord_sample, write.path("ord", "ord_sample_%03d.rds"))

cat("Ordinations saved...\n")

post_ord_ace <- Map(function(rate_matrix, rate_tree) {
  Map(function(fossil_matrix, fossil_tree) {
      multi.ace(fossil_matrix, fossil_tree, models = "BM", output = "multi.ace")
  }, rate_matrix, rate_tree)
}, ord_no_ace, fossil_trees)

trait_normal  <-  list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))

point_post_ord_ace <- lapply(post_ord_ace, lapply,  multi.ace, output = "combined.matrix")

sample_post_ord_ace <- lapply(post_ord_ace, lapply, function(x){
  multi.ace(x, sample = 100, sample.fun = trait_normal, output = "combined.matrix")
})

cat("Post ordination estimates completed...\n")


saveRDS(point_post_ord_ace, write.path("ord", "post_ord_point_%03d.rds"))
saveRDS(sample_post_ord_ace, write.path("ord", "post_ord_sample_%03d.rds"))


cat("Post ordination estimates saved...\n")

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

cat("chrono subsets completed...\n")


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

cat("Sum variances completed...\n")


saveRDS(strict_sum_var_change, write.path("disp", "strict_sum_var_change_%03d.rds"))
saveRDS(rel_sum_var_change, write.path("disp", "rel_sum_var_change_%03d.rds"))
saveRDS(sample_sum_var_change, write.path("disp", "sample_sum_var_change_%03d.rds"))
saveRDS(no_ace_sum_var_change, write.path("disp", "no_ace_sum_var_change_%03d.rds"))
saveRDS(true_sum_var_change, write.path("disp", "true_sum_var_change_%03d.rds"))
saveRDS(sample_post_ord_sum_var_change, write.path("disp", "sample_postord_sum_var_change_%03d.rds"))
saveRDS(point_post_ord_sum_var_change, write.path("disp", "point_postord_sum_var_change_%03d.rds"))

################################################################################
#                                                                              #
#                           SUM OF QUANTILES                                   #
#                                                                              #
################################################################################

true_sum_quant_change <- lapply(chrono_subsets_true, function(mat) {
  disp <- dispRity(mat, metric = c(sum, quantiles))
  values <- get.disparity(disp)
  change <- (values[[2]] - values[[1]]) / values[[1]]
  return(change)
})

strict_sum_quant_change <- lapply(chrono_subsets_strict, lapply,  function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(sum, quantiles))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)
    }, warning = function(w) {
      NULL
    })
})


rel_sum_quant_change <- lapply(chrono_subsets_rel, lapply, function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(sum, quantiles))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)
    }, warning = function(w) {
      NULL
    })
})

sample_sum_quant_change <- lapply(chrono_subsets_sample, lapply, lapply,  function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(sum, quantiles))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)
    }, warning = function(w) {
      NULL
    })
})

no_ace_sum_quant_change <- lapply(chrono_subsets_no_ace, lapply,  function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(sum, quantiles))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)}, warning = function(w) {
      NULL
    })
})

point_post_ord_sum_quant_change <- lapply(chrono_subsets_point_post_ord, lapply,  function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(sum, quantiles))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)}, warning = function(w) {
      NULL
    })
})

sample_post_ord_sum_quant_change <- lapply(chrono_subsets_sample_post_ord, lapply, lapply,  function(chrono) {
  if(is.null(chrono)) {
    return(NULL)
  }
  tryCatch({
    disp <- dispRity(chrono, metric = c(sum, quantiles))
    values <- get.disparity(disp)
    change <- (values[[2]] - values[[1]]) / values[[1]]
    return(change)
    }, warning = function(w) {
      NULL
    })
})

saveRDS(strict_sum_quant_change, write.path("disp", "strict_sum_quant_change_%03d.rds"))
saveRDS(rel_sum_quant_change, write.path("disp", "rel_sum_quant_change_%03d.rds"))
saveRDS(sample_sum_quant_change, write.path("disp", "sample_sum_quant_change_%03d.rds"))
saveRDS(no_ace_sum_quant_change, write.path("disp", "no_ace_sum_quant_change_%03d.rds"))
saveRDS(true_sum_quant_change, write.path("disp", "true_sum_quant_change_%03d.rds"))
saveRDS(sample_post_ord_sum_quant_change, write.path("disp", "sample_postord_sum_quant_change_%03d.rds"))
saveRDS(point_post_ord_sum_quant_change, write.path("disp", "point_postord_sum_quant_change_%03d.rds"))

cat("Sum quantiles completed...\n")


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
    }, error = function(e) {
      cat("Error produced with sample post ord neighbours:", e$message, "\n")
      NULL
    })
})

saveRDS(true_neighbours_change, write.path("disp", "true_neighbours_change_%03d.rds"))
saveRDS(strict_neighbours_change, write.path("disp", "strict_neighbours_change_%03d.rds"))
saveRDS(rel_neighbours_change, write.path("disp", "rel_neighbours_change_%03d.rds"))
saveRDS(sample_neighbours_change, write.path("disp", "sample_neighbours_change_%03d.rds"))
saveRDS(no_ace_neighbours_change, write.path("disp", "no_ace_neighbours_change_%03d.rds"))
saveRDS(point_post_ord_neighbours_change, write.path("disp", "point_postord_neighbours_change_%03d.rds"))
saveRDS(sample_post_ord_neighbours_change, write.path("disp", "sample_postord_neighbours_change_%03d.rds"))

cat("Neighbours completed...\n")



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
    }, error = function(e){
      cat("Error produced with sample post ord displacement:", e$message, "\n")
      NULL
    })
})

cat("Displacement completed...\n")

saveRDS(true_displacement_change, write.path("disp", "true_displacement_change_%03d.rds"))
saveRDS(strict_displacement_change, write.path("disp", "strict_displacement_change_%03d.rds"))
saveRDS(rel_displacement_change, write.path("disp", "rel_displacement_change_%03d.rds"))
saveRDS(sample_displacement_change, write.path("disp", "sample_displacement_change_%03d.rds"))
saveRDS(no_ace_displacement_change, write.path("disp", "no_ace_displacement_change_%03d.rds"))
saveRDS(point_post_ord_displacement_change, write.path("disp", "point_postord_displacement_change_%03d.rds"))
saveRDS(sample_post_ord_displacement_change, write.path("disp", "sample_postord_displacement_change_%03d.rds"))

