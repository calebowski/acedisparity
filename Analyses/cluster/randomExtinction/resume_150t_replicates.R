# Resume script for post-ordination ACE - replicates 82, 96, 98
args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])

library(treats)
library(dispRity)
library(parallel)

source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")
source("/users/bip24cns/acedisparity/randomExtinction/scripts/find.extinction.time.R")

base_path <- "/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/150t/"
job_id <- 8153183

write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}

cat("Resuming post-ordination ACE for replicate", replicate_id, "\n")

# Load required saved data
ord_no_ace <- readRDS(write.path("ord", "ord_no_ace_%03d.rds"))
fossil_trees <- readRDS(write.path("trees", "fossil_trees_%03d.rds"))
tree <- read.tree(write.path("trees", "rand_ext_tree_%03d.tre"))

tree <- drop.singles(tree)
tree <- fix.zero.branches(tree)
tree <- set.root.time(tree)

# Resume from post-ordination ACE
# cat("Starting post-ordination ACE...\n")
cat("Starting post-ordination ACE...\n")

# Set up parallel processing for 5 fossil levels
# Function to flatten and parallelize all 15 tasks
run_ace_all_parallel <- function(matrices, trees, models = "BM", n_cores = 15) {
  
  # Create task list for all 15 combinations
  tasks <- expand.grid(
    rate = names(matrices),
    fossil = names(matrices[[1]]),
    stringsAsFactors = FALSE
  )
  
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl))
  
  clusterEvalQ(cl, library(treats))
  clusterExport(cl, c("matrices", "trees"), envir = environment())
  
  cat("Processing", nrow(tasks), "ACE tasks with", length(cl), "cores...\n")
  
  # Run all 15 tasks in parallel
  results <- parLapply(cl, 1:nrow(tasks), function(i) {
    rate <- tasks$rate[i]
    fossil <- tasks$fossil[i]
    multi.ace(matrices[[rate]][[fossil]], 
              trees[[rate]][[fossil]], 
              models = models, 
              output = "multi.ace")
  })
  
  # Reconstruct nested structure
  output <- setNames(vector("list", length(names(matrices))), names(matrices))
  for(i in 1:nrow(tasks)) {
    rate <- tasks$rate[i]
    fossil <- tasks$fossil[i]
    if(is.null(output[[rate]])) output[[rate]] <- list()
    output[[rate]][[fossil]] <- results[[i]]
  }
  
  return(output)
}

# Replace your current post-ordination ACE section with:
n_cores <- 15  # Match your 3 models × 5 fossil levels
cat("Starting post-ordination ACE with", n_cores, "cores (15 total tasks)...\n")

post_ord_ace <- run_ace_all_parallel(ord_no_ace, fossil_trees, models = "BM", n_cores = n_cores)

cat("Post-ordination ACE completed with parallel processing...\n")


trait_normal <- list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))

point_post_ord_ace <- lapply(post_ord_ace, lapply, multi.ace, output = "combined.matrix")

sample_post_ord_ace <- lapply(post_ord_ace, lapply, function(x){
  multi.ace(x, sample = 100, sample.fun = trait_normal, output = "combined.matrix")
})

cat("Post ordination estimates completed...\n")

saveRDS(point_post_ord_ace, write.path("ord", "post_ord_point_%03d.rds"))
saveRDS(sample_post_ord_ace, write.path("ord", "post_ord_sample_%03d.rds"))

cat("Post ordination estimates saved...\n")

# Continue with chrono subsets and disparity analysis
tip_ages <- tip.ages(tree)
extinction_time <- find.extinction.time(tip_ages)

tree_height <- max(node.depth.edgelength(tree))
stage <- tree_height / 6
time_slices <- c(extinction_time + stage, extinction_time - 0.01, (extinction_time - 0.01) - stage)

# Load previously computed ordinations for chrono subsets
ord_strict <- readRDS(write.path("ord", "ord_strict_%03d.rds"))
ord_rel <- readRDS(write.path("ord", "ord_rel_%03d.rds"))
ord_sample <- readRDS(write.path("ord", "ord_sample_%03d.rds"))
ord_true <- readRDS(write.path("ord", "ord_true_%03d.rds"))


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



cat("Finished replicate", replicate_id, "\n")