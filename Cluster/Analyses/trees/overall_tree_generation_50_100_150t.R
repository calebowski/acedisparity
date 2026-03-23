args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
set.seed(100 + replicate_id)
base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/trees/overallDisparity/")
source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")

## 50t
library(treats)
library(ape)

sizes <- c(50, 100, 150)
bd_params <- make.bd.params(speciation = 1, extinction = 0.90)
seed_val <- 100 + replicate_id

# stop_rule <- list(max.living = 50)
# set.seed(123)
metadata_list <- lapply(sizes, function(living_tips) {
  tree_size <- as.character(living_tips)
  max_total_size <- living_tips * 10
  
  # Keep generating trees until we get one within size limit
  repeat {
    # set.seed(seed_val)
    tree <- treats(stop.rule = list(max.living = living_tips), bd.params = bd_params, null.error = 1e6)
    tree <- drop.singles(tree)
    tree <- fix.zero.branches(tree)
    
    total_tips <- length(tree$tip.label)
    
    # Check if tree is within acceptable size
    if (total_tips <= max_total_size) {
      break  # Accept this tree
    } else {
      cat("Tree too large (", total_tips, " > ", max_total_size, "), regenerating...\n")
      # seed_val <- seed_val + 100  # Change seed for next attempt
    }
  }
  
  write.tree(tree, paste0(base_path, "tree_", tree_size, "t_", sprintf("%03d.tre", replicate_id)))

  est <- crude.bd.est(tree, "estimate")

  data.frame(
    replicate_id = replicate_id,
    living_size = living_tips,
    total_size = total_tips,
    seed = seed_val,
    speciation = est$call$speciation, 
    extinction = est$call$extinction,
    stringsAsFactors = FALSE
  )
})

metadata_df <- do.call(rbind, metadata_list)

metadata_file <- sprintf(file.path(base_path, "metadata", "metadata_%03d.csv"), replicate_id)
write.csv(metadata_df, metadata_file, row.names = FALSE)
