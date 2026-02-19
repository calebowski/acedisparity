## install R packages 
if (!require(devtools, quietly = TRUE)) {
  install.packages("devtools")
  library(devtools)
}


install_github("TGuillerme/treats", ref = "map.traits.events")
library(treats)


dir.create("../Data/trees", recursive = TRUE, showWarnings = FALSE)
dir.create("../Data/trees/metadata", recursive = TRUE, showWarnings = FALSE)

source("../Functions/utility.R")

sizes <- c(50, 100, 150) ## three tree sizes
bd_params <- make.bd.params(speciation = 1, extinction = 0.90)
# seed_val <- 100 + replicate_id

# stop_rule <- list(max.living = 50)
# set.seed(123)

n_replicates <- 100
for(i in 1:n_replicates) {
  replicate_id <- as.numeric(i)
  seed <- 100 + replicate_id
  set.seed(seed)

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
    tree_file <- sprintf("../Data/trees/tree_%s_%03d.tre", tree_size, replicate_id)
    write.tree(tree, tree_file)

    est <- crude.bd.est(tree, "estimate")

    data.frame(
      replicate_id = replicate_id,
      living_size = living_tips,
      total_size = total_tips,
      seed = seed,
      speciation = est$call$speciation, 
      extinction = est$call$extinction,
      stringsAsFactors = FALSE
    )
  })
  metadata_df <- do.call(rbind, metadata_list)

  metadata_file <- sprintf(file.path("../Data/trees", "metadata", "metadata_%03d.csv"), replicate_id)
  write.csv(metadata_df, metadata_file, row.names = FALSE)
}

# Aggregate all metadata files
all_metadata <- do.call(rbind, lapply(1:n_replicates, function(i) {
  read.csv(sprintf("../Data/trees/metadata/metadata_%03d.csv", i))
}))

# Save combined metadata
write.csv(all_metadata, "../Data/trees/metadata/metadata_all.csv", row.names = FALSE)

cat("\nTree generation complete!\n")