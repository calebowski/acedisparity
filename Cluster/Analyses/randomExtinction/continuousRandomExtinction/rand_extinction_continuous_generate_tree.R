args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
tree_size <- args[2]

library(treats)

base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/randomExtinction/continuous/", tree_size, "t/")
job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")

write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}

set.seed(100 + replicate_id) # set seed to change for each replicate

bd_params <- make.bd.params(speciation = 1.0, extinction = 0.7)

stop_rule <- list(max.living = as.numeric(tree_size)) # different tree sizes

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
  
  if(n_tips <= 200) {
    cat("Tree found with", n_tips, "tips\n")
    break
  }
  
  if(attempt >= max_attempts) {
    cat("Warning: Reached maximum attempts, accepting tree with", n_tips, "tips\n")
    break
  }
  
  attempt <- attempt + 1
}

tree <- drop.singles(tree)
write.tree(tree, write.path("trees", "tree_%03d.tre"))

metadata_df <- data.frame(
  replicate_id = replicate_id,
  tree_size = stop_rule$max.living,
  extinction_intensity = extinction_intensity,
  seed = 100 + replicate_id
)

write.csv(metadata_df, 
          write.path("metadata", "metadata_%03d.csv"),
          row.names = FALSE)
