args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])
model_name <- args[2] 
tree_size <- args[3]

cat("=== PROCESS TEST START ===\n")
cat("Replicate:", replicate_id, "\n")
cat("Model:", model_name, "\n")
cat("Tree size:", tree_size, "\n")
cat("Process ID:", Sys.getpid(), "\n")
cat("Time:", Sys.time(), "\n")

library(ape)
library(treats)

base_path <- paste0("/mnt/parscratch/users/bip24cns/acedisparity/continuous/", tree_size, "/")
job_id <- 8469104

write.path <- function(subfolder, filename) {
  return(paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id)))
}

# Test 1: Read tree
cat("\n=== STEP 1: Reading tree ===\n")
tree_file <- write.path("trees", "tree_%03d.tre")
cat("Tree file:", tree_file, "\n")

tryCatch({
  tree <- read.tree(tree_file)
  cat("SUCCESS: Read tree with", length(tree$tip.label), "tips\n")
}, error = function(e) {
  cat("ERROR reading tree:", e$message, "\n")
  quit(status = 1)
})

# Test 2: Simple trait simulation (just 5 traits for speed)
cat("\n=== STEP 2: Testing trait simulation ===\n")
set.seed(replicate_id * 1000 + which(c("bm", "bm_t", "ou_st", "ou_w", "ou_sh") == model_name))

tryCatch({
  if(model_name == "bm") {
    traits <- treats(bd.params = NULL, traits = list(n = 5, process = "BM", process.args = list(Sigma = 1)))
  } else if(model_name == "bm_t") {
    traits <- treats(bd.params = NULL, traits = list(n = 5, process = "BM", process.args = list(Sigma = 1, trends = rnorm(5))))
  } else {
    # Simplified OU for testing
    traits <- treats(bd.params = NULL, traits = list(n = 5, process = "OU", process.args = list(A = diag(5))))
  }
  
  cat("SUCCESS: Created traits with", length(traits), "dimensions\n")
}, error = function(e) {
  cat("ERROR creating traits:", e$message, "\n")
  quit(status = 1)
})

# Test 3: Map traits to tree
cat("\n=== STEP 3: Testing trait mapping ===\n")
tryCatch({
  trait_data <- map.traits(traits, tree)
  cat("SUCCESS: Mapped traits to tree\n")
  cat("Trait matrix dimensions:", dim(trait_data), "\n")
}, error = function(e) {
  cat("ERROR mapping traits:", e$message, "\n")
  quit(status = 1)
})

cat("\n=== PROCESS TEST COMPLETE ===\n")
cat("Model", model_name, "completed successfully for replicate", replicate_id, "\n")
cat("End time:", Sys.time(), "\n")