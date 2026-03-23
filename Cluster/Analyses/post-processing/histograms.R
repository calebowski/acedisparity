library(treestats)
library(apTreeshape)
library(ape)
job_ids <- list(
  "50t" = "8558401",
  "100t" = "8561556",  # Replace with actual job ID
  "150t" = "8561712"   # Replace with actual job ID
)

tree_sizes <- c("50t", "100t", "150t")


max_colless <- function(n) {
  return( ((n-1)*(n-2)) / 2 )
}

fossil_matrices <- list()
for (size in tree_sizes) {
  fossil_matrices[[size]] <- list()
  job_id <- job_ids[[size]]  #  Get job ID for this tree size
    for(i in 1:100) {
      file_path <- file.path("..", "Data", "cluster", "discrete_crown", size, 
                            "matrices",
                            sprintf("%s_fossil_matrices_%03d.rds", job_id, i))
      if(file.exists(file_path)) {
        fossil_matrices[[size]][[i]] <- readRDS(file_path)
      } else {
        warning("Missing file: ", file_path)
        fossil_matrices[[size]][[i]] <- NULL
      }
    }
}

fossil_trees <- lapply(fossil_matrices, lapply, lapply, lapply, function(x){x$tree})
ts <- lapply(fossil_trees, lapply, lapply, lapply, as.treeshape)

get.max.colless <- function(n) {
  return(as.numeric((n-1)*(n-2)) / 2)
}


normalised.colless <- function(tree) {
    n <- length(tree$tip.label)
    max_colless <- get.max.colless(n)
    raw_score <- treestats::colless(tree, normalization = "none")
    norm_score <- raw_score / max_colless
}

colless_indices <- lapply(fossil_trees, lapply, lapply, lapply, normalised.colless)
colless_indices_2 <- lapply(fossil_trees, lapply, lapply, lapply, colless_corr)
pdf("../Manuscript/draft/figures/symmetry.pdf")
hist(unlist(colless_indices_2))
dev.off()


# ## PATTERNS IN TREE BALANCE AMONG CLADISTIC,
# PHENETIC, AND RANDOMLY GENERATED
# PHYLOGENETIC TREES Stephen Heard is the paper to cite for this