args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])

library(treats)

base_path <- "/mnt/parscratch/users/bip24cns/acedisparity/continuous/100t/"
job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")

write.path <- function(subfolder, filename) {
  paste0(base_path, subfolder, "/", job_id, "_", sprintf(filename, replicate_id))
}

set.seed(100 + replicate_id)

bd_params <- make.bd.params(speciation = 1, extinction = 0.7)
stop_rule <- list(max.living = 100)

tree <- treats(stop.rule = stop_rule, bd.params = bd_params, null.error = 100)
tree <- drop.singles(tree)

write.tree(tree, write.path("trees", "tree_%03d.tre"))

b_d_est <- crude.bd.est(tree, "estimate")
metadata_df <- data.frame(
  replicate_id = replicate_id,
  tree_size = length(tree$tip.label),
  seed = 100 + replicate_id,
  speciation = b_d_est$call$speciation,
  extinction = b_d_est$call$extinction
)
write.csv(metadata_df, write.path("metadata", "metadata_%03d.csv"), row.names = FALSE)

cat("Tree generated for replicate", replicate_id, "\n")