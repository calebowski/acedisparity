metadata <- read.csv("/home/caleb/Documents/PhD/acedisparity/Data/cluster/trees/large/metadata_all.csv")

metadata$net_div <- metadata$speciation - metadata$extinction


pdf("/home/caleb/Documents/PhD/acedisparity/Manuscript/draft/figures/net_div_hist_overall_trees.pdf")
hist(metadata$net_div, 
     xaxp = c(0, 1, 10), 
     main = "", 
     xlab = expression("Net Diversification (" * lambda * " - " * mu * ")"))
dev.off()

library(treestats)
library(apTreeshape)
library(ape)

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


colless_indices_2 <- lapply(fossil_trees, lapply, lapply, lapply, colless_corr)
pdf("../Manuscript/draft/figures/symmetry.pdf")
hist(unlist(colless_indices_2))
dev.off()
