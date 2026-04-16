metadata <- read.csv("../Data/trees/metadata/metadata_all.csv")

metadata$net_div <- metadata$speciation - metadata$extinction


pdf("../Manuscript/draft/figures/net_div_hist_overall_trees.pdf")
hist(metadata$net_div, 
     xaxp = c(0, 1, 10), 
     main = "", 
     xlab = expression("Net Diversification (" * lambda * " - " * mu * ")"))
dev.off()

library(treestats)
library(ape)

tree_sizes <- c("50t", "100t", "150t")

fossil_trees <- list()
for (size in tree_sizes) {
  fossil_trees[[size]] <- list()
    for(i in 1:100) {
      file_path <- file.path("..", "Data", "trees", 
                            sprintf("fossil_trees_%st_%03d.rds", size, i))
      if(file.exists(file_path)) {
        fossil_trees[[size]][[i]] <- readRDS(file_path)
      } else {
        warning("Missing file: ", file_path)
        fossil_trees[[size]][[i]] <- NULL
      }
    }
}


colless_indices_2 <- lapply(fossil_trees, lapply, lapply, lapply, colless_corr)







pdf("../Manuscript/draft/figures/symmetry.pdf")
par(mfrow = c(3,2))
hist(unlist(colless_indices_2))
dev.off()
