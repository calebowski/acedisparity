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
      file_path <- file.path("Data", "trees", 
                            sprintf("fossil_trees_%s_%03d.rds", size, i))
      if(file.exists(file_path)) {
        fossil_trees[[size]][[i]] <- readRDS(file_path)
      } else {
        warning("Missing file: ", file_path)
        fossil_trees[[size]][[i]] <- NULL
      }
    }
}



pdf("../Manuscript/draft/figures/tree_size_hist.pdf", width = 166 / 25.4, height = 83 / 25.4)
par(mfrow = c(1,3))

small_trees <- metadata[metadata$living_size == 50,]
hist(small_trees$total_size, main = "Tip size distribution (50 extant tips)", cex.main = 0.8, xlab = "Number of tips")

medium_trees <- metadata[metadata$living_size == 100,]
hist(medium_trees$total_size, main = "Tip size distribution (100 extant tips)", cex.main = 0.8, xlab = "Number of tips")

large_trees <- metadata[metadata$living_size == 150,]
hist(large_trees$total_size, main = "Tip size distribution (150 extant tips)", cex.main = 0.8, xlab = "Number of tips")
dev.off()

colless_indices_2 <- lapply(fossil_trees, lapply, lapply, lapply, colless_corr)

sampling_levels <- c("all", "fossil_high", "fossil_med", "fossil_low", "living")
level_titles <- c("all" = "100%", "fossil_high" = "50%", "fossil_med" = "15%", "fossil_low" = "5%", "living" = "0%")

pdf("Manuscript/draft/figures/symmetry.pdf", width = 15, height = 9)
par(mfrow = c(length(tree_sizes), length(sampling_levels)), mar = c(4, 4, 3, 1))

for (size in tree_sizes) {
  for (level in sampling_levels) {
    vals <- unlist(lapply(colless_indices_2[[size]], function(rep) {
      rep[[1]][[level]]
    }))
    
    hist(vals, 
          main = paste(size, "-", level_titles[level]), 
          xlab = "Corrected Colless Index",
          # col = "grey80", border = "white",
          cex.main = 1.8,  
          cex.lab  = 1.6, 
          cex.axis = 1.4)
  }
}
dev.off()
