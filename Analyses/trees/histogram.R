metadata <- read.csv("/home/caleb/Documents/PhD/acedisparity/Data/cluster/trees/large/metadata_all.csv")

metadata$net_div <- metadata$speciation - metadata$extinction


pdf("/home/caleb/Documents/PhD/acedisparity/Manuscript/draft/figures/net_div_hist_overall_trees.pdf")
hist(metadata$net_div, 
     xaxp = c(0, 1, 10), 
     main = "", 
     xlab = expression("Net Diversification (" * lambda * " - " * mu * ")"))
dev.off()
