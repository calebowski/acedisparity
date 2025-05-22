png("Manuscript/images/phylo_tree.png", width=800, height=600)
# par(mar = c(5, 4, 6, 2)) 
tree_newick <- "((A:0.2,B:0.2)internal:0.4,C:0.6);"
 
    # Read tree
    tree <- read.tree(text = tree_newick)

    tree$node.label <- c("root", "AB")
 
    # Plot the tree
    plot(tree, show.tip.label = TRUE, show.node.label = TRUE, edge.width = 2, cex = 2, label.offset = 0.02)
 
    # Add branch lengths as text labels
    edgelabels(round(tree$edge.length, 2), frame = "none", adj = c(1, -0.5), cex =1.2)
 
    # Define tip states: A=1 (black), B=0 (white), C=0 (white)
    tip_states <- c(A=1, B=0, C=0)
 
    # Add colored squares for tip states
    tiplabels(pch=22, col="black", bg=ifelse(tip_states == 1, "black", "white"), cex=4)

    legend(
  "topleft",                      
  legend = c("1", "0"),           
  pch = 22,                        
  pt.bg = c("black", "white"),   
  col = "black",                   
  pt.cex = 4,                      
  title = "Tip state"
)



# Close the device to save the file
dev.off()



# mat <- matrix(c(-0.5, 0.5, 0.5, -0.5), ncol = 2, nrow = 2, byrow = TRUE )

# exp <- expm(mat * 0.4)
