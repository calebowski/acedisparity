library(treats)

trait <- make.traits(BM.process, trait.names = "Trait value")
bd_params <- make.bd.params(speciation = 1, extinction = 0.7)
stop_rule <- list(max.living = 20)

# tree <- rtree(n = )

trees <- treats(traits = trait, stop.rule = stop_rule, bd.params = bd_params, null.error = 100, replicates = 5)

## First plot is just tips and fossils, no edges or 
png(filename = "../besMacro/BMnoedge.png", width = 780, height = 650)
plot.treats(trees[[3]], legend = TRUE, cex = 2, cex.axis = 1.2, legend.cex = 1.5, cex.lab = 1.4, col = c("nodes" = "NA"), edges = NULL)
dev.off()

png(filename = "../besMacro/BMnonode.png", width = 780, height = 650)
plot.treats(trees[[2]], legend = TRUE, cex = 2, cex.axis = 1.2, legend.cex = 1.5, cex.lab = 1.4, col = c("nodes" = "NA"), edges = "black")
dev.off()


png(filename = "../besMacro/BMfull.png", width = 780, height = 650)
plot.treats(trees[[5]], legend = TRUE, cex = 2, cex.axis = 1.2, legend.cex = 1.5, cex.lab = 1.4, edges = "black")
dev.off()


# plot.treats(trees[[2]], legend = TRUE, cex = 2, cex.axis = 1.2, legend.cex = 1.5, cex.lab = 1.4, col = c("nodes" = "NA"), edges = NULL)

discrete <- make.traits(discrete.process)
discrete_trees <- treats(traits = discrete, stop.rule = stop_rule, bd.params = bd_params, null.error = 100, replicates = 5)
plot.treats(discrete_trees[[1]], legend = TRUE, cex = 2, cex.axis = 1.2, legend.cex = 1.5, cex.lab = 1.4, tips.nodes = "blue")



###########################################################################################
# Make figures for different fossil sampling levels
bd_params <- make.bd.params(speciation = 1, extinction = 0.7)

stop_rule <- list(max.living = 8)
set.seed(123)
trees <- treats(stop.rule = stop_rule, bd.params = bd_params, null.error = 100, replicates = 3)
full_fossil_tree <- trees[[3]]
ages <- tree.age(tree)
living <- ages$element[ages$ages == 0]
living_tree <- keep.tip(full_fossil_tree, c(living, "t1"))




## What if use a similar logic, but try and project the estimated values onto the real treats plot to compare how far they are.

