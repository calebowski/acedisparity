library(treats)
library(phytools)
sizes <- c(50, 100, 150)
bd_params <- make.bd.params(speciation = 1, extinction = 0.90)
traits <- make.traits(process = BM.process)
tree <- treats(traits = traits, stop.rule = list(max.living = 150), bd.params = bd_params, null.error = 100)

tip_mat <- tree$data[grepl("^t", rownames(tree$data)),]

Sys.time()
test <- multi.ace(as.data.frame(tip_mat), tree$tree, models = "BM", output = "multi.ace")
Sys.time()
test <- fastAnc(tree$tree, as.data.frame(tip_mat), CI = TRUE)
Sys.time()
## takes 3 minutes on a 150 large tree
## for one trait
3 * 1498
Sys.time()


Sys.time()
test <- multi.ace(as.data.frame(tip_mat), tree$tree, models = "ML", output = "multi.ace")
Sys.time()