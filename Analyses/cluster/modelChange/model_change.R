library(treats)
library(devtools)
devtools::load_all("/home/caleb/Documents/PhD/treats/treats")
tree <- read.tree("/home/caleb/Documents/PhD/acedisparity/Data/cluster/trees/large/tree_150t_001.tre")
tree <- set.root.time(tree)



# make.traits(process = BM.process, n = 20, 
#                      process.args = list(Sigma = diag(0.25, 20)))


tree_height <- max(node.depth.edgelength(tree))

ou_switch <- make.events(target = "traits", condition = age.condition(tree_height/2), modification = traits.update(process = OU.process, process.args = list(alpha = (log(2) / (tree_height * 0.10)), 
                                        Sigma = diag(0.25, 20))))

bm_switch <- make.events(target = "traits", condition = age.condition(tree_height/2), modification = traits.update(process = BM.process, process.args = list(Sigma = diag(0.25, 20))))

bm <- make.traits(process = BM.process, n = 20, process.args = list(Sigma = diag(0.25, 20)))                                  


ou_w  <- make.traits(process = OU.process, n = 20, 
                      process.args = list(alpha = (log(2) / (tree_height * 0.75)), 
                                        Sigma = diag(0.25, 20)))


switch_mapped <- map.traits(traits = ou_w, tree = tree, events = ou_switch)
plot(switch_mapped)


bm_t <- make.traits(process = BM.trend.process, n = 20, process.args = list(Sigma = diag(0.25, 20), trend = 0.3))




traits <- switch(model_name,
  "bm" = make.traits(process = BM.process, n = 20, 
                     process.args = list(Sigma = diag(0.25, 20))),
  
  "bm_t" = make.traits(process = BM.trend.process, n = 20, 
                      process.args = list(Sigma = diag(0.25, 20), trend = 0.3)),
  
  "ou_w" = make.traits(process = OU.process, n = 20, 
                      process.args = list(alpha = (log(2) / (tree_height * 0.75)), 
                                        Sigma = diag(0.25, 20))),
  
  "ou_st" = make.traits(process = OU.process, n = 20, 
                       process.args = list(alpha = (log(2) / (tree_height * 0.25)), 
                                         Sigma = diag(0.25, 20))),
  
  "ou_sh" = make.traits(process = OU.process, n = 20, 
                       process.args = list(optimum = 2, 
                                         alpha = (log(2) / (tree_height * 0.75)), 
                                         Sigma = diag(0.25, 20))),
  
  stop("Unknown model: ", model_name)
)