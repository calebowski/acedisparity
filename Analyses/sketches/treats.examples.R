library(treats)
## Making asymmetric trees by modifying lineage selection:

## Base stuff
bd_params <- make.bd.params(speciation = 1, extinction = 0.9)
stop_rule <- list("max.living" = 50)


## Biased selection
bias.select <- function(lineage) {
    ## Sample one lineage among the existing lineages
    ## "lineage" is an internal treats object details in manual.
    ## The sample proportion is a decreasing log exponential distribution
    ## i.e. the last lineage has always more changes to be selected
    ## You can modify the rate parameter internally. Bigger = more ladderised
    probs <- rev(dexp(seq(1, lineage$n, by = 1), rate = 0.5))
    return(sample(1:lineage$n, 1, prob = probs))
}



## A modifier for selection
biased_modifier <- make.modifiers(selection = bias.select)

## Generating a new tree with this modifier
set.seed(1)
ladderised_tree <- treats(bd.params  = bd_params,
                          stop.rule  = stop_rule,
                          modifiers  = biased_modifier,
                      	  null.error = 100)
## Displaying the results
plot(ladderize(ladderised_tree), main = "Asymmetric-ish tree")





## Making fossils as near-sampled nodes
bd_params <- make.bd.params(speciation = 1, extinction = 0.9)
stop_rule <- list("max.living" = 50)
my_trait <- make.traits()

## Generating a new tree with this modifier
set.seed(2)
tree_with_steps <- treats(bd.params  = bd_params,
                          stop.rule  = stop_rule,
                          traits     = my_trait,
                      	  null.error = 100,
                      	  save.steps = 1)
## Displaying the results
plot(tree_with_steps, main = "Saved_steps")
plot(drop.singles(tree_with_steps), main = "Saved_steps (but dropped internal branches)")

## This one needs a bit more thinking but basically you'd want to collapse the fossil branches to saved steps the nearest to their node.
## If you plot the tree with all the internal nodes:
plot(ladderize(tree_with_steps$tree))



nodelabels(".", cex = 0.2)
# mapped <- map.traits(traits = make.traits(BM.process),tree_with_steps$tree)

## You can see for example t1 has three ticks before reaching it's node.
## There's probably some smart way to do this by modify the tree$edge table but I don't have time to focus on that right now.
## Would be fun to hackathon though!





