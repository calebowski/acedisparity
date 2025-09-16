args <- commandArgs(trailingOnly = TRUE)
replicate_id <- as.numeric(args[1])

################################################################################
#                                                                              #
#                           Random Extinction SIMULATION                       #
#                                                                              #
################################################################################

library(treats)
source("/users/bip24cns/acedisparity/discrete/scripts/utility.R")
source("/users/bip24cns/acedisparity/randomExtinction/scripts/find.extinction.time.R")


cat("Starting replicate", replicate_id, "\n")

set.seed(100 + replicate_id) # set seed to change for each replicate

bd_params <- make.bd.params(speciation = 0.1, extinction = 0.07)

stop_rule <- list(max.living = sample(x = c(50, 100, 150), size = 1)) # different tree sizes

slow_multi_transitions <- matrix(c(
    0.99, 0.025, 0.025,
    0.025, 0.99, 0.025,
    0.025, 0.025, 0.99
), nrow = 3, byrow = TRUE)


med_multi_transitions <- matrix(c(
  0.9, 0.05, 0.05,
  0.05, 0.9, 0.05,
  0.05, 0.05, 0.9
), nrow = 3, byrow = TRUE)


fast_multi_transitions <- matrix(c(
    0.1, 0.45, 0.45, 
    0.45, 0.1, 0.45,
    0.45, 0.45, 0.1 
), nrow = 3, byrow = TRUE)


slow_cor <- correlate.traits(correlation = 0.3, n.traits = 20, transition.matrix = slow_multi_transitions, n.correlated.traits = 5) ## high phylogenetic signal in correlation

med_cor <- correlate.traits(correlation = 0.3, n.traits = 50, transition.matrix = med_multi_transitions, n.correlated.traits = 5) 

fast_cor <- correlate.traits(correlation = 0.3, n.traits = 50, transition.matrix = fast_multi_transitions, n.correlated.traits = 5) ## low phylogenetic signal in the correlation



selective_extinction <- make.events(
                      target       = "taxa",
                      condition    = taxa.condition(30),
                      modification = trait.extinction(x = 2, trait = c(seq(from =  1, to = 7, by = 1)),
                                                      condition = `<`, severity = 0.8, threshold = 0.5))


selective <- lapply(multi_traits, function(rate) {
  treats(traits = rate, stop.rule = stop_rule, bd.params = bd_params, null.error = 10000, replicates = 30, events = selective_extinction)
})