# correlate.traits.binary <- function(correlation = 0.9, transition.matrix, n.traits = 5) {
#   trait1 <- make.traits( # create first parent trait
#     process = discrete.process,
#     process.args = list(transitions = transition.matrix),
#     trait.names = "trait1",
#     n = n.traits
#   )
  
#   # trait2 depends on trait1 with specified correlation strength
#   trait2_0 <- make.traits(
#     process = function(x0 = 0, edge.length = 1) {
#       # When trait1 = 0, trait2 follows correlation pattern
#       rbinom(1, 1, prob = 1 - correlation)  # Low correlation = more likely to be different
#     }, 
#     trait.names = "trait2", 
#     n = n.traits
#   )
  
#   trait2_1 <- make.traits(
#     process = function(x0 = 0, edge.length = 1) {
#       # When trait1 = 1, trait2 follows correlation pattern  
#       rbinom(1, 1, prob = correlation)  # High correlation = more likely to be same
#     }, 
#     trait.names = "trait2", 
#     n = n.traits
#   )
  
#   link_args <- list(
#     choose_0 = function(x) x == 0,
#     choose_1 = function(x) x == 1
#   )
  
#   linked_traits <- link.traits(
#     base.trait = trait1,
#     next.trait = list(trait2_0, trait2_1),
#     link.type = "conditional",
#     link.args = link_args
#   )
  
#   return(linked_traits)
# }


# correlate.traits.multi <- function(correlation = 0.9, transition.matrix, n.traits = 5) { ## note that currently this works for multistate, can adapt it for binary.
#   trait1 <- make.traits( # create first parent trait
#     process = discrete.process,
#     process.args = list(transitions = transition.matrix),
#     n = n.traits
#   )
  
#   # trait2 depends on trait1 with specified correlation strength
#   # When trait1 = 0, trait2 should mostly be 0
#   trait2_0 <- make.traits(
#     process = function(x0 = 0, edge.length = 1) {
#       probs <- c(correlation, (1-correlation)/2, (1-correlation)/2)  # Favor state 0
#       sample(0:2, 1, prob = probs)
#     }, 
#     n = n.traits
#   )
  
#   # When trait1 = 1, trait2 should mostly be 1
#   trait2_1 <- make.traits(
#     process = function(x0 = 0, edge.length = 1) {
#       probs <- c((1-correlation)/2, correlation, (1-correlation)/2)  # Favor state 1
#       sample(0:2, 1, prob = probs)
#     }, 
#     n = n.traits
#   )
  
#   # When trait1 = 2, trait2 should mostly be 2
#   trait2_2 <- make.traits(
#     process = function(x0 = 0, edge.length = 1) {
#       probs <- c((1-correlation)/2, (1-correlation)/2, correlation)  # Favor state 2
#       sample(0:2, 1, prob = probs)
#     }, 
#     n = n.traits
#   )
  
#   link_args <- list(
#     choose_0 = function(x) x == 0,
#     choose_1 = function(x) x == 1,
#     choose_2 = function(x) x == 2
#   )
  
#   linked_traits <- link.traits(
#     base.trait = trait1,
#     next.trait = list(trait2_0, trait2_1, trait2_2),
#     link.type = "conditional",
#     link.args = link_args
#   )
  
#   return(linked_traits)
# }

correlate.traits <- function(correlation = 0.9, transition.matrix, n.traits = 6, n.correlated.traits = 3) {
  
  # Create a custom process that applies correlation more realistically
  correlate.process <- function(x0 = rep(0, n.traits), edge.length = 1, 
                                          correlation = correlation, 
                                          transitions = transition.matrix) {
    
    # If x0 is a single value, expand it to vector
    if(length(x0) == 1) {
      x0 <- rep(x0, n.traits)
    }
    
    # Generate all traits independently first
    independent_traits <- sapply(1:n.traits, function(i) {
      discrete.process(x0 = x0[i], edge.length = edge.length, transitions = transitions)
    })
    
    # Apply correlation by modifying some traits to match the first one
    trait_values <- independent_traits
    
    for(i in 2:n.correlated.traits) {
      # Decide if this trait should be correlated
      if(runif(1) < correlation) {
        # Make it the same as trait 1
        trait_values[i] <- trait_values[1]
      }
      # Otherwise keep the independent evolution
    }
    
    return(trait_values)
  }
  
  # Create traits object
  make.traits(
    process = correlate.process,
    process.args = list(correlation = correlation, transitions = transition.matrix),
    n = n.traits
  )
}
