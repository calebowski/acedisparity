create.matched.data <- function(true_vals, ace_vals_clean) {
  true_matched <- c()
  ace_matched <- c()
  fossil_levels <- c()  # Add this line
  
  for(method_name in names(ace_vals_clean)) {
    for(rate_name in names(ace_vals_clean[[method_name]])) {
      for(rep_idx in seq_along(ace_vals_clean[[method_name]][[rate_name]])) {
        
        true_val <- true_vals[[rate_name]][[rep_idx]]
        ace_rep_vals <- ace_vals_clean[[method_name]][[rate_name]][[rep_idx]]
        
        # Get fossil level names and repeat true values
        fossil_names <- names(ace_rep_vals)
        n_fossil_levels <- length(ace_rep_vals)
        
        true_matched <- c(true_matched, rep(true_val, n_fossil_levels))
        ace_matched <- c(ace_matched, unlist(ace_rep_vals))
        fossil_levels <- c(fossil_levels, fossil_names)  # Add this line
      }
    }
  }
  
  return(list(true = true_matched, ace = ace_matched, fossil_level = fossil_levels))
}

create.correlated.traits <- function(correlation = 0.9, transition.matrix, n.traits = 5) {
  trait1 <- make.traits( # create first parent trait
    process = discrete.process,
    process.args = list(transitions = transition.matrix),
    trait.names = "trait1",
    n = n.traits
  )
  
  # trait2 depends on trait1 with specified correlation strength
  trait2_0 <- make.traits(
    process = function(x0 = 0, edge.length = 1) {
      # When trait1 = 0, trait2 follows correlation pattern
      rbinom(1, 1, prob = 1 - correlation)  # Low correlation = more likely to be different
    }, 
    trait.names = "trait2", 
    n = n.traits
  )
  
  trait2_1 <- make.traits(
    process = function(x0 = 0, edge.length = 1) {
      # When trait1 = 1, trait2 follows correlation pattern  
      rbinom(1, 1, prob = correlation)  # High correlation = more likely to be same
    }, 
    trait.names = "trait2", 
    n = n.traits
  )
  
  link_args <- list(
    choose_0 = function(x) x == 0,
    choose_1 = function(x) x == 1
  )
  
  linked_traits <- link.traits(
    base.trait = trait1,
    next.trait = list(trait2_0, trait2_1),
    link.type = "conditional",
    link.args = link_args
  )
  
  return(linked_traits)
}