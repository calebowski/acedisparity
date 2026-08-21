remove.fossil <- function(trees, matrices, type = c("discrete", "continuous")) {
  
  # Base function for processing one tree and one matrix (discrete or continuous)
  process.living <- function(tree, matrix, type) {
    ages <- tree.age(tree) # get tip ages
    tips <- ages$element[ages$ages == 0] # find only living species
    
    living_matrix <- matrix[rownames(matrix) %in% tips, , drop = FALSE]
    
    # Check if we have any living species
    if (nrow(living_matrix) == 0) {
      warning("No living species found")
      return(list(matrix = living_matrix, tree = NULL))
    }
    
    if (type == "discrete") {
      living_matrix <- apply(living_matrix, c(1, 2), as.character)
    }
    
    # keep only living species in tree
    living_tree <- keep.tip(tree, tips)
    # tree <- set.root.time(tree)
    # old_root <- tree$root.time
    # living_tree$root.time <- old_root
    
    return(list(matrix = living_matrix, tree = living_tree))
  }
  
  # Adapt code for when there is nested list structures
  if (is.list(trees) && is.list(matrices)) {
    
    living_matrices <- Map(function(tree, matrix) process.living(tree, matrix, type), trees, matrices)
  } else {
    
    living_matrices <- process.living(trees, matrices, type)
  }
  
  return(living_matrices)
}


fossil.pres <- function(trees, matrices, preservation = c(0.05, 0.15, 0.5, 1.0), type = c("discrete", "continuous"), seed = NULL) {
  process.fossil <- function(tree, matrix, type, seed) {
    ages <- tree.age(tree)
    tips <- ages$element[ages$ages == 0]  # Keep living species
    fossils <- ages$element[ages$ages > 0 & grepl("^t", ages$element)]
    set.seed(seed)
    max_attempts <- 20
    attempt <- 1
    repeat {

      if(length(fossils) == 0){
        sample_fossil  <- character(0) 
        cat("No fossils in tree... \n")
        break

      } else {

      keep_vector <- as.logical(rbinom(length(fossils), size = 1, prob = preservation)) ## uses Bernoulli: samples fossils independently
      sample_fossil  <- fossils[keep_vector]
        if (length(sample_fossil) > 0) break

      }

      if(attempt >= max_attempts) {
        sample_fossil <- sample(fossils, 1) ## if reach max attempts just use one fossil
      cat("Warning: Reached maximum attempts, accepting 1 fossil... \n")
      break
      }

      attempt <- attempt + 1

    }

    kept <- c(unlist(unname(sample_fossil)), tips)
    fossil_matrix <- matrix[rownames(matrix) %in% kept, ]

    # Discrete needs characters returned for ace, continuous needs numeric
    if (type == "discrete") {
      fossil_matrix <- apply(fossil_matrix, c(1, 2), as.character)
    } 

    if(type == "continuous") {
      fossil_matrix <- apply(fossil_matrix, c(1, 2), as.numeric)
      # fossil_matrix <- as.data.frame(fossil_matrix)
    }


    pruned <- keep.tip(tree, kept)
    return(list(matrix = fossil_matrix, tree = pruned))
  }

  # Check if inputs are lists or single objects
  if (is.list(trees) && is.list(matrices)) {
    # Use Map for lists of trees and matrices
    fossil_matrices <- Map(function(tree, matrix) process.fossil(tree, matrix, type, seed), trees, matrices)
  } else {
    # Process a single tree and matrix
    fossil_matrices <- process.fossil(trees, matrices, type, seed)
  }

  return(fossil_matrices)
}


## preservation should be a vector of two values
# non.random.fossil.pres <- function(trees, matrices, preservation = c(0.05, 0.5), type = c("discrete", "continuous"), seed = NULL) {
#   process.fossil <- function(tree, matrix, type, seed) {
    
#     root_node <- tree$node.label[1]
#     clade_1 <- extract.clade(tree, tree$node.label[2])
#     clade_2 <- extract.clade(tree, tree$node.label[3])
#     clades <- list(clade_1, clade_2)
#     tips <- lapply(clades, function(clade){
#       ages <- tree.age(clade)
#       tips <- ages$element[grepl("^t", ages$element)]
#       # fossils <- ages$element[grepl("^t", ages$element)]
#       return(tips)
#     }) ## get ages of trees
#     # living_tips_clade_1 <- sample()
#     # living_tips <- subset(tree.age(tree), ages == 0)$elements

#     set.seed(seed)
#     max_attempts <- 20
#     attempt <- 1
#     repeat {

#       if(length(tips[[1]]) == 0 || length(tips[[2]]) == 0){
#         sample_fossil  <- character(0) 
#         cat("No fossils in tree... \n")
#         break

#       } else {

#       clade_1_keep <- as.logical(rbinom(length(tips[[1]]), size = 1, prob = preservation[1])) ## uses Bernoulli: samples tips independently
#       clade_2_keep <- as.logical(rbinom(length(tips[[2]]), size = 1, prob = preservation[2])) ## uses Bernoulli: samples tips independently

#       sample_fossil_1  <- tips[[1]][clade_1_keep]
#       sample_fossil_2 <- tips[[2]][clade_2_keep]
#         if (length(sample_fossil_1) > 0 && length(sample_fossil_2) > 0) break
#       }

#       # if(attempt >= max_attempts) {
#       #   sample_fossil <- sample(fossils, 1) ## if reach max attempts just use one fossil
#       # cat("Warning: Reached maximum attempts, accepting 1 fossil... \n")
#       # break
#       # }
#       # attempt <- attempt + 1
#     }

#     kept <- c(unlist(c(sample_fossil_1, sample_fossil_2)))
#     fossil_matrix <- matrix[rownames(matrix) %in% kept, ]

#     # Discrete needs characters returned for ace, continuous needs numeric
#     if (type == "discrete") {
#       fossil_matrix <- apply(fossil_matrix, c(1, 2), as.character)
#     } 

#     if(type == "continuous") {
#       fossil_matrix <- apply(fossil_matrix, c(1, 2), as.numeric)
#       # fossil_matrix <- as.data.frame(fossil_matrix)
#     }


#     pruned <- keep.tip(tree, kept)
#     return(list(matrix = fossil_matrix, tree = pruned))
#   }

#   # Check if inputs are lists or single objects
#   if (is.list(trees) && is.list(matrices)) {
#     # Use Map for lists of trees and matrices
#     fossil_matrices <- Map(function(tree, matrix) process.fossil(tree, matrix, type, seed), trees, matrices)
#   } else {
#     # Process a single tree and matrix
#     fossil_matrices <- process.fossil(trees, matrices, type, seed)
#   }

#   return(fossil_matrices)
# }




fossil.pres.alt<- function(trees, matrices, preservation = c(0.05, 0.15, 0.5, 1.0), type = c("discrete", "continuous"), seed = NULL) {
  process.fossil.alt <- function(tree, matrix, type, seed) {
    ages <- tree.age(tree)
    tips <- ages$element[ages$ages == 0]  # Keep living species
    fossils <- ages$element[ages$ages > 0]
    set.seed(seed)
    max_attempts <- 20
    attempt <- 1
    repeat {

      # if(length(fossils) == 0){
      #   sample_fossil  <- character(0) 
      #   cat("No fossils (extinct tips) found in tree... \n")
      #   break

      # } else {

      keep_vector <- as.logical(rbinom(length(fossils), size = 1, prob = preservation)) ## uses Bernoulli: samples fossils independently
      sample_fossil  <- fossils[keep_vector]
        if (length(sample_fossil) > 0) break

      

      if(attempt >= max_attempts) {
        sample_fossil <- sample(fossils, 1) ## if reach max attempts just use one fossil
      cat("Warning: Reached maximum attempts, accepting 1 fossil... \n")
      break
      }

      attempt <- attempt + 1

    }

    kept <- c(unlist(unname(sample_fossil)), tips)
    fossil_matrix <- matrix[rownames(matrix) %in% kept, ]

    # Discrete needs characters returned for ace, continuous needs numeric
    if (type == "discrete") {
      fossil_matrix <- apply(fossil_matrix, c(1, 2), as.character)
    } 

    if(type == "continuous") {
      fossil_matrix <- apply(fossil_matrix, c(1, 2), as.numeric)
      # fossil_matrix <- as.data.frame(fossil_matrix)
    }

    pruned <- bind.nodes.tips(tree,kept)
    rownames(fossil_matrix) <- ifelse(grepl("^n", rownames(fossil_matrix)), paste0("f_" , rownames(fossil_matrix)), rownames(fossil_matrix)) ## rename matrix as well
    # pruned <- keep.tip(tree, kept)
    return(list(matrix = fossil_matrix, tree = pruned))
  }

  # Check if inputs are lists or single objects
  if (is.list(trees) && is.list(matrices)) {
    # Use Map for lists of trees and matrices
    fossil_matrices <- Map(function(tree, matrix) process.fossil.alt(tree, matrix, type, seed), trees, matrices)
  } else {
    # Process a single tree and matrix
    fossil_matrices <- process.fossil.alt(trees, matrices, type, seed)
  }

  return(fossil_matrices)
}


bind.nodes.tips <- function(tree, kept) {
  kept_tips  <- kept[kept %in% tree$tip.label]
  kept_nodes <- kept[kept %in% tree$node.label]

  tree_with_fossils <- tree

  for (n_label in kept_nodes) {
    

    current_node_idx <- which(tree_with_fossils$node.label == n_label) + Ntip(tree_with_fossils)
    
    tree_with_fossils <- bind.tip(tree_with_fossils, 
                                  tip.label = paste0("f_",n_label), ## give new name (f_x) to separate from node label
                                  where = current_node_idx, 
                                  edge.length = 1e-6)
    }

  kept_nodes <- paste0("f_", kept_nodes) ## add "f_" logic to kept vector
  kept <- c(kept_nodes, kept_tips)
  pruned_tree <- keep.tip(tree_with_fossils, kept, collapse.singles = TRUE)
  pruned_tree <- multi2di(pruned_tree)
  pruned_tree$edge.length[pruned_tree$edge.length == 0] <- 1e-6
  return(pruned_tree)
}
