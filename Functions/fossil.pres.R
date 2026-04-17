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