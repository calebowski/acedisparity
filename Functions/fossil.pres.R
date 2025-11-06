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
    tree <- set.root.time(tree)
    old_root <- tree$root.time
    living_tree$root.time <- old_root
    
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

## function to remove fossils at varying preservation levels
## idea - split time into time bins and vary preservation level per bin
fossil.pres.stages <- function(trees, matrices, preservation = c(0.05, 0.15, 0.5, 1.0), type = c("discrete", "continuous")) {
  # if (!preservation %in% c(0.05, 0.15, 0.5, 1.0)) {
  #   stop("Invalid preservation value. Must be one of 1.0, 0.05, 0.15, or 0.5.")
  # } ## commented this out so that preservation level can be any value between >0 - 1.

  # Helper function to process a single tree and matrix
  process.fossil.stages <- function(tree, matrix, type) {
    ages <- tree.age(tree)
    tips <- ages$element[ages$ages == 0]  # Keep living species

    ## Create time bins
    range <- range(ages$ages) 
    bins <- seq(range[1], range[2], length.out = 4)
    bin_ids <- cut(ages$ages, breaks = bins, labels = FALSE, include.lowest = TRUE)
    ages$bin <- bin_ids

    fossilb1 <- ages$element[ages$ages > 0 & grepl("^t", ages$element) & ages$bin == 1]
    fossilb2 <- ages$element[ages$ages > 0 & grepl("^t", ages$element) & ages$bin == 2]
    fossilb3 <- ages$element[ages$ages > 0 & grepl("^t", ages$element) & ages$bin == 3]

    fossils <- list(bin1 = fossilb1, bin2 = fossilb2, bin3 = fossilb3)  # Find fossils

    sample <- lapply(fossils, function(bin) {
      sample(bin, size = ceiling(length(bin) * preservation))
    })

    kept <- c(unlist(unname(sample)), tips)
    fossil_matrix <- matrix[rownames(matrix) %in% kept, ]

    # Discrete needs characters returned for ace, continuous needs numeric
    if (type == "discrete") {
      fossil_matrix <- apply(fossil_matrix, c(1, 2), as.character)
    } 

    if(type == "continuous") {
      fossil_matrix <- apply(fossil_matrix, c(1, 2), as.numeric)
      # fossil_matrix <- as.data.frame(fossil_matrix)
    }


    pruned <- keep.tip(tree, kept) # Rescale tree so that tree height is maintained
    tree <- set.root.time(tree)
    old_root <- tree$root.time
    pruned$root.time <- old_root
    return(list(matrix = fossil_matrix, tree = pruned))
  }

  # Check if inputs are lists or single objects
  if (is.list(trees) && is.list(matrices)) {
    # Use Map for lists of trees and matrices
    fossil_matrices <- Map(function(tree, matrix) process.fossil(tree, matrix, type), trees, matrices)
  } else {
    # Process a single tree and matrix
    fossil_matrices <- process.fossil(trees, matrices, type)
  }

  return(fossil_matrices)
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


    pruned <- keep.tip(tree, kept) # Rescale tree so that tree height is maintained
    tree <- set.root.time(tree)
    old_root <- tree$root.time
    pruned$root.time <- old_root
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


# With Bernoulli sampling, the probability that any given species is preserved doesn’t depend on how many species existed at the same time.
# That avoids the bias that trees with low diversity won’t automatically get overrepresented.