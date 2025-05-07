## function to create trees of only extant species
remove.fossil <- function(trees, matrices, type = c("discrete", "continuous")) {
  living_matrices <- Map(function(tree, matrix) {
    ages <- tree.age(tree) # get tip ages
    tips <- ages$element[ages$ages == 0]# find only living species
    # root <- tree$node.label[1] 
    living_matrix <- matrix[rownames(matrix) %in% c(tips, "t1"), ]
    living_matrix["t1", ] <- "?"
    if (type == "discrete") {
    living_matrix <- apply(living_matrix, c(1, 2), as.character)
    }
 # keep only living species in matrix
    living_tree <- keep.tip(tree, c(tips, "t1")) # keep only living species in tree
    return(list(matrix = living_matrix, tree = living_tree)) # return as list
  }, trees, matrices)
  
  return(living_matrices)
}

## function to remove fossils at varying preservation levels
## idea - split time into time bins and vary preservation level per bin
fossil.pres <- function(trees, matrices, preservation = c(0.05, 0.15, 0.5, 1.0), type = c("discrete", "continuous")) {
  if (!preservation %in% c(0.05, 0.15, 0.5, 1.0)) {
    stop("Invalid preservation value. Must be one of 1.0, 0.05, 0.15, or 0.5.")
  }

  # Helper function to process a single tree and matrix
  process.fossil <- function(tree, matrix, type) {
    # print(paste("Processing with type:", type))
    
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

    kept <- c(unlist(unname(sample)), tips, "t1")
    fossil_matrix <- matrix[rownames(matrix) %in% kept, ]

    # if (preservation %in% c(0.05, 0.15)) {
    #   fossil_matrix["t1", ] <- "?"
    # }

    # Ensure type is respected
    if (type == "discrete") {  
      fossil_matrix <- apply(fossil_matrix, c(1, 2), as.character)
    } 

    if(type == "continuous") {
      fossil_matrix <- apply(fossil_matrix, c(1, 2), as.numeric)
      # fossil_matrix <- as.data.frame(fossil_matrix)
    }

    if (preservation %in% c(0.05, 0.15) && type == "discrete") {
      fossil_matrix["t1", ] <- "?"
    }

    pruned <- keep.tip(tree, kept)  # Prune tree
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