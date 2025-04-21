## function to create trees of only extant species
remove.fossil <- function(trees, matrices, type = c("discrete", "continuous")) {
  living_matrices <- Map(function(tree, matrix) {
    ages <- tree.age(tree) # get tip ages
    tips <- ages$element[ages$ages == 0]# find only living species
    # root <- tree$node.label[1] 
    living_matrix <- matrix[rownames(matrix) %in% c(tips, "t1"), ]
    if (type == "discrete"){
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
fossil.pres <- function(trees, matrices, preservation = c(0.05, 0.15, 0.5), type = c("discrete", "continuous")){
  if (!preservation %in% c(0.05, 0.15, 0.5)) {
       stop("Invalid preservation value. Must be one of 0.05, 0.15, or 0.5.")
   }
   fossil_matrices <- Map(function(tree, matrix){
    ages <- tree.age(tree)
    tips <- ages$element[ages$ages == 0]     # keep living species

    ## create time bins
    range <- range(ages$ages)
    bins <- seq(range[1], range[2], length.out = 4)
    bin_ids <- cut(ages$ages, breaks = bins, labels = FALSE, include.lowest = TRUE)
    ages$bin <- bin_ids
    # root <- tree$node.label[1]

    fossilb1 <- ages$element[ages$ages > 0 & grepl("^t", ages$element) & ages$bin == 1]
    fossilb2 <- ages$element[ages$ages > 0 & grepl("^t", ages$element) & ages$bin == 2]
    fossilb3 <- ages$element[ages$ages > 0 & grepl("^t", ages$element) & ages$bin == 3]

    fossils <- list(bin1 = fossilb1, bin2 = fossilb2, bin3 = fossilb3)  # find fossils

    sample <- lapply(fossils, function(bin) {
      sample(bin, size = ceiling(length(bin) * preservation))
    })

    kept <- c(unlist(unname(sample)), tips, "t1")
    fossil_matrix <- matrix[rownames(matrix) %in% kept, ]
    if (type == "discrete") { # keep taxa that are fossil preserved
    fossil_matrix <- apply(fossil_matrix, c(1, 2), as.character)
    } # make character for discrete ace
    pruned  <- keep.tip(tree, kept) # prune tree
    return(list(matrix = fossil_matrix, tree = pruned))
   }, trees, matrices)
   return(fossil_matrices)
}
