## function to create trees of only extant species
remove.fossil <- function(trees, matrices) {
  living_matrices <- Map(function(tree, matrix) {
    ages <- tree.age(tree) # get tip ages
    tips <- ages$element[ages$ages == 0] # find only living species
    living_matrix <- matrix[rownames(matrix) %in% tips, ]
    living_matrix <- apply(living_matrix, c(1, 2), as.character)
 # keep only living species in matrix
    living_tree <- keep.tip(tree, tips) # keep only living species in tree
    return(list(matrix = living_matrix, tree = living_tree)) # return as list
  }, trees, matrices)
  
  return(living_matrices)
}

## function to remove fossils at varying preservation levels
fossil.pres <- function(trees, matrices, preservation = c(0.05, 0.15, 0.5)){
  if (!preservation %in% c(0.05, 0.15, 0.5)) {
       stop("Invalid preservation value. Must be one of 0.05, 0.15, or 0.5.")
   }
   fossil_matrices <- Map(function(tree, matrix){
    ages <- tree.age(tree)
    tips <- ages$element[ages$ages == 0] # keep living species
    fossils <- ages$element[ages$ages > 0 & grepl("^t", ages$element)]  # find fossils
    if(preservation == 0.5){
        sample <- sample(fossils, size = length(fossils)/2)
        } ## sample fossils by a specified preservation level
    if(preservation == 0.15){
        sample <- sample(fossils, size = ceiling(length(fossils) * 0.15))
        }
    if(preservation == 0.05){
        sample <- sample(fossils, size = ceiling(length(fossils) * 0.05))
        }
    kept <- c(sample, tips)
    fossil_matrix <- matrix[rownames(matrix) %in% kept, ] # keep taxa that are fossil preserved
    fossil_matrix <- apply(fossil_matrix, c(1, 2), as.character) # make character for ace
    pruned  <- keep.tip(tree, kept) # prune tree
    return(list(matrix = fossil_matrix, tree = pruned))
   }, trees, matrices)
   return(fossil_matrices)
}