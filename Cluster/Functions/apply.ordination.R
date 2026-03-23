# for pre-ordination ace
apply.ordination <- function(distance_matrices) {
    # Apply function over each topology 
    ordination_results <- lapply(distance_matrices, function(topology) {
        #Loop over each distance matrix
        lapply(topology, function(distance_matrix) {
            # PCoA
            cmdscale(
                distance_matrix, 
                k = ncol(distance_matrix) - 2, 
                add = TRUE  # add caillez correction
            )
        })
    })
    return(ordination_results)
}


## for post-ordination ace
apply.ordination.two <- function(distance_matrices) {
    # Apply function over each topology 
    ordination_results <- lapply(distance_matrices, function(metric) {
            # PCoA
            cmdscale(
                metric, 
                k = ncol(metric) - 2, 
                add = TRUE  # add caillez correction
            )
        })
    return(ordination_results)
    }




pca.func <- function(matrices) {
    pca_results <- lapply(matrices, function(matrix) {
            # Perform Principal Coordinates Analysis 
          pca_output <-   prcomp(matrix, scale. = TRUE)
        pca_output$x    
})
    # Return the list of ordination results
    return(pca_results)
}


ord.func <- function (distances){
    ordination_matrices <- lapply(distances, function(distance_matrix){
    ordinated_matrix <- cmdscale(distance_matrix, k = ncol(distance_matrix) - 2, add = TRUE)
    traitspace <- ordinated_matrix$points
    })

    return(ordination_matrices)
    }