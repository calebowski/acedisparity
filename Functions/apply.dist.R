apply.distance.metrics <- function(anc_states_list, metrics = c("mord", "euclidean", "manhattan", "hamming")) {
    ## Loop over each topology 
    distances <- lapply(anc_states_list, function(anc_states) {
        # Call char.diff to apply each metric to the matrix of ancestral states
        metric_distances <- lapply(metrics, function(metric) {
            char.diff(matrix = anc_states, method = metric, by.col = FALSE)
        })
        names(metric_distances) <- metrics  # Name each matrix by metric
        return(metric_distances)
    })
    names(distances) <- names(anc_states_list)  # Name by topology
    return(distances)
}

apply.distance.metrics.two <- function(matrix, metrics = c("mord", "euclidean", "manhattan", "hamming")) {
        metric_distances <- lapply(metrics, function(metric) {
            char.diff(matrix = matrix, method = metric, by.col = FALSE)
        })
        names(metric_distances) <- metrics  # Name each matrix by metric
        return(metric_distances)
    }



dist.func <- function (matrix_list){
    distance_matrices <- lapply(matrix_list, char.diff, method = "mord", by.col = FALSE)
    return(distance_matrices)
}