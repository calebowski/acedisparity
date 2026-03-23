handle.missing.data <- function(matrix, missing_value = "?") {
#   all_states <- sort(unique(as.vector(matrix[matrix != missing_value])))

    output <- lapply(1:ncol(matrix), function(y) {
    characters <- matrix[, y]
    
    chars_scored <- characters[which(!characters == missing_value)]
    chars_missing <- characters[which(characters == missing_value)]
    
    # Use full set of states to build matrix
    char_states <- sort(unique(chars_scored))
    # Build scored matrix
    chars_scored_matrix <- to.matrix(chars_scored, seq = char_states)
    
    # Build missing matrix with same columns
    chars_missing_matrix <- matrix(1 / length(char_states),
                                   nrow = length(chars_missing),
                                   ncol = length(char_states))
    rownames(chars_missing_matrix) <- names(chars_missing)
    colnames(chars_missing_matrix) <- colnames(chars_scored_matrix)
    
    # Combine
    combined_matrix <- rbind(chars_scored_matrix, chars_missing_matrix)
    
    # Return matrix with consistent row order
    combined_matrix[rownames(matrix), , drop = FALSE]
  })
  
  return(output)
}