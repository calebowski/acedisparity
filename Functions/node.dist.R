# Distance metric function with lower penalty for uncertainty
lenient.distance <- function(true_state, ace_states) {
  # If the ace state contains true state, return 0
  if (true_state %in% ace_states) {
    return(0)
  }

  # Otherwise, the true state and estimated state are different (i.e., distance is 1)
  return(1)
}

# Distance metric with strict penalty against uncertainty
strict.distance <- function(true_state, ace_states) {
  # If the state is exact match, distance is 0
  if (true_state %in% ace_states && length(ace_states) == 1) {
    return(0)
  }
  
  # anything else is returned as 1
  return(1)
}


node.dist <- function(matrices, distance = c("strict", "uncertain")) {
  true <- matrices$true
  ace <- matrices$ace

  if (inherits(ace, "list")) {
    return(lapply(ace, function(ace_matrix) {
      node.dist(list(true = true, ace = ace_matrix), distance = distance)
    }))
  }

  n_rows <- nrow(ace)
  n_cols <- ncol(ace)


  parsed_ace <- vector("list", n_rows)
  for (i in 1:n_rows) {
    parsed_ace[[i]] <- vector("list", n_cols)
    for (j in 1:n_cols) {
      cell <- ace[i, j]
      if (is.na(cell) || cell == "") {
        parsed_ace[[i]][[j]] <- NA
      } else {
        parsed_ace[[i]][[j]] <- as.numeric(unlist(strsplit(cell, "/")))
      }
    }
  }
  distance_matrix <- matrix(NA, nrow = n_rows, ncol = n_cols)

  if (distance == "lenient"){
  # Loop through nodes and characters
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      true_states <- as.numeric(true[i, j])
      ace_states <- parsed_ace[[i]][[j]]
      distance_matrix[i, j] <- lenient.distance(true_states, ace_states)
    }
  }
  }

  if (distance == "strict"){
    for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      true_states <- as.numeric(true[i, j])
      ace_states <- parsed_ace[[i]][[j]]
      distance_matrix[i, j] <- strict.distance(true_states, ace_states)
    }
  }
  }

  return(distance_matrix)

}