find.extinction.time <- function(extinction_ages, exclude_zero = TRUE, tolerance = 0.1) {
  ages <- extinction_ages$ages
  
  # Remove ages at or very close to 0 (terminal tips)
  if(exclude_zero) {
    ages <- ages[ages > tolerance]
  }
  
  # Round ages to handle slight numerical differences
  # rounded_ages <- round(ages, 1)
  
  # Count frequency of each age
  age_counts <- table(ages)
  
  # Find the most common age
  if(length(age_counts) > 0) {
    extinction_time <- as.numeric(names(age_counts)[which.max(age_counts)])
    extinction_count <- max(age_counts)
    

    
    return(extinction_time)
  } else {
    cat("No mass extinction found\n")
    return(list(time = NA, count = 0))
  }
}



tip.ages <- function(tree) {
  ages <- tree.age(tree)
  tips <- grepl("^t", ages$elements)
  tip_ages <- ages[tips, ]  
  return(tip_ages)
}