rate_names <- c("slow", "med", "fast")

# Load sample data (keep as is)
sum_var_change_sample <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        tryCatch({
            file_path <- sprintf("../Data/cluster/randomExtinction/disp/rand_ext_change_sum_var_sample_%03d.rds", rep)
            data <- readRDS(file_path)
            return(data[[rate]])
        }, warning = function(w) {
            NULL
        })
    })
}), rate_names)

# Load true data (keep as is)
sum_var_change_true_raw <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        tryCatch({
            file_path <- sprintf("../Data/cluster/randomExtinction/disp/rand_ext_change_sum_var_true_%03d.rds", rep)
            data <- readRDS(file_path)
            return(data[[rate]])
        }, warning = function(w) {
            NULL
        })
    })
}), rate_names)

# FILTER BOTH DATASETS TO MATCHING NON-NULL REPLICATES
sum_var_change_sample_filtered <- lapply(rate_names, function(rate) {
    # Find which replicates are complete in BOTH datasets
    sample_complete <- !sapply(sum_var_change_sample[[rate]], is.null)
    true_complete <- !sapply(sum_var_change_true_raw[[rate]], is.null)
    both_complete <- sample_complete & true_complete
    
    # Keep only replicates that are complete in both
    return(sum_var_change_sample[[rate]][both_complete])
})
names(sum_var_change_sample_filtered) <- rate_names

sum_var_change_true_filtered <- lapply(rate_names, function(rate) {
    # Find which replicates are complete in BOTH datasets
    sample_complete <- !sapply(sum_var_change_sample[[rate]], is.null)
    true_complete <- !sapply(sum_var_change_true_raw[[rate]], is.null)
    both_complete <- sample_complete & true_complete
    
    # Keep only replicates that are complete in both
    return(sum_var_change_true_raw[[rate]][both_complete])
})
names(sum_var_change_true_filtered) <- rate_names

# Calculate median on filtered sample data
sum_var_change_sample_median <- lapply(sum_var_change_sample_filtered, function(rate_data) {
    lapply(rate_data, function(replicate_data) {
        # Apply median directly to each fossil level
        result <- lapply(replicate_data, function(fossil_level_values) {
            values <- unlist(fossil_level_values, use.names = FALSE)
            return(median(values, na.rm = TRUE))
        })
        return(result)  # This ensures we return the fossil level list directly
    })
})

# Use filtered true data
sum_var_change_true <- sum_var_change_true_filtered

# Your function should now work perfectly
create.matched.data <- function(true, estimated) {
    true_matched <- c()
    est_matched <- c()
    fossil_levels <- c()

    for(rate_name in names(estimated)) {
        for(rep in seq_along(estimated[[rate_name]])) {
            true_val <- true[[rate_name]][[rep]]
            est_rep_val <- estimated[[rate_name]][[rep]]
            fossil_names <- names(est_rep_val)
            n_levels <- length(est_rep_val)

            true_matched <- c(true_matched, rep(true_val, n_levels))
            est_matched <- c(est_matched, unlist(est_rep_val, use.names = FALSE))
            fossil_levels <- c(fossil_levels, fossil_names)
        }
    }
    return(list(true = true_matched, estimated = est_matched, fossil_level = fossil_levels))
}

# Test the match
matched_sample <- create.matched.data(sum_var_change_true, sum_var_change_sample_median)

cat("True length:", length(matched_sample$true), "\n")
cat("Estimated length:", length(matched_sample$estimated), "\n")
cat("Should be equal now!\n")

# Plot
library(viridis)
fossil_colors <- viridis(5, option = "plasma", direction = -1)
names(fossil_colors) <- c("all", "fossil_high", "fossil_med", "fossil_low", "living")

plot(matched_sample$true, matched_sample$estimated,
     xlab = "True Disparity Change",
     ylab = "Estimated Disparity Change",
     main = "Disparity Change: No majority Rule vs True",
     pch = 16, 
     col = fossil_colors[matched_sample$fossil_level],
     cex = 1.7
)

# Add 1:1 reference line
abline(0, 1, col = "red", lty = 2, lwd = 2)