rate_names <- c("slow", "med", "fast")

# Load sample data (keep as is)
sum_var_change_sample <- setNames(lapply(rate_names, function(rate) {
    # Load all replicates first
    all_reps <- lapply(1:100, function(rep) {
            file_path <- sprintf("../Data/cluster/randomExtinction/disp/8092287_sample_sum_var_change_%03d.rds", rep)
            data <- readRDS(file_path)
            return(data[[rate]])
    })
}), rate_names
)

sum_var_change_no_ace <- setNames(lapply(rate_names, function(rate) {
    # Load all replicates first
    all_reps <- lapply(1:100, function(rep) {
            file_path <- sprintf("../Data/cluster/randomExtinction/disp/8092287_no_ace_sum_var_change_%03d.rds", rep)
            data <- readRDS(file_path)
            return(data[[rate]])
    })
}), rate_names
)



# Load true data (keep as is)
sum_var_change_true <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        tryCatch({
            file_path <- sprintf("../Data/cluster/randomExtinction/disp/8092287_true_sum_var_change_%03d.rds", rep)
            data <- readRDS(file_path)
            return(data[[rate]])
        }, warning = function(w) {
            NULL
        })
    })
}), rate_names)




# Calculate median on filtered sample data
sum_var_change_sample_median <- lapply(sum_var_change_sample, function(rate_data) {
    lapply(rate_data, function(replicate_data) {
        # Apply median directly to each fossil level
        result <- lapply(replicate_data, function(fossil_level_values) {
            values <- unlist(fossil_level_values, use.names = FALSE)
            return(median(values, na.rm = TRUE))
        })
        return(result)  # This ensures we return the fossil level list directly
    })
})

clean_sample <- lapply(sum_var_change_sample_median, lapply, function(x){
    x[!sapply(x, is.null)]  # Keep only non-NULL elements
})

clean_no_ace <- lapply(sum_var_change_no_ace, lapply, function(x){
    x[!sapply(x, is.null)]  # Keep only non-NULL elements
})

# Use filtered true data
# sum_var_change_true <- sum_var_change_true_filtered

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
matched_sample <- create.matched.data(sum_var_change_true, clean_sample)
matched_no_ace <- create.matched.data(sum_var_change_true, clean_no_ace)

all_ace <- c(matched_sample$estimated, matched_no_ace$estimated)



library(viridis)
fossil_colors <- viridis(5, option = "plasma", direction = -1)
names(fossil_colors) <- c("all", "fossil_high", "fossil_med", "fossil_low", "living")

plot(matched_sample$true, matched_sample$estimated,
    xlab = "True Disparity Change",
    ylab = "Est Disparity Change",
    main = "Disparity Change: No majority Rule vs True",
    pch = 16, 
    col = fossil_colors[matched_sample$fossil_level],
    cex = 1.7, 
    ylim = range(all_ace)
)

# Add 1:1 reference line
abline(0, 1, col = "red", lty = 2, lwd = 2)

plot(matched_no_ace$true, matched_no_ace$estimated,
    xlab = "True Disparity Change",
    ylab = "Estimated Disparity Change",
    main = "Disparity Change: No ace vs True",
    pch = 16, 
    col = fossil_colors[matched_no_ace$fossil_level],
    cex = 1.7,
    ylim = range(all_ace)
)

abline(0, 1, col = "red", lty = 2, lwd = 2)


