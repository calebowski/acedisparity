## read in the disp diffs


# Read and rearrange in one step
rate_names <- c("slow", "med", "fast")  

disp_diff_strict <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/dispdiffs/rawoutput/diff_disp_strict_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

disp_diff_sample <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/dispdiffs/rawoutput/diff_disp_sample_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

disp_diff_rel <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/dispdiffs/rawoutput/diff_disp_rel_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)

disp_diff_no_ace <- setNames(lapply(rate_names, function(rate) {
    lapply(1:100, function(rep) {
        file_path <- sprintf("../Data/cluster/discrete/dispdiffs/rawoutput/diff_disp_no_ace_%03d.rds", rep)
        data <- readRDS(file_path)
        return(data[[rate]])
    })
}), rate_names)



saveRDS(disp_diff_strict, "../Data/cluster/discrete/dispdiffs/disp_diff_strict.rds")
saveRDS(disp_diff_sample, "../Data/cluster/discrete/dispdiffs/disp_diff_sample.rds")
saveRDS(disp_diff_rel, "../Data/cluster/discrete/dispdiffs/disp_diff_rel.rds")
saveRDS(disp_diff_no_ace, "../Data/cluster/discrete/dispdiffs/disp_diff_no_ace.rds")
