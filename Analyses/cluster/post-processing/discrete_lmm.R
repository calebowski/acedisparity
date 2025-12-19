library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest) 
# library(emmeans)
# library(multcomp)
# Define job IDs for each tree size
job_ids <- list(
  "50t" = "8558401",
  "100t" = "8561556",
  "150t" = "8561712"
)

tree_sizes <- c("50t", "100t", "150t")
methods <- c("pre_ord_point_tiebreaker", "pre_ord_sample", "post_ord_point", 
             "post_ord_sample", "no_ace")


args <- commandArgs(trailingOnly = TRUE)

if(args[1] == "raw") {
    cat("Starting Raw LMM...\n")
    # Load all data with different job IDs
    raw_results <- list()
    for (size in tree_sizes) {
    raw_results[[size]] <- list()
    job_id <- job_ids[[size]]  #  Get job ID for this tree size
    
    for (method in methods) {
        raw_results[[size]][[method]] <- list()
        for(i in 1:100) {
        file_path <- file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "discrete", "crown", size,
                                "disparity", "raw",
                                sprintf("%s_%s_%03d.rds", job_id, method, i))
        if(file.exists(file_path)) {
            raw_results[[size]][[method]][[i]] <- readRDS(file_path)
        } else {
            warning("Missing file: ", file_path)
            raw_results[[size]][[method]][[i]] <- NULL
        }
        }
    }
    }

    metric_names <- names(raw_results[["50t"]][[1]][[1]])  # Get from first available


    # Transpose structure: metric -> tree_size -> method -> replicate
    raw_results_relisted <- setNames(lapply(metric_names, function(metric) {
    setNames(lapply(tree_sizes, function(size) {
        setNames(lapply(methods, function(method) {
        is_sample_method <- grepl("sample", method)
        
        lapply(raw_results[[size]][[method]], function(rep_data) {
            metric_data <- rep_data[[metric]]
            
            if(is_sample_method) {
            lapply(metric_data, function(rate_data) {
                lapply(rate_data, function(fossil_data) {
                if(is.list(fossil_data) && !is.data.frame(fossil_data)) {
                    unlist(fossil_data)
                } else {
                    fossil_data
                }
                })
            })
            } else {
            metric_data
            }
        })
        }), methods)
    }), tree_sizes)
    }), metric_names)

    # Build dataframe with tree_size
    raw_results_df <- do.call(rbind, lapply(names(raw_results_relisted), function(metric_name) {
    metric_data <- raw_results_relisted[[metric_name]]
    
    do.call(rbind, lapply(names(metric_data), function(size_name) {
        size_data <- metric_data[[size_name]]
        
        do.call(rbind, lapply(names(size_data), function(method_name) {
        method_data <- size_data[[method_name]]
        
        do.call(rbind, lapply(seq_along(method_data), function(rep_idx) {
            rep_data <- method_data[[rep_idx]]
            
            if(is.null(rep_data)) return(NULL)
            
            do.call(rbind, lapply(names(rep_data), function(rate_name) {
            rate_data <- rep_data[[rate_name]]
            
            if(is.null(rate_data)) return(NULL)
            
            data.frame(
                replicate = rep_idx,
                all = I(list(rate_data$all)),
                high = I(list(rate_data$fossil_high)),
                mid = I(list(rate_data$fossil_med)),
                low = I(list(rate_data$fossil_low)),
                living = I(list(rate_data$living)),
                rate = rate_name,
                method = method_name,
                metric = metric_name,
                tree_size = size_name,  #  Add tree size
                stringsAsFactors = FALSE
            )
            }))
        }))
        }))
    }))
    }))

    cat("raw_results_df created...\n")

    rm(raw_results)
    gc()

    # Pivot to long format
    raw_results_df_long <- raw_results_df %>%
    pivot_longer(
        cols = c("living", "low", "mid", "high", "all"),
        names_to = "preservation_level",
        values_to = "error"
    ) %>%
    unnest(error)#

    rm(raw_results_df)
    gc()

    # Factor variables
    raw_results_df_long$tree_size <- factor(
    raw_results_df_long$tree_size,
    levels = c("50t", "100t", "150t"),
    labels = c("50 taxa", "100 taxa", "150 taxa")
    )

    raw_results_df_long$method <- factor(
    raw_results_df_long$method,
    levels = c("pre_ord_point_tiebreaker", "pre_ord_sample", 
                "post_ord_point", "post_ord_sample", "no_ace"),
    labels = c("Pre-ord ASE\n(point + tiebreaker)",
                "Pre-ord ASE\n(dist)", "Post-ord ASE\n(point)", 
                "Post-ord ASE\n(dist)", "No ASE")
    )

    raw_results_df_long$preservation_level <- factor(
    raw_results_df_long$preservation_level,
    levels = c("living", "low", "mid", "high", "all"),
    labels = c("0%", "5%", "15%", "50%", "100%")
    )

    raw_results_df_long$rate <- factor(
    raw_results_df_long$rate,
    levels = c("slow", "med", "fast"),
    labels = c("Slow", "Medium", "Fast")
    )

    raw_results_df_long$metric <- factor(
    raw_results_df_long$metric,
    levels = c("sum_var", "sum_quant", "pairwise"),
    labels = c("Sum of Variances", "Sum of Quantiles", "Mean Pairwise Distance")
    )

    cat("raw_results_df_long created...\n")


    raw_results_df_long$abs_error <- abs(raw_results_df_long$error)
    raw_results_df_long$log_abs_error <- log(raw_results_df_long$abs_error + 0.001)


    lmm_model <- lmer(log_abs_error ~ rate * method * preservation_level * metric + 
                        (1|replicate) + (1|tree_size), 
                        data = raw_results_df_long,
                        REML = TRUE)


    saveRDS(lmm_model, file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "discrete", "crown", "lmm", "lmm_discrete_raw.rds"))

    filtered_model <- lmerTest::step(lmm_model) ## backwards step selection to find best model

    final_formula <- lmerTest::get_model(filtered_model)

    lmm_model_final <- lmer(final_formula, data = results_df_long_filtered, REML = FALSE)

    saveRDS(lmm_model_final, file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "discrete", "crown","lmm", "lmm_discrete_raw_final.rds"))

    cat("Finished Raw LMM...\n")
}


if(args[1] == "rmaxes") {

    cat("Starting rmaxes LMM...\n")

    # Load all data with different job IDs
    rmaxes_results <- list()
    for (size in tree_sizes) {
    rmaxes_results[[size]] <- list()
    job_id <- job_ids[[size]]  #  Get job ID for this tree size
    
    for (method in methods) {
        rmaxes_results[[size]][[method]] <- list()
        for(i in 1:100) {
        file_path <- file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "discrete", "crown", size,
                            "disparity", "rm_axes",
                                sprintf("%s_%s_%03d.rds", job_id, method, i))
        if(file.exists(file_path)) {
            rmaxes_results[[size]][[method]][[i]] <- readRDS(file_path)
        } else {
            warning("Missing file: ", file_path)
            rmaxes_results[[size]][[method]][[i]] <- NULL
        }
        }
    }
    }

    metric_names <- names(rmaxes_results[["50t"]][[1]][[1]])  # Get from first available


    # Transpose structure: metric -> tree_size -> method -> replicate
    rmaxes_results_relisted <- setNames(lapply(metric_names, function(metric) {
    setNames(lapply(tree_sizes, function(size) {
        setNames(lapply(methods, function(method) {
        is_sample_method <- grepl("sample", method)
        
        lapply(rmaxes_results[[size]][[method]], function(rep_data) {
            metric_data <- rep_data[[metric]]
            
            if(is_sample_method) {
            lapply(metric_data, function(rate_data) {
                lapply(rate_data, function(fossil_data) {
                if(is.list(fossil_data) && !is.data.frame(fossil_data)) {
                    unlist(fossil_data)
                } else {
                    fossil_data
                }
                })
            })
            } else {
            metric_data
            }
        })
        }), methods)
    }), tree_sizes)
    }), metric_names)

    # Build dataframe with tree_size
    rmaxes_results_df <- do.call(rbind, lapply(names(rmaxes_results_relisted), function(metric_name) {
    metric_data <- rmaxes_results_relisted[[metric_name]]
    
    do.call(rbind, lapply(names(metric_data), function(size_name) {
        size_data <- metric_data[[size_name]]
        
        do.call(rbind, lapply(names(size_data), function(method_name) {
        method_data <- size_data[[method_name]]
        
        do.call(rbind, lapply(seq_along(method_data), function(rep_idx) {
            rep_data <- method_data[[rep_idx]]
            
            if(is.null(rep_data)) return(NULL)
            
            do.call(rbind, lapply(names(rep_data), function(rate_name) {
            rate_data <- rep_data[[rate_name]]
            
            if(is.null(rate_data)) return(NULL)
            
            data.frame(
                replicate = rep_idx,
                all = I(list(rate_data$all)),
                high = I(list(rate_data$fossil_high)),
                mid = I(list(rate_data$fossil_med)),
                low = I(list(rate_data$fossil_low)),
                living = I(list(rate_data$living)),
                rate = rate_name,
                method = method_name,
                metric = metric_name,
                tree_size = size_name,  #  Add tree size
                stringsAsFactors = FALSE
            )
            }))
        }))
        }))
    }))
    }))

    # Pivot to long format
    rmaxes_results_df_long <- rmaxes_results_df %>%
    pivot_longer(
        cols = c("living", "low", "mid", "high", "all"),
        names_to = "preservation_level",
        values_to = "error"
    ) %>%
    unnest(error)

    # Factor variables
    rmaxes_results_df_long$tree_size <- factor(
    rmaxes_results_df_long$tree_size,
    levels = c("50t", "100t", "150t"),
    labels = c("50 taxa", "100 taxa", "150 taxa")
    )

    rmaxes_results_df_long$method <- factor(
    rmaxes_results_df_long$method,
    levels = c("pre_ord_point_tiebreaker", "pre_ord_sample", 
                "post_ord_point", "post_ord_sample", "no_ace"),
    labels = c("Pre-ord ASE\n(point + tiebreaker)",
                "Pre-ord ASE\n(dist)", "Post-ord ASE\n(point)", 
                "Post-ord ASE\n(dist)", "No ASE")
    )

    rmaxes_results_df_long$preservation_level <- factor(
    rmaxes_results_df_long$preservation_level,
    levels = c("living", "low", "mid", "high", "all"),
    labels = c("0%", "5%", "15%", "50%", "100%")
    )

    rmaxes_results_df_long$rate <- factor(
    rmaxes_results_df_long$rate,
    levels = c("slow", "med", "fast"),
    labels = c("Slow", "Medium", "Fast")
    )

    rmaxes_results_df_long$metric <- factor(
    rmaxes_results_df_long$metric,
    levels = c("sum_var", "sum_quant", "pairwise"),
    labels = c("Sum of Variances", "Sum of Quantiles", "Mean Pairwise Distance")
    )

    rmaxes_results_df_long$abs_error <- abs(rmaxes_results_df_long$error)
    rmaxes_results_df_long$log_abs_error <- log(rmaxes_results_df_long$abs_error + 0.001)


    lmm_model <- lmer(log_abs_error ~ rate * method * preservation_level * metric + 
                        (1|replicate) + (1|tree_size), 
                        data = rmaxes_results_df_long,
                        REML = FALSE)


    saveRDS(lmm_model, file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "discrete", "crown", "lmm", "lmm_discrete_rmaxes.rds"))

}

if(args[1] == "combined_raw_rmaxes") {
    cat("Starting combined raw and rmaxes LMM...\n")

    # Load all data with different job IDs
    raw_results <- list()
    for (size in tree_sizes) {
    raw_results[[size]] <- list()
    job_id <- job_ids[[size]]  #  Get job ID for this tree size
    
    for (method in methods) {
        raw_results[[size]][[method]] <- list()
        for(i in 1:100) {
        file_path <- file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "discrete", "crown", size,
                            "disparity", "raw",
                                sprintf("%s_%s_%03d.rds", job_id, method, i))
        if(file.exists(file_path)) {
            raw_results[[size]][[method]][[i]] <- readRDS(file_path)
        } else {
            warning("Missing file: ", file_path)
            raw_results[[size]][[method]][[i]] <- NULL
        }
        }
    }
    }

    metric_names <- names(raw_results[["50t"]][[1]][[1]])  # Get from first available


    # Transpose structure: metric -> tree_size -> method -> replicate
    raw_results_relisted <- setNames(lapply(metric_names, function(metric) {
    setNames(lapply(tree_sizes, function(size) {
        setNames(lapply(methods, function(method) {
        is_sample_method <- grepl("sample", method)
        
        lapply(raw_results[[size]][[method]], function(rep_data) {
            metric_data <- rep_data[[metric]]
            
            if(is_sample_method) {
            lapply(metric_data, function(rate_data) {
                lapply(rate_data, function(fossil_data) {
                if(is.list(fossil_data) && !is.data.frame(fossil_data)) {
                    unlist(fossil_data)
                } else {
                    fossil_data
                }
                })
            })
            } else {
            metric_data
            }
        })
        }), methods)
    }), tree_sizes)
    }), metric_names)

    rm(raw_results)
    gc()


    # Build dataframe with tree_size
    raw_results_df <- do.call(rbind, lapply(names(raw_results_relisted), function(metric_name) {
    metric_data <- raw_results_relisted[[metric_name]]
    
    do.call(rbind, lapply(names(metric_data), function(size_name) {
        size_data <- metric_data[[size_name]]
        
        do.call(rbind, lapply(names(size_data), function(method_name) {
        method_data <- size_data[[method_name]]
        
        do.call(rbind, lapply(seq_along(method_data), function(rep_idx) {
            rep_data <- method_data[[rep_idx]]
            
            if(is.null(rep_data)) return(NULL)
            
            do.call(rbind, lapply(names(rep_data), function(rate_name) {
            rate_data <- rep_data[[rate_name]]
            
            if(is.null(rate_data)) return(NULL)
            
            data.frame(
                replicate = rep_idx,
                all = I(list(rate_data$all)),
                high = I(list(rate_data$fossil_high)),
                mid = I(list(rate_data$fossil_med)),
                low = I(list(rate_data$fossil_low)),
                living = I(list(rate_data$living)),
                rate = rate_name,
                method = method_name,
                metric = metric_name,
                tree_size = size_name,  #  Add tree size
                stringsAsFactors = FALSE
            )
            }))
        }))
        }))
    }))
    }))

    rm(raw_results_relisted)
    gc()

    # Pivot to long format
    raw_results_df_long <- raw_results_df %>%
    pivot_longer(
        cols = c("living", "low", "mid", "high", "all"),
        names_to = "preservation_level",
        values_to = "error"
    ) %>%
    unnest(error)

    rm(raw_results_df)
    gc()

    # Factor variables
    raw_results_df_long$tree_size <- factor(
    raw_results_df_long$tree_size,
    levels = c("50t", "100t", "150t"),
    labels = c("50 taxa", "100 taxa", "150 taxa")
    )

    raw_results_df_long$method <- factor(
    raw_results_df_long$method,
    levels = c("pre_ord_point_tiebreaker", "pre_ord_sample", 
                "post_ord_point", "post_ord_sample", "no_ace"),
    labels = c("Pre-ord ASE\n(point + tiebreaker)",
                "Pre-ord ASE\n(dist)", "Post-ord ASE\n(point)", 
                "Post-ord ASE\n(dist)", "No ASE")
    )

    raw_results_df_long$preservation_level <- factor(
    raw_results_df_long$preservation_level,
    levels = c("living", "low", "mid", "high", "all"),
    labels = c("0%", "5%", "15%", "50%", "100%")
    )

    raw_results_df_long$rate <- factor(
    raw_results_df_long$rate,
    levels = c("slow", "med", "fast"),
    labels = c("Slow", "Medium", "Fast")
    )

    raw_results_df_long$metric <- factor(
    raw_results_df_long$metric,
    levels = c("sum_var", "sum_quant", "pairwise"),
    labels = c("Sum of Variances", "Sum of Quantiles", "Mean Pairwise Distance")
    )

    # Load all data with different job IDs
    rmaxes_results <- list()
    for (size in tree_sizes) {
    rmaxes_results[[size]] <- list()
    job_id <- job_ids[[size]]  #  Get job ID for this tree size
    
    for (method in methods) {
        rmaxes_results[[size]][[method]] <- list()
        for(i in 1:100) {
        file_path <- file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "discrete", "crown", size,
                            "disparity", "rm_axes",
                                sprintf("%s_%s_%03d.rds", job_id, method, i))
        if(file.exists(file_path)) {
            rmaxes_results[[size]][[method]][[i]] <- readRDS(file_path)
        } else {
            warning("Missing file: ", file_path)
            rmaxes_results[[size]][[method]][[i]] <- NULL
        }
        }
    }
    }

    metric_names <- names(rmaxes_results[["50t"]][[1]][[1]])  # Get from first available


    # Transpose structure: metric -> tree_size -> method -> replicate
    rmaxes_results_relisted <- setNames(lapply(metric_names, function(metric) {
    setNames(lapply(tree_sizes, function(size) {
        setNames(lapply(methods, function(method) {
        is_sample_method <- grepl("sample", method)
        
        lapply(rmaxes_results[[size]][[method]], function(rep_data) {
            metric_data <- rep_data[[metric]]
            
            if(is_sample_method) {
            lapply(metric_data, function(rate_data) {
                lapply(rate_data, function(fossil_data) {
                if(is.list(fossil_data) && !is.data.frame(fossil_data)) {
                    unlist(fossil_data)
                } else {
                    fossil_data
                }
                })
            })
            } else {
            metric_data
            }
        })
        }), methods)
    }), tree_sizes)
    }), metric_names)

    rm(rmaxes_results)
    gc()

    # Build dataframe with tree_size
    rmaxes_results_df <- do.call(rbind, lapply(names(rmaxes_results_relisted), function(metric_name) {
    metric_data <- rmaxes_results_relisted[[metric_name]]
    
    do.call(rbind, lapply(names(metric_data), function(size_name) {
        size_data <- metric_data[[size_name]]
        
        do.call(rbind, lapply(names(size_data), function(method_name) {
        method_data <- size_data[[method_name]]
        
        do.call(rbind, lapply(seq_along(method_data), function(rep_idx) {
            rep_data <- method_data[[rep_idx]]
            
            if(is.null(rep_data)) return(NULL)
            
            do.call(rbind, lapply(names(rep_data), function(rate_name) {
            rate_data <- rep_data[[rate_name]]
            
            if(is.null(rate_data)) return(NULL)
            
            data.frame(
                replicate = rep_idx,
                all = I(list(rate_data$all)),
                high = I(list(rate_data$fossil_high)),
                mid = I(list(rate_data$fossil_med)),
                low = I(list(rate_data$fossil_low)),
                living = I(list(rate_data$living)),
                rate = rate_name,
                method = method_name,
                metric = metric_name,
                tree_size = size_name,  #  Add tree size
                stringsAsFactors = FALSE
            )
            }))
        }))
        }))
    }))
    }))

    rm(rmaxes_results_relisted)
    gc()

    # Pivot to long format
    rmaxes_results_df_long <- rmaxes_results_df %>%
    pivot_longer(
        cols = c("living", "low", "mid", "high", "all"),
        names_to = "preservation_level",
        values_to = "error"
    ) %>%
    unnest(error)

    rm(rmaxes_results_df)
    gc()

    # Factor variables
    rmaxes_results_df_long$tree_size <- factor(
    rmaxes_results_df_long$tree_size,
    levels = c("50t", "100t", "150t"),
    labels = c("50 taxa", "100 taxa", "150 taxa")
    )

    rmaxes_results_df_long$method <- factor(
    rmaxes_results_df_long$method,
    levels = c("pre_ord_point_tiebreaker",  "pre_ord_sample", 
                "post_ord_point", "post_ord_sample", "no_ace"),
    labels = c("Pre-ord ASE\n(point + tiebreaker)",
                "Pre-ord ASE\n(dist)", "Post-ord ASE\n(point)", 
                "Post-ord ASE\n(dist)", "No ASE")
    )

    rmaxes_results_df_long$preservation_level <- factor(
    rmaxes_results_df_long$preservation_level,
    levels = c("living", "low", "mid", "high", "all"),
    labels = c("0%", "5%", "15%", "50%", "100%")
    )

    rmaxes_results_df_long$rate <- factor(
    rmaxes_results_df_long$rate,
    levels = c("slow", "med", "fast"),
    labels = c("Slow", "Medium", "Fast")
    )

    rmaxes_results_df_long$metric <- factor(
    rmaxes_results_df_long$metric,
    levels = c("sum_var", "sum_quant", "pairwise"),
    labels = c("Sum of Variances", "Sum of Quantiles", "Mean Pairwise Distance")
    )

    raw_results_df_long$axes_treatment <- "Raw"
    rmaxes_results_df_long$axes_treatment <- "Axes Removed"

    # Now combine
    combined_rm_raw_df <- rbind(raw_results_df_long, rmaxes_results_df_long)

    rm(raw_results_df_long)
    rm(rmaxes_results_df_long)
    gc()

    # Optionally make it a factor with specific order
    combined_rm_raw_df$axes_treatment <- factor(
    combined_rm_raw_df$axes_treatment,
    levels = c("Raw", "Axes Removed")
    )

    combined_rm_raw_df$abs_error <- abs(combined_rm_raw_df$error)
    combined_rm_raw_df$log_abs_error <- log(combined_rm_raw_df$abs_error + 0.001)


    lmm_model <- lmer(log_abs_error ~ rate * method * preservation_level * metric + 
                        axes_treatment +
                        axes_treatment:method + axes_treatment:preservation_level +
                        (1|replicate) + (1|tree_size), 
                        data = combined_rm_raw_df,
                        REML = FALSE)


    saveRDS(lmm_model, file.path("/mnt", "parscratch", "users", "bip24cns", "acedisparity", "discrete", "crown", "lmm", "lmm_discrete_combined_raw_rmaxes.rds"))
    cat("Finished combined raw and rmaxes LMM...\n")

}



