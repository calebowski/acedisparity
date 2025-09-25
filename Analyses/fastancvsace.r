#!/usr/bin/env Rscript
# Test script to isolate multi.ace issues
# filepath: test_multi_ace.R

library(treats)
library(ape)

cat("=== TESTING multi.ace function ===\n")

# Test 1: Simple synthetic data
cat("\n--- Test 1: Simple synthetic tree and data ---\n")
set.seed(123)

# Create simple tree
simple_tree <- rtree(20)
simple_tree$node.label <- paste0("n", 1:simple_tree$Nnode)

# Create simple trait matrix (20 tips x 10 traits)
simple_matrix <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10)
rownames(simple_matrix) <- simple_tree$tip.label
colnames(simple_matrix) <- paste0("trait_", 1:10)

cat("Simple tree tips:", length(simple_tree$tip.label), "\n")
cat("Simple matrix dims:", dim(simple_matrix), "\n")
cat("Simple matrix NAs:", sum(is.na(simple_matrix)), "\n")

tryCatch({
  simple_result <- multi.ace(data = simple_matrix, 
                            tree = simple_tree, 
                            models = "ML", 
                            output = "multi.ace")
  cat("✅ Simple test PASSED - no warnings\n")
}, warning = function(w) {
  cat("⚠️ Simple test WARNING:", w$message, "\n")
}, error = function(e) {
  cat("❌ Simple test ERROR:", e$message, "\n")
})

# Test 2: Replicate your exact data structure
cat("\n--- Test 2: Simulating your data structure ---\n")
set.seed(456)

# Simulate similar to your workflow
bd_params <- make.bd.params(speciation = 1, extinction = 0.7)
stop_rule <- list(max.living = 20)  # Smaller for testing

test_tree <- treats(stop.rule = stop_rule, bd.params = bd_params, null.error = 100)
test_tree <- drop.singles(test_tree)
cat("Test tree tips:", length(test_tree$tip.label), "\n")

# Create traits like yours
bm <- make.traits(process = BM.process, n = 50)  # Smaller n for testing
test_matrix_full <- map.traits(bm, test_tree)$data

cat("Full matrix dims (with nodes):", dim(test_matrix_full), "\n")

# *** FIX: Remove node states - keep only tip states ***
test_matrix <- test_matrix_full[test_tree$tip.label, , drop = FALSE]

cat("Fixed matrix dims (tips only):", dim(test_matrix), "\n")
cat("Test matrix NAs:", sum(is.na(test_matrix)), "\n")
cat("Test matrix range:", range(test_matrix), "\n")
cat("Matrix rows match tree tips:", all(rownames(test_matrix) %in% test_tree$tip.label), "\n")

cat("Test matrix dims:", dim(test_matrix), "\n")
cat("Test matrix NAs:", sum(is.na(test_matrix)), "\n")
cat("Test matrix range:", range(test_matrix), "\n")

# Test the anc.states function exactly as you use it
anc.states.test <- function(matrix, tree) {
  anc_states <- multi.ace(data = matrix, 
                          tree = tree, 
                          models = "ML", 
                          output = "multi.ace")
  return(anc_states)
}

tryCatch({
  test_result <- anc.states.test(test_matrix, test_tree)
  cat("✅ Your structure test PASSED - no warnings\n")
}, warning = function(w) {
  cat("⚠️ Your structure test WARNING:", w$message, "\n")
}, error = function(e) {
  cat("❌ Your structure test ERROR:", e$message, "\n")
})

# Test 3: Edge cases that might cause issues
cat("\n--- Test 3: Testing edge cases ---\n")

# Test with very small values
cat("Testing very small values...\n")
small_matrix <- test_matrix * 1e-10
tryCatch({
  small_result <- multi.ace(data = small_matrix, tree = test_tree, models = "ML", output = "multi.ace")
  cat("✅ Small values test PASSED\n")
}, warning = function(w) {
  cat("⚠️ Small values WARNING:", w$message, "\n")
})

# Test with very large values
cat("Testing very large values...\n")
large_matrix <- test_matrix * 1e10
tryCatch({
  large_result <- multi.ace(data = large_matrix, tree = test_tree, models = "ML", output = "multi.ace")
  cat("✅ Large values test PASSED\n")
}, warning = function(w) {
  cat("⚠️ Large values WARNING:", w$message, "\n")
})

# Test with extreme outliers
cat("Testing with outliers...\n")
outlier_matrix <- test_matrix
outlier_matrix[1, 1] <- 1000  # Add extreme outlier
tryCatch({
  outlier_result <- multi.ace(data = outlier_matrix, tree = test_tree, models = "ML", output = "multi.ace")
  cat("✅ Outliers test PASSED\n")
}, warning = function(w) {
  cat("⚠️ Outliers WARNING:", w$message, "\n")
})

# Test 4: Test the sampling function that might cause NAs
cat("\n--- Test 4: Testing sampling function ---\n")
trait_normal = list(fun = rnorm, param = list(mean = mean, sd = function(x)return(diff(range(x))/4)))

# Test sampling with your exact function
tryCatch({
  # First get a basic result
  basic_result <- multi.ace(data = test_matrix, tree = test_tree, models = "ML", output = "multi.ace")
  
  # Then test sampling
  sample_result <- multi.ace(basic_result, output = "combined.matrix", sample = 10, sample.fun = trait_normal)
  cat("✅ Sampling test PASSED\n")
}, warning = function(w) {
  cat("⚠️ Sampling WARNING:", w$message, "\n")
}, error = function(e) {
  cat("❌ Sampling ERROR:", e$message, "\n")
})

# Test 5: Check if it's specific to fossil preservation
cat("\n--- Test 5: Testing fossil preservation workflow ---\n")

# Load your fossil preservation functions
tryCatch({
  source("/users/bip24cns/acedisparity/discrete/scripts/fossil.pres.R")
  
  # Test fossil preservation
  test_matrices <- list(bm = test_matrix)
  fossilised_test <- lapply(test_matrices, fossil.pres, trees = test_tree, preservation = 0.5, type = "continuous")
  
  cat("Fossilised matrix dims:", dim(fossilised_test$bm$matrix), "\n")
  cat("Fossilised tree tips:", length(fossilised_test$bm$tree$tip.label), "\n")
  
  # Test ACE on fossilised data
  fossil_result <- multi.ace(data = fossilised_test$bm$matrix, 
                            tree = fossilised_test$bm$tree, 
                            models = "ML", 
                            output = "multi.ace")
  cat("✅ Fossil preservation test PASSED\n")
  
}, warning = function(w) {
  cat("⚠️ Fossil preservation WARNING:", w$message, "\n")
}, error = function(e) {
  cat("❌ Fossil preservation ERROR:", e$message, "\n")
})

cat("\n=== TEST SUMMARY ===\n")
cat("If all tests pass but your real data fails, the issue is in your specific data\n")
cat("If tests fail, the issue is with the multi.ace function or setup\n")
cat("Run this script to isolate the problem!\n")



# Test with fossil preservation and larger scale to trigger warnings
cat("\n--- Test 6: Large scale with fossil preservation (to trigger warnings) ---\n")
set.seed(789)

# Use your EXACT parameters from main script
bd_params <- make.bd.params(speciation = 1, extinction = 0.7)
stop_rule <- list(max.living = 50)  # Same as your main script

big_test_tree <- treats(stop.rule = stop_rule, bd.params = bd_params, null.error = 100)
big_test_tree <- drop.singles(big_test_tree)
cat("Big test tree tips:", length(big_test_tree$tip.label), "\n")

# Create traits with SAME parameters as your main script
bm_big <- make.traits(process = BM.process, n = 100)  # Same as your script
ou_strong_big <- make.traits(process = OU.process, n = 100, 
                             process.args = list(alpha = (log(2) / (max(node.depth.edgelength(big_test_tree)) / 10))))
ou_weak_big <- make.traits(process = OU.process, n = 100, 
                           process.args = list(alpha = (log(2) / max(node.depth.edgelength(big_test_tree)))))

# Map traits (this creates the full matrix with nodes)
bm_big_full <- map.traits(bm_big, big_test_tree)$data
ou_s_big_full <- map.traits(ou_strong_big, big_test_tree)$data
ou_w_big_full <- map.traits(ou_weak_big, big_test_tree)$data

cat("Full matrices dims (with nodes):\n")
cat("  BM:", dim(bm_big_full), "\n")
cat("  OU strong:", dim(ou_s_big_full), "\n")
cat("  OU weak:", dim(ou_w_big_full), "\n")

# Apply your fix: Remove node states
bm_big_matrices <- bm_big_full[big_test_tree$tip.label, , drop = FALSE]
ou_s_big_matrices <- ou_s_big_full[big_test_tree$tip.label, , drop = FALSE]
ou_w_big_matrices <- ou_w_big_full[big_test_tree$tip.label, , drop = FALSE]

cat("Fixed matrices dims (tips only):\n")
cat("  BM:", dim(bm_big_matrices), "\n")
cat("  OU strong:", dim(ou_s_big_matrices), "\n")
cat("  OU weak:", dim(ou_w_big_matrices), "\n")

# Test with fossil preservation (load your function first)
cat("Loading fossil preservation functions...\n")
tryCatch({
source("../Functions/fossil.pres.R")
  
  # Create matrices list like in your script
  big_matrices <- list(bm = bm_big_matrices, ou_w = ou_w_big_matrices, ou_s = ou_s_big_matrices)
  
  # Apply fossil preservation with SAME parameters as your script
  cat("Applying fossil preservation...\n")
  living <- lapply(big_matrices, remove.fossil, trees = big_test_tree, type = "continuous")
  fossilised_high <- lapply(big_matrices, fossil.pres, trees = big_test_tree, preservation = 0.5, type = "continuous")
  all_fossil <- lapply(big_matrices, fossil.pres, trees = big_test_tree, preservation = 1.0, type = "continuous")
  fossilised_med <- lapply(big_matrices, fossil.pres, trees = big_test_tree, preservation = 0.15, type = "continuous")
  fossilised_low <- lapply(big_matrices, fossil.pres, trees = big_test_tree, preservation = 0.05, type = "continuous")
  
  # Create fossil_matrices structure like your script
  big_fossil_matrices <- lapply(names(big_matrices), function(level) {
    list(
      all = all_fossil[[level]],
      fossil_high = fossilised_high[[level]],
      fossil_med = fossilised_med[[level]],
      fossil_low = fossilised_low[[level]],
      living = living[[level]]
    )
  })
  names(big_fossil_matrices) <- names(big_matrices)
  
  cat("Fossil matrices created:\n")
  for(rate in names(big_fossil_matrices)) {
    for(fossil_level in names(big_fossil_matrices[[rate]])) {
      dims <- dim(big_fossil_matrices[[rate]][[fossil_level]]$matrix)
      cat("  ", rate, "-", fossil_level, ": matrix", dims[1], "x", dims[2], 
          ", tree tips:", length(big_fossil_matrices[[rate]][[fossil_level]]$tree$tip.label), "\n")
    }
  }
  
  # Now test ACE on these problematic datasets
  cat("\n=== TESTING ACE ON FOSSIL MATRICES (should trigger warnings) ===\n")
  
  anc.states.big <- function(x) {
    anc_states <- multi.ace(data = x$matrix, 
                            tree = x$tree, 
                            models = "ML", 
                            output = "multi.ace")
    return(anc_states)
  }
  
  # Test each rate and fossil level
  for(rate_name in names(big_fossil_matrices)) {
    cat("Testing rate:", rate_name, "\n")
    for(fossil_name in names(big_fossil_matrices[[rate_name]])) {
      cat("  Testing fossil level:", fossil_name, "\n")
      tryCatch({
        result <- anc.states.big(big_fossil_matrices[[rate_name]][[fossil_name]])
        cat("    ✅ PASSED for", rate_name, "-", fossil_name, "\n")
      }, warning = function(w) {
        cat("    ⚠️ WARNING for", rate_name, "-", fossil_name, ":", w$message, "\n")
      }, error = function(e) {
        cat("    ❌ ERROR for", rate_name, "-", fossil_name, ":", e$message, "\n")
      })
    }
  }
  
}, error = function(e) {
  cat("❌ Could not load fossil preservation functions:", e$message, "\n")
  cat("Skipping fossil preservation test\n")
})

cat("\n=== TESTING EXTREME CONDITIONS ===\n")

# Test conditions that are likely to cause hessian warnings
cat("Testing with many traits and small sample...\n")
extreme_matrix <- test_matrix[1:5, , drop = FALSE]  # Very few tips, many traits
extreme_tree <- drop.tip(test_tree, test_tree$tip.label[6:length(test_tree$tip.label)])

tryCatch({
  extreme_result <- multi.ace(data = extreme_matrix, tree = extreme_tree, models = "ML", output = "multi.ace")
  cat("✅ Extreme conditions test PASSED\n")
}, warning = function(w) {
  cat("⚠️ Extreme conditions WARNING (expected):", w$message, "\n")
  cat("    This is likely the same warning you see in your main script!\n")
})

cat("\n=== SCALE TEST SUMMARY ===\n")
cat("This test uses your exact parameters and should reproduce the warnings\n")
cat("If you see 'sqrt(1/out$hessian): NaNs produced' here, that's the source!\n")cat("\n=== COMPARING multi.ace vs fastAnc ===\n")

# Load phytools for fastAnc
if(!require(phytools, quietly = TRUE)) {
  cat("Installing phytools...\n")
  install.packages("phytools")
  library(phytools)
}

# Function to test both methods on the same data
compare_ace_methods <- function(matrix_data, tree_data, name) {
  cat("\n--- Comparing methods for", name, "---\n")
  cat("Data: matrix", dim(matrix_data), ", tree tips:", length(tree_data$tip.label), "\n")
  
  # Test multi.ace (your current method)
  cat("Testing multi.ace...\n")
  multi_ace_warnings <- 0
  multi_ace_errors <- 0
  multi_ace_success <- FALSE
  
  tryCatch({
    multi_result <- multi.ace(data = matrix_data, 
                             tree = tree_data, 
                             models = "ML", 
                             output = "multi.ace")
    multi_ace_success <- TRUE
    cat("  ✅ multi.ace PASSED\n")
  }, warning = function(w) {
    multi_ace_warnings <<- multi_ace_warnings + 1
    cat("  ⚠️ multi.ace WARNING:", w$message, "\n")
    if(grepl("NaNs produced", w$message)) {
      cat("    *** This is your problematic warning! ***\n")
    }
  }, error = function(e) {
    multi_ace_errors <<- multi_ace_errors + 1
    cat("  ❌ multi.ace ERROR:", e$message, "\n")
  })
  
  # Test fastAnc on each trait individually (it only handles single traits)
  cat("Testing fastAnc (single traits)...\n")
  fastanc_warnings <- 0
  fastanc_errors <- 0
  fastanc_success <- 0
  
  # Test first 5 traits to avoid spam
  test_traits <- min(5, ncol(matrix_data))
  
  for(i in 1:test_traits) {
    single_trait <- matrix_data[, i]
    names(single_trait) <- rownames(matrix_data)
    
    tryCatch({
      fastanc_result <- fastAnc(tree_data, single_trait, vars = TRUE, CI = TRUE)
      fastanc_success <- fastanc_success + 1
      cat("    ✅ fastAnc trait", i, "PASSED\n")
    }, warning = function(w) {
      fastanc_warnings <<- fastanc_warnings + 1
      cat("    ⚠️ fastAnc trait", i, "WARNING:", w$message, "\n")
    }, error = function(e) {
      fastanc_errors <<- fastanc_errors + 1
      cat("    ❌ fastAnc trait", i, "ERROR:", e$message, "\n")
    })
  }
  
  # Test ape::ace on single trait for comparison
  cat("Testing ape::ace (single trait)...\n")
  single_trait <- matrix_data[, 1]
  names(single_trait) <- rownames(matrix_data)
  
  tryCatch({
    ace_result <- ape::ace(single_trait, tree_data, type = "continuous", method = "ML", CI = TRUE)
    cat("  ✅ ape::ace PASSED\n")
    
    # Check if CIs are reliable
    if(!is.null(ace_result$CI95) && any(is.na(ace_result$CI95))) {
      cat("    ⚠️ ape::ace has NAs in confidence intervals\n")
    }
    
  }, warning = function(w) {
    cat("  ⚠️ ape::ace WARNING:", w$message, "\n")
    if(grepl("NaNs produced", w$message)) {
      cat("    *** Same problem as multi.ace! ***\n")
    }
  }, error = function(e) {
    cat("  ❌ ape::ace ERROR:", e$message, "\n")
  })
  
  # Summary for this dataset
  cat("\n--- Summary for", name, "---\n")
  cat("multi.ace: warnings =", multi_ace_warnings, ", errors =", multi_ace_errors, "\n")
  cat("fastAnc: warnings =", fastanc_warnings, ", errors =", fastanc_errors, ", successes =", fastanc_success, "/", test_traits, "\n")
  
  return(list(
    multi_ace_warnings = multi_ace_warnings,
    fastanc_warnings = fastanc_warnings,
    name = name
  ))
}

# Test on your problematic datasets
results <- list()

# Test 1: Simple case (should work for all)
if(exists("simple_matrix") && exists("simple_tree")) {
  results$simple <- compare_ace_methods(simple_matrix, simple_tree, "Simple case")
}

# Test 2: BM high variance (should show problems)
if(exists("bm_big_matrices") && exists("big_test_tree")) {
  # Use subset to make it manageable
  bm_subset <- bm_big_matrices[, 1:10, drop = FALSE]
  results$bm <- compare_ace_methods(bm_subset, big_test_tree, "BM high variance")
}

# Test 3: OU_strong (should be better)
if(exists("ou_s_big_matrices") && exists("big_test_tree")) {
  ou_s_subset <- ou_s_big_matrices[, 1:10, drop = FALSE]
  results$ou_s <- compare_ace_methods(ou_s_subset, big_test_tree, "OU strong")
}

# Test 4: OU_weak (should have problems like BM)
if(exists("ou_w_big_matrices") && exists("big_test_tree")) {
  ou_w_subset <- ou_w_big_matrices[, 1:10, drop = FALSE]
  results$ou_w <- compare_ace_methods(ou_w_subset, big_test_tree, "OU weak")
}

# Test 5: Extreme case (very few tips, many traits)
if(exists("test_matrix") && exists("test_tree")) {
  extreme_tips <- min(8, length(test_tree$tip.label))
  extreme_matrix <- test_matrix[1:extreme_tips, 1:20, drop = FALSE]  # 8 tips, 20 traits
  extreme_tree_subset <- drop.tip(test_tree, test_tree$tip.label[(extreme_tips+1):length(test_tree$tip.label)])
  results$extreme <- compare_ace_methods(extreme_matrix, extreme_tree_subset, "Extreme case")
}

cat("\n=== FINAL COMPARISON SUMMARY ===\n")
cat("Method comparison across all test cases:\n")

total_multi_warnings <- sum(sapply(results, function(x) x$multi_ace_warnings))
total_fastanc_warnings <- sum(sapply(results, function(x) x$fastanc_warnings))

cat("Total multi.ace warnings:", total_multi_warnings, "\n")
cat("Total fastAnc warnings:", total_fastanc_warnings, "\n")

if(total_multi_warnings > 0 && total_fastanc_warnings == 0) {
  cat("🎯 CONCLUSION: multi.ace has specific issues that fastAnc doesn't\n")
  cat("   Consider using fastAnc for individual traits if multi.ace fails\n")
} else if(total_multi_warnings > 0 && total_fastanc_warnings > 0) {
  cat("🎯 CONCLUSION: Both methods struggle with high-variance data\n")
  cat("   This is a fundamental numerical issue, not specific to multi.ace\n")
} else {
  cat("🎯 CONCLUSION: No major differences detected in this test\n")
}

cat("\n=== PRACTICAL RECOMMENDATIONS ===\n")
if(total_multi_warnings > 0) {
  cat("Since you're getting warnings:\n")
  cat("1. ✅ Your point estimates are still valid\n")
  cat("2. ❌ Don't trust confidence intervals when NaNs are produced\n")
  cat("3. 🔄 Consider fastAnc for single-trait analyses\n")
  cat("4. 🔄 Consider reducing the number of traits processed simultaneously\n")
  cat("5. 🔄 Consider using more constrained evolutionary models (stronger OU)\n")
}







cat("\n=== SIMPLE COMPARISON: multi.ace vs fastAnc on BM data ===\n")

if(exists("bm_big_matrices") && exists("big_test_tree")) {
  
  # Test 1: multi.ace on all traits at once
  cat("Testing multi.ace on BM data (all 100 traits)...\n")
  tryCatch({
    multi_ace_result <- multi.ace(data = bm_big_matrices, 
                                 tree = big_test_tree, 
                                 models = "ML", 
                                 output = "multi.ace")
    cat("✅ multi.ace COMPLETED\n")
    cat("Result structure:", names(multi_ace_result), "\n")
    if("ace" %in% names(multi_ace_result)) {
      cat("Ancestral states dimensions:", dim(multi_ace_result$ace), "\n")
    }
  }, warning = function(w) {
    cat("⚠️ multi.ace WARNING:", w$message, "\n")
  }, error = function(e) {
    cat("❌ multi.ace ERROR:", e$message, "\n")
  })
  
  # Test 2: fastAnc on each trait individually, then bind
  cat("\nTesting fastAnc on BM data (trait by trait)...\n")
  
  # Check if phytools is available
  if(!require(phytools, quietly = TRUE)) {
    cat("Installing phytools...\n")
    install.packages("phytools")
    library(phytools)
  }
  
  fastanc_results <- list()
  fastanc_warnings <- 0
  fastanc_errors <- 0
  
  for(i in 1:ncol(bm_big_matrices)) {
    single_trait <- bm_big_matrices[, i]
    names(single_trait) <- rownames(bm_big_matrices)
    
    tryCatch({
      fastanc_result <- fastAnc(big_test_tree, single_trait, vars = FALSE, CI = FALSE)
      fastanc_results[[i]] <- fastanc_result
      
      if(i %% 20 == 0) {  # Progress indicator
        cat("  Completed", i, "/", ncol(bm_big_matrices), "traits\n")
      }
      
    }, warning = function(w) {
      fastanc_warnings <<- fastanc_warnings + 1
      fastanc_results[[i]] <<- rep(NA, length(big_test_tree$tip.label) - 1)
      if(i <= 5) cat("    ⚠️ WARNING trait", i, ":", w$message, "\n")
    }, error = function(e) {
      fastanc_errors <<- fastanc_errors + 1
      fastanc_results[[i]] <<- rep(NA, length(big_test_tree$tip.label) - 1)
      if(i <= 5) cat("    ❌ ERROR trait", i, ":", e$message, "\n")
    })
  }
  
  # Bind fastAnc results into a matrix
  if(length(fastanc_results) > 0) {
    fastanc_matrix <- do.call(cbind, fastanc_results)
    rownames(fastanc_matrix) <- paste0("Node_", (length(big_test_tree$tip.label) + 1):(length(big_test_tree$tip.label) + big_test_tree$Nnode))
    colnames(fastanc_matrix) <- colnames(bm_big_matrices)
    
    cat("✅ fastAnc COMPLETED\n")
    cat("Combined matrix dimensions:", dim(fastanc_matrix), "\n")
    cat("fastAnc warnings:", fastanc_warnings, "/ errors:", fastanc_errors, "\n")
    
    # Quick comparison
    cat("\n--- COMPARISON ---\n")
    if(exists("multi_ace_result") && "ace" %in% names(multi_ace_result)) {
      cat("multi.ace ancestral states range:", range(multi_ace_result$ace, na.rm = TRUE), "\n")
    }
    cat("fastAnc ancestral states range:", range(fastanc_matrix, na.rm = TRUE), "\n")
    cat("fastAnc NAs:", sum(is.na(fastanc_matrix)), "\n")
    
    # Check if results are similar (first few traits)
    if(exists("multi_ace_result") && "ace" %in% names(multi_ace_result) && 
       nrow(multi_ace_result$ace) == nrow(fastanc_matrix)) {
      
      cat("Correlation between methods (first 5 traits):\n")
      for(i in 1:min(5, ncol(fastanc_matrix))) {
        if(!any(is.na(fastanc_matrix[, i]))) {
          corr <- cor(multi_ace_result$ace[, i], fastanc_matrix[, i])
          cat("  Trait", i, "correlation:", round(corr, 4), "\n")
        }
      }
    }
    
  } else {
    cat("❌ fastAnc failed completely\n")
  }
  
} else {
  cat("❌ BM matrices or tree not found\n")
}

cat("\n=== SIMPLE SUMMARY ===\n")
cat("This shows whether fastAnc avoids the warnings that multi.ace produces\n")
cat("and whether the ancestral state estimates are similar between methods\n")



# Get multi.ace output
multi_ace_result <- multi.ace(data = bm_big_matrices, 
                             tree = big_test_tree, 
                             models = "ML", 
                             output = "multi.ace")

# Get fastAnc output for all traits
fastanc_results <- list()
for(i in 1:ncol(bm_big_matrices)) {
  single_trait <- bm_big_matrices[, i]
  names(single_trait) <- rownames(bm_big_matrices)
  fastanc_results[[i]] <- fastAnc(big_test_tree, single_trait, vars = FALSE, CI = TRUE)
}

# Bind into matrix
fastanc_matrix <- do.call(cbind, fastanc_results)
rownames(fastanc_matrix) <- paste0("Node_", (length(big_test_tree$tip.label) + 1):(length(big_test_tree$tip.label) + big_test_tree$Nnode))
colnames(fastanc_matrix) <- colnames(bm_big_matrices)

# Now you have:
# multi_ace_result - the multi.ace output
# fastanc_matrix - the fastAnc output matrix
#







# Compare CIs for all traits
compare_all_CIs <- function(multi_result, fastanc_list) {
  
  # Extract multi.ace CIs (assuming structure like multi_result$continuous$estimates[[1]]$A1$CI95)
  multi_CIs <- list()
  for(i in 1:length(fastanc_list)) {
    trait_name <- paste0("A", i)
    if(!is.null(multi_result$continuous$estimates[[1]][[trait_name]]$CI95)) {
      multi_CIs[[i]] <- multi_result$continuous$estimates[[1]][[trait_name]]$CI95
    }
  }
  
  # Extract fastAnc CIs
  fastanc_CIs <- lapply(fastanc_list, function(x) x$CI95)
  
  # Compare each trait
  correlations <- numeric(length(fastanc_list))
  mean_differences <- numeric(length(fastanc_list))
  max_differences <- numeric(length(fastanc_list))
  
  for(i in 1:length(fastanc_list)) {
    if(!is.null(multi_CIs[[i]]) && !is.null(fastanc_CIs[[i]])) {
      
      # Flatten both CI matrices for comparison
      multi_flat <- as.vector(multi_CIs[[i]])
      fast_flat <- as.vector(fastanc_CIs[[i]])
      
      # Check for NAs
      if(!any(is.na(multi_flat)) && !any(is.na(fast_flat)) && length(multi_flat) == length(fast_flat)) {
        correlations[i] <- cor(multi_flat, fast_flat)
        differences <- multi_flat - fast_flat
        mean_differences[i] <- mean(abs(differences))
        max_differences[i] <- max(abs(differences))
      } else {
        correlations[i] <- NA
        mean_differences[i] <- NA
        max_differences[i] <- NA
      }
    } else {
      correlations[i] <- NA
      mean_differences[i] <- NA
      max_differences[i] <- NA
    }
  }
  
  return(data.frame(
    trait = 1:length(fastanc_list),
    correlation = correlations,
    mean_abs_diff = mean_differences,
    max_abs_diff = max_differences
  ))
}

# Run the comparison
CI_comparison <- compare_all_CIs(multi_ace_result, fastanc_results)

# Summary statistics
summary(CI_comparison$correlation, na.rm = TRUE)
summary(CI_comparison$mean_abs_diff, na.rm = TRUE)
summary(CI_comparison$max_abs_diff, na.rm = TRUE)

# Count how many traits have reliable CI comparisons
reliable_CIs <- sum(!is.na(CI_comparison$correlation))
cat("Traits with reliable CI comparisons:", reliable_CIs, "out of", nrow(CI_comparison), "\n")

# Show the comparison dataframe
CI_comparison


