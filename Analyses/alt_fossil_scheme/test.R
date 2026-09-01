# Load required packages
library(ape)
library(phytools)
library(geiger)
library(cluster)

# ==============================================================================
# 1. SETUP & MOCK DATA GENERATION
# ==============================================================================
set.seed(42)

# Simulate a tree with deep time fossils
full_tree <- rtree(50)
t0 <- max(nodeHeights(full_tree)) * 0.7 # Intervention at 70% of tree height

# Simulate a discrete morphological matrix (e.g., 20 binary characters)
n_chars <- 20
Q_true <- matrix(c(-0.1, 0.1, 0.1, -0.1), 2, 2)
rownames(Q_true) <- colnames(Q_true) <- c("0", "1")
full_data <- sim.char(full_tree, Q_true, nsim = n_chars, model = "discrete")[,,]

# Induce taphonomic missingness (30% missing data)
full_data[sample(length(full_data), 0.3 * length(full_data))] <- NA

# Identify pre- and post-intervention tips
node_heights <- nodeHeights(full_tree)
tip_heights <- node_heights[full_tree$edge[, 2] <= Ntip(full_tree), 2]
names(tip_heights) <- full_tree$tip.label

pre_tips <- names(tip_heights[tip_heights <= t0])
post_tips <- names(tip_heights[tip_heights > t0])

# ==============================================================================
# 2. PRE-INTERVENTION PARAMETER ESTIMATION
# ==============================================================================
pre_tree <- drop.tip(full_tree, post_tips)
pre_data <- full_data[pre_tips, ]

# A. Estimate Taphonomic Missingness Rate
p_missing <- sum(is.na(pre_data)) / length(pre_data)

# B. Estimate Q-matrices for each character
# (In a real dataset, you would loop through chars and fit Mk models)
Q_estimates <- list()
for(i in 1:n_chars) {
  # Simplified rate estimation (using a uniform ER model for demonstration)
  # Replace with fitMk(pre_tree, pre_data[,i], model="ER") in real data
  est_rate <- 0.1 # Mock estimated rate
  Q_est <- matrix(c(-est_rate, est_rate, est_rate, -est_rate), 2, 2)
  rownames(Q_est) <- colnames(Q_est) <- c("0", "1")
  Q_estimates[[i]] <- Q_est
}

# C. Ancestral State Reconstruction at t0
# We need the states of lineages that cross t0. 
# Here, we extract the pre-intervention root/nodes to seed the simulation.
# (Simplified: using the root of the post-intervention clades)
post_tree <- drop.tip(full_tree, pre_tips)
anc_states_t0 <- sim.char(post_tree, Q_estimates[[1]], nsim = n_chars, model = "discrete")[1,,]
# ==============================================================================
# 3. FORWARD SIMULATION (THE NULL HYPOTHESIS) - CORRECTED
# ==============================================================================
n_simulations <- 100 # Number of counterfactual universes
simulated_matrices <- list()

for(s in 1:n_simulations) {
  # Create an empty matrix to hold the simulated traits for this run
  sim_post_mat <- matrix(NA, nrow = Ntip(post_tree), ncol = n_chars)
  rownames(sim_post_mat) <- post_tree$tip.label
  
  # Simulate each character independently using its specific Q matrix
  for(i in 1:n_chars) {
    sim_char_i <- sim.char(
      phy = post_tree, 
      par = Q_estimates[[i]], # Pass the specific single matrix for character i
      nsim = 1, 
      model = "discrete",     # Explicitly tell it to use Markov transition
      root = anc_states_t0[i] # Seed with the specific starting state for character i
    )
    sim_post_mat[, i] <- sim_char_i[,,1]
  }
  
  # Inject empirical taphonomic missingness
  sim_post_mat[sample(length(sim_post_mat), p_missing * length(sim_post_mat))] <- NA
  
  simulated_matrices[[s]] <- sim_post_mat
}

# ==============================================================================
# 4. ROBUST TRAITSPACE DISTANCE & NOVELTY TESTING (GOWER SPACE)
# ==============================================================================

# 1. Calculate Mean D_min for EMPIRICAL post-event data
emp_combined <- rbind(pre_data, emp_post_data)
# Suppress NA warnings from daisy
emp_gower <- suppressWarnings(as.matrix(cluster::daisy(as.data.frame(emp_combined), metric = "gower")))

# Extract only the rectangle of distances: Post-event (rows) vs Pre-event (cols)
idx_pre <- 1:nrow(pre_data)
idx_post <- (nrow(pre_data) + 1):nrow(emp_combined)
dist_emp_to_pre <- emp_gower[idx_post, idx_pre]

# Find the minimum distance to ANY pre-event tip for EACH post-event tip
dmin_emp_tips <- apply(dist_emp_to_pre, 1, min, na.rm = TRUE)
mean_dmin_emp <- mean(dmin_emp_tips)


# 2. Build the NULL DISTRIBUTION from all 100 simulations
null_mean_dmins <- numeric(n_simulations)

for(s in 1:n_simulations) {
  sim_combined <- rbind(pre_data, simulated_matrices[[s]])
  sim_gower <- suppressWarnings(as.matrix(cluster::daisy(as.data.frame(sim_combined), metric = "gower")))
  
  dist_sim_to_pre <- sim_gower[idx_post, idx_pre]
  dmin_sim_tips <- apply(dist_sim_to_pre, 1, min, na.rm = TRUE)
  
  null_mean_dmins[s] <- mean(dmin_sim_tips)
}

# 3. Calculate Results & True P-value
# P-value: Proportion of null simulations that explored LESS or EQUAL space than empirical
p_val_greater <- sum(mean_dmin_emp > null_mean_dmins) / n_simulations

cat("--- RESULTS ---\n")
cat("Empirical Mean D_min:", round(mean_dmin_emp, 4), "\n")
cat("Null Expectation Mean D_min:", round(mean(null_mean_dmins), 4), "\n")
cat("95% Null Envelope:", round(quantile(null_mean_dmins, 0.025), 4), "to", round(quantile(null_mean_dmins, 0.975), 4), "\n")
cat("P-value (Empirical expanded more than expected):", p_val_greater, "\n")