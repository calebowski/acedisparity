library(paleotree)
library(ape)

# 1. Parameters
p <- 0.75             # Speciation rate
q <- 0.95 * p         # Extinction rate
nTotalTaxa <- 200     # Fixed total taxa for the single true tree
rho_levels <- c(rho_50 = 0.50, rho_15 = 0.15, rho_05 = 0.05)

nchars <- 136
nstates <- c(
  rep(2, nchars / 2),
  rep(3, nchars / 4),
  rep(4, nchars / 4)
)

# 2. Simulate 1 Full True Phylogeny (r = 0: no sampling during diversification)
simtrees <- simFossilRecord(
  p = p,
  q = q,
  r = 0,
  nTotalTaxa = nTotalTaxa
)

# Convert to fossilTaxa object and ape phylo tree
taxa_full <- fossilRecord2fossilTaxa(simtrees)
tree_full <- taxa2phylo(taxa_full)

# Extract taxon information
taxon_info <- data.frame(
  sp = seq_len(nTotalTaxa),
  anc = taxa_full[, "ancestor.id"],
  birth = taxa_full[, "orig.time"],
  death = taxa_full[, "ext.time"]
)

# 3. Simulate Traits Once Across the Full True Tree
sim_chmatrix <- matrix(0, nrow = nTotalTaxa, ncol = nchars)

for (taxon in 2:nTotalTaxa) {
  ancestor <- taxon_info$anc[taxon]
  sim_chmatrix[taxon, ] <- sim_chmatrix[ancestor, ]

  n_changes <- 0
  while (n_changes == 0) {
    n_changes <- max(1, rbinom(1, nchars, prob = 0.06))
  }

  changed_characters <- sample(
    seq_len(nchars),
    size = n_changes,
    replace = FALSE
  )

  for (character in changed_characters) {
    ancestor_state <- sim_chmatrix[ancestor, character]
    sim_chmatrix[taxon, character] <- sample(
      setdiff(seq_len(nstates[character]) - 1, ancestor_state),
      size = 1
    )
  }
}


rownames(sim_chmatrix) <- taxon_info$sp

# 4. Post-Hoc Sampling at Different Rho Levels
sampled_results <- lapply(rho_levels, function(rho) {
  # Convert target sampling probability (rho) to continuous sampling rate (r)
  r_val <- -q / (1 - (1 / rho))
  
  # Simulate fossil recovery on true lifespans using sampleRanges
  sampled_ranges <- sampleRanges(taxa_full, r = r_val, ranges.only = TRUE)
  
  # Identify which taxa were sampled at least once (non-NA appearance dates)
  sampled_flags <- !is.na(sampled_ranges[, 1])
  
  # Prune unsampled lineages from the true tree
  unsampled_tips <- tree_full$tip.label[!sampled_flags]
  tree_sampled <- drop.tip(tree_full, tip = unsampled_tips)
  
  list(
    rho = rho,
    r_rate = r_val,
    sampled_flags = sampled_flags,
    sampled_tree = tree_sampled,
    trait_matrix_sampled = sim_chmatrix[sampled_flags, ]
  )
})

# Final output object containing 1 true tree and sub-sampled results
output <- list(
  true_tree = tree_full,
  taxon_info = taxon_info,
  trait_matrix_full = sim_chmatrix,
  sampled_results = sampled_results
)


