library(paleotree)

p <- 0.75 ## speciation
q <- 0.95 * p   ## extinction
rho <- 1 / 8 # intended sampling probability
r <- -q / (1 - (1 / rho)) ## this is the sampling rate

nchars <- 136
nstates <- c(
  rep(2, nchars / 2),
  rep(3, nchars / 4),
  rep(4, nchars / 4)
)

nTotalTaxa <- 100 / rho

# Simulate fossil phylogeny and taxon histories
simtrees <- simFossilRecord(
  p = p,
  q = q,
  r = r,
  nTotalTaxa = nTotalTaxa
)

# 2. Extract Full True Phylogeny
taxa_full <- fossilRecord2fossilTaxa(simtrees)
tree_full <- taxa2phylo(taxa_full)

# 3. Extract Fossil Sampled Phylogeny
# Identify which taxa were sampled at least once
sampled_flags <- sapply(simtrees, function(x) length(x$sampling.times) > 0)
unsampled_taxa <- tree_full$tip.label[!sampled_flags]

# Prune unsampled lineages to leave only fossil-sampled taxa
tree_sampled <- drop.tip(tree_full, tip = unsampled_taxa)

# Extract taxon information
taxon_info <- data.frame(
  sp = 1:nTotalTaxa,
  anc = NA_integer_,
  birth = NA_real_,
  death = NA_real_
)

ancestor_id <- match("ancestor.id", names(simtrees$t1$taxa.data))
origin_time <- match("orig.time", names(simtrees$t1$taxa.data))
extinction_time <- match("ext.time", names(simtrees$t1$taxa.data))

for (taxon in seq_len(nTotalTaxa)) {
  taxon_info$anc[taxon] <-
    simtrees[[taxon]]$taxa.data[ancestor_id]

  taxon_info$birth[taxon] <-
    simtrees[[taxon]]$taxa.data[origin_time]

  taxon_info$death[taxon] <-
    simtrees[[taxon]]$taxa.data[extinction_time]
}

# Simulate discrete traits along the phylogeny
sim_chmatrix <- matrix(
  0,
  nrow = nTotalTaxa,
  ncol = nchars
)

for (taxon in 2:nTotalTaxa) {
  ancestor <- taxon_info$anc[taxon]

  sim_chmatrix[taxon, ] <- sim_chmatrix[ancestor, ]

  n_changes <- 0
  while (n_changes == 0) {
    n_changes <- max(1, rbinom(1, nchars, prob = 0.06 * rho))
  }

  changed_characters <- sample(
    seq_len(nchars),
    size = n_changes,
    replace = FALSE
  )

  for (character in changed_characters) {
    ancestor_state <- sim_chmatrix[ancestor, character]

    sim_chmatrix[taxon, character] <-
      sample(
        setdiff(seq_len(nstates[character]) - 1, ancestor_state),
        size = 1
      )
  }
}

rownames(sim_chmatrix) <- taxon_info$sp

output <- list(
  true_tree= tree_full,
  sampled_tree = tree_sampled,
  taxon_info = taxon_info,
  trait_matrix = sim_chmatrix
)