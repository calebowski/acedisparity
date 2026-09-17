library(dispRity)


tree <- read.tree("../../Data/trees/tree_50t_001.tre")


extract.crown.tree <- function(tree){
  ages <- tree.age(tree) # get tip ages
  extant <- ages$element[ages$ages == 0]
  living_tree <- keep.tip(tree, extant) ## root at mrca - alike to crown group analyses.

  crown_tree  <- extract.clade(tree, living_tree$node.label[1]) ## this is the crown tree
  return(crown_tree)
}




ages <- tree.age(tree) # get tip ages
extant <- ages$element[ages$ages == 0]
living_tree <- keep.tip(tree, extant) ## root at mrca - alike to crown group analyses.

crown_tree  <- extract.clade(tree, living_tree$node.label[1]) ## this is the crown tree

fossil_matrices <- readRDS("..//../Data/revisions/wagner/discrete/matrices/11429393_fossil_matrices_001.rds")
ord_no_ace <- readRDS("..//../Data/revisions/wagner/discrete/ord/11429393_ord_no_ace_001.rds")

ord_true <- readRDS("..//../Data/revisions/wagner/discrete/ord/11429393_ord_true_001.rds")

ord_point <- readRDS("..//../Data/revisions/wagner/discrete/ord/11429393_ord_point_001.rds")

post_ord_point_fossil_med <- readRDS("..//../Data/revisions/wagner/discrete/ord/11429393_fossil_med_point_post_ord_ace.rds")

# par(mfrow = c(2,2))
# plot(ord_no_ace$slow$fossil_high)
# plot(ord_point$slow$fossil_high)
# plot(ord_true$slow)



# get.disparity(dispRity(ord_no_ace$slow$fossil_high, metric = c(sum, variances)))
# get.disparity(dispRity(ord_point$slow$fossil_high, metric = c(sum, variances)))
# get.disparity(dispRity(ord_true$slow, metric = c(sum, variances)))


# lapply(ords, function(x) get.disparity(dispRity(x, metric = c(sum, variances))))

# Use the same two axes and limits for every panel
spaces <- list(
  true = ord_true$slow,
  no_ace = ord_no_ace$slow$fossil_med,
  ace = ord_point$slow$fossil_med
)

fossil_tree <- fossil_matrices$slow$fossil_med$tree

# all_points <- do.call(rbind, spaces)
# xlim <- range(all_points[, 1], na.rm = TRUE)
# ylim <- range(all_points[, 2], na.rm = TRUE)


plot_phylogeny <- function(x, tree, col = "grey70", lwd = 1) {
  coordinate_names <- rownames(x)

  tree_names <- c(tree$tip.label, tree$node.label)

  edge_names <- cbind(
    parent = tree_names[tree$edge[, 1]],
    child = tree_names[tree$edge[, 2]]
  )

  parent_index <- match(edge_names[, "parent"], coordinate_names)
  child_index <- match(edge_names[, "child"], coordinate_names)

  valid <- !is.na(parent_index) & !is.na(child_index)

  message(
    "Matched edges: ", sum(valid), "/", nrow(edge_names),
    "; matched tips: ",
    sum(tree$tip.label %in% coordinate_names), "/", length(tree$tip.label),
    "; matched nodes: ",
    sum(tree$node.label %in% coordinate_names), "/", length(tree$node.label)
  )

  if (!any(valid)) {
    warning("No tree edges have both endpoints in the ordination.")
    return(invisible(NULL))
  }

  for (i in which(valid)) {
    segments(
      x[parent_index[i], 1], x[parent_index[i], 2],
      x[child_index[i], 1], x[child_index[i], 2],
      col = col,
      lwd = lwd
    )
  }

  invisible(NULL)
}


plot_space <- function(x, title, tree = NULL) {
  row_names <- rownames(x)


  metrics <- list(
    sum_variances = c(sum, variances),
    sum_quantiles = c(sum, quantiles),
    mean_pairwise = c(mean, pairwise.dist.na.rm)
  )

  disparity_values <- sapply(
    metrics,
    function(metric) {
      get.disparity(dispRity(x, metric = metric))[[1]]
    }
  )

  sampled_nodes <- grepl("^f_n", row_names)
  ordinary_tips <- grepl("^t", row_names)
  ordinary_nodes <- grepl("^n", row_names)

    plot(
    x[, 1], x[, 2],
    type = "n",
    main = title,
    xlab = "PCoA axis 1",
    ylab = "PCoA axis 2",
    cex = 1.5
  )

  if (!is.null(tree)) {
    plot_phylogeny(x, tree)
  }

  # Ordinary tips include living and extinct tips labelled t...
  points(
    x[ordinary_tips, 1],
    x[ordinary_tips, 2],
    pch = 19,
    col = "steelblue"
  )

  # Sampled ancestral nodes are represented as fossil tips f_n...
  points(
    x[sampled_nodes, 1],
    x[sampled_nodes, 2],
    pch = 19,
    col = "darkorange"
  )

  # Reconstructed or true ancestral states n...
  points(
    x[ordinary_nodes, 1],
    x[ordinary_nodes, 2],
    pch = 19,
    col = "firebrick"
  )

    sampled <- which(sampled_nodes)
  ancestral <- which(ordinary_nodes)

  if (length(sampled) > 0 && length(ancestral) > 0) {
    distances <- as.matrix(dist(rbind(x[sampled, , drop = FALSE],
                                       x[ancestral, , drop = FALSE])))

    n_sampled <- length(sampled)
    overlap_distances <- distances[
      seq_len(n_sampled),
      n_sampled + seq_along(ancestral),
      drop = FALSE
    ]

    overlap <- which(overlap_distances < 0.01, arr.ind = TRUE)

    if (nrow(overlap) > 0) {
      overlap_sampled <- sampled[overlap[, 1]]
      overlap_ancestral <- ancestral[overlap[, 2]]

      # Highlight coordinates occupied by both types of node
      points(
        x[overlap_sampled, 1],
        x[overlap_sampled, 2],
        pch = 21,
        bg = "yellow",
        col = "black",
        cex = 2
      )

      text(
        x[overlap_sampled, 1],
        x[overlap_sampled, 2],
        labels = paste(
          rownames(x)[overlap_sampled],
          rownames(x)[overlap_ancestral],
          sep = " / "
        ),
        pos = 3,
        cex = 0.7
      )

      message("Overlapping nodes:")
      print(data.frame(
        sampled_node = rownames(x)[overlap_sampled],
        ancestral_node = rownames(x)[overlap_ancestral],
        x = x[overlap_sampled, 1],
        y = x[overlap_sampled, 2]
      ))
    }
  }

  legend(
    "topleft",
    legend = c(
      "Tips (living & fossil)",
      "Sampled fossil nodes",
      "Ancestral nodes"
    ),
    pch = 19,
    col = c("steelblue", "darkorange", "firebrick"),
    bty = "n",
    cex = 1.8
  )
  legend(
  "topright",
  legend = c(
    sprintf("Sum variances: %.3f", disparity_values["sum_variances"]),
    sprintf("Sum quantiles: %.3f", disparity_values["sum_quantiles"]),
    sprintf("Mean pairwise: %.3f", disparity_values["mean_pairwise"])
  ),
  bty = "n",
  cex = 1)
}
  
par(mfrow = c(2, 2))

plot_space(spaces$true, "True", crown_tree)
plot_space(spaces$no_ace, "No ASE: fossil_med")
plot_space(spaces$ace, "ASE: fossil_med")



rm_axes <- lapply(spaces, function(x) x[,1:48])


lapply(rm_axes, function(x) dispRity(x, metric = c(sum, quantiles))$disparity)


plot_space(spaces$post_ord_ace, "post ord ACE: fossil_med")


get.disparity(dispRity(spaces$true[grepl("^t", rownames(spaces$true)),], metric = c(sum, ranges))) ## tips and nodes space have the same disparity value
get.disparity(dispRity(spaces$true[grepl("^n", rownames(spaces$true)),], metric = c(sum, ranges)))

get.disparity(dispRity(spaces$no_ace[grepl("^f", rownames(spaces$no_ace)),], metric = c(sum, ranges)))
get.disparity(dispRity(spaces$no_ace[grepl("^t", rownames(spaces$no_ace)),], metric = c(sum, ranges)))


get.disparity(dispRity(spaces$ace[grepl("^t", rownames(spaces$ace)),], metric = c(sum, ranges)))
get.disparity(dispRity(spaces$ace[grepl("^f", rownames(spaces$ace)),], metric = c(sum, ranges))) ## sampling nodes has lower disparity than ancestral state estimation, which is likely the reason why wagner had it wrong?#

## do i need to do a contrast of disparity in wagner's way, which is using the ancestral nodes by range extension?
get.disparity(dispRity(spaces$ace[grepl("^n", rownames(spaces$ace)),], metric = c(sum, ranges)))



## do this over all 100, check to see if there is any difference in disparity space using an anova??









### TIP TO PARENT TEST
ord_true_list <- list()
trees <- list()
for (i in 1:100){
  # ord_true_list[[i]] <- readRDS(paste0(sprintf("../../Data/revisions/wagner2/discrete/ord/11532243_ord_true_%03d.rds", i)))
  trees[[i]] <- extract.crown.tree(read.tree(sprintf("../../Data/trees/tree_50t_%03d.tre", i)))
}

get.parent.child <- function(tree) {
  tree_names <- c(tree$tip.label, tree$node.label)

  parent_child <- data.frame(
    parent = tree_names[tree$edge[, 1]],
    child = tree_names[tree$edge[, 2]],
    branch_length = tree$edge.length,
    tip_vs_node = ifelse(grepl("^n", (tree_names[tree$edge[, 2]])), "node", "tip"),
    stringsAsFactors = FALSE
  )

  return(parent_child)
}

parent_children <- lapply(trees, get.parent.child)
parent_children_df <- do.call(rbind, parent_children)

boxplot(
    (branch_length) ~ tip_vs_node,
    data = parent_children_df,
    xlab = "",
    ylab = "Branch length"
)





parent.child.dist <- function(parent_child, ordination){
  distance_mat <- as.matrix(dist(ordination, method = "euclidean"))
  for (i in seq_len(nrow(parent_child))) {
    parent_child[i, "morphological_distance"] <- distance_mat[parent_child[i,"parent"], parent_child[i,"child"]]
  }
  return(parent_child)
}

parent_child_dist <- list()
for (i in 1:100){
  pc <- parent_children[[i]]
  ord <- ord_true_list[[i]]
  parent_child_dist[[i]] <- lapply(ord, function(rate){
    parent.child.dist(pc, rate)
  })
}


slow_comparison <- do.call(rbind, lapply(parent_child_dist, function(x) x$slow))
plot(slow_comparison$morphological_distance ~sqrt(slow_comparison$branch_length))

med_comparison <- do.call(rbind, lapply(parent_child_dist, function(x) x$med))
plot(med_comparison$morphological_distance ~sqrt(med_comparison$branch_length))


fast_comparison <- do.call(rbind, lapply(parent_child_dist, function(x) x$fast))
plot(fast_comparison$morphological_distance ~sqrt(fast_comparison$branch_length))


slow_comparison_internal <- subset(slow_comparison, grepl("^n", parent) & grepl("^n", child))
slow_comparison_terminal <- subset(slow_comparison,  grepl("^t", child))



### should we not be looking at distance from centroid, tips vs nodes???

distance.from.centroid <- function(ordination) {
  coordinates <- ordination[, , drop = FALSE]

  centroid <- colMeans(coordinates, na.rm = TRUE)

  data.frame(
    name = rownames(coordinates),
    type = ifelse(
      grepl("^t", rownames(coordinates)), "tip",
      ifelse(grepl("^n", rownames(coordinates)), "node", "other")
    ),
    distance_from_centroid = sqrt(
      rowSums((sweep(coordinates, 2, centroid, "-"))^2)
    )
  )
}

# centroid_distances <- do.call(rbind, do.call(rbind, lapply(ord_true_list, lapply, distance.from.centroid)))


# boxplot(
#   distance_from_centroid ~ type,
#   data = subset(centroid_distances, type %in% c("tip", "node")),
#   xlab = "",
#   ylab = "Distance from centroid",
#   cex.axis = 1.4,
#   cex.lab = 1.5
# )


centroid_distances <- do.call(
  rbind,
  lapply(seq_along(ord_true_list), function(i) {
    do.call(
      rbind,
      lapply(names(ord_true_list[[i]]), function(rate) {
        data <- distance.from.centroid(ord_true_list[[i]][[rate]])
        data$transition_rate <- rate
        data$replicate <- i
        data
      })
    )
  })
)

centroid_distances <- subset(
  centroid_distances,
  type %in% c("tip", "node")
)

boxplot(
  distance_from_centroid ~ interaction(transition_rate, type),
  data = centroid_distances,
  xlab = "Transition rate and tip or node",
  ylab = "Distance from centroid",
  cex.axis = 1.2,
  cex.lab = 1.5
)



centroid_summary <- do.call(
  rbind,
  lapply(seq_along(ord_true_list), function(i) {
    replicate_data <- do.call(
      rbind,
      lapply(ord_true_list[[i]], distance.from.centroid)
    )

    aggregate(
      distance_from_centroid ~ type,
      data = subset(replicate_data, type %in% c("tip", "node")),
      FUN = mean
    ) |>
      transform(replicate = i)
  })
)

tips <- subset(centroid_summary, type == "tip")
nodes <- subset(centroid_summary, type == "node")

wilcox.test(
  nodes$distance_from_centroid,
  tips$distance_from_centroid,
  paired = TRUE
)


nearest_tip_distance <- function(ordination) {
  tip_index <- grepl("^t", rownames(ordination))
  node_index <- grepl("^n", rownames(ordination))

  distances <- as.matrix(dist(ordination))

  data.frame(
    node = rownames(ordination)[node_index],
    nearest_tip_distance = apply(
      distances[node_index, tip_index, drop = FALSE],
      1,
      min
    )
  )
}
nearest_tip <- do.call(
  rbind,
  lapply(seq_along(ord_true_list), function(i) {
    data <- nearest_tip_distance(ord_true_list[[i]]$slow)
    data$replicate <- i
    data
  })
)

summary(nearest_tip$nearest_tip_distance)


## show nodes that descended from another node, vs nodes that descended from tips. this can prove the point that derived nodes, clsoer to end of tree, are important to sample because they fill in the gaps.




## TIP TO PARENT TEST WITH PUNCTUATED

ord_true_list <- list()
trees <- list()
for (i in 1:100){
  ord_true_list[[i]] <- readRDS(paste0(sprintf("../../Data/revisions/wagner2/discrete/ord/11532243_ord_true_%03d.rds", i)))
  trees[[i]] <- extract.crown.tree(read.tree(sprintf("../../Data/trees/tree_50t_%03d.tre", i)))
}

slow <- lapply(ord_true_list, function(x) x$slow)







distance.from.centroid <- function(ordination) {
  coordinates <- ordination[, , drop = FALSE]
  names <- rownames(coordinates)

  centroid <- colMeans(coordinates, na.rm = TRUE)

  data.frame(
    name = names,
    type = ifelse(
      grepl("^t", names), "tip",
      ifelse(
        grepl("^f_n", names), "node",
        ifelse(grepl("^n", names), "node", "other")
      )
    ),
    distance_from_centroid = sqrt(
      rowSums(
        sweep(coordinates, 2, centroid, "-")^2,
        na.rm = TRUE
      )
    )
  )
}

centroid_distances <- do.call(rbind, do.call(rbind, lapply(ord_true_list, lapply, distance.from.centroid)))


boxplot(
  distance_from_centroid ~ type,
  data = subset(centroid_distances, type %in% c("tip", "node")),
  xlab = "",
  ylab = "Distance from shared centroid"
)


centroid_summary <- do.call(
  rbind,
  lapply(seq_along(ord_true_list), function(i) {
    replicate_data <- do.call(
      rbind,
      lapply(ord_true_list[[i]], distance.from.centroid)
    )

    aggregate(
      distance_from_centroid ~ type,
      data = subset(replicate_data, type %in% c("tip", "node")),
      FUN = mean
    ) |>
      transform(replicate = i)
  })
)

# centroid_distances_slow <- do.call(rbind, lapply(slow, distance.from.centroid))


# boxplot(
#   distance_from_centroid ~ type,
#   data = subset(centroid_distances_slow, type %in% c("tip", "node")),
#   xlab = "",
#   ylab = "Distance from shared centroid"
# )


##############################################################################


# THIS IS IMPORTANT BIT ####################################################

##############################################################################

matrices_true_list <- list()

for (i in 1:2){
  matrices_true_list[[i]] <- readRDS(paste0(sprintf("../../Data/revisions/wagner2/discrete/matrices/11532243_matrices_%03d.rds", i)))
}


distances_true <- lapply(matrices_true_list, lapply, char.diff, method = "hamming", by.col = FALSE)






distances_true_slow <- lapply(distances_true, function(x) x$slow)

node_dist_slow <- lapply(distances_true_slow, function(x){
  mat <- x[grepl("^n",rownames(x))]
  return(mean(mat))
})



tip_dist_slow <- lapply(distances_true_slow, function(x){
  mat <- x[grepl("^t",rownames(x))]
  return(mean(mat))
})
node_values <- unlist(node_dist_slow)
tip_values <- unlist(tip_dist_slow)

boxplot(
  c(node_values, tip_values) ~
    factor(rep(c("Node", "Tip"), each = length(node_values))),
  xlab = "",
  ylab = "Mean pairwise distance",
  col = c("firebrick", "steelblue")
)


disparity_vals_true <- lapply(
  distances_true_slow,mean
)



pre_ord_ace_list <- list()
for (i in 1:100){
  pre_ord_ace_list[[i]] <- readRDS(paste0(sprintf("../../Data/revisions/wagner2/discrete/anc/11532243_pre_ord_point_%03d.rds", i)))
}


pre_ord_ace_list <- list()
trees <- list()
for (i in 1:100){
  pre_ord_ace_list[[i]] <- readRDS(paste0(sprintf("../../Data/revisions/wagner/discrete/ord/11429393_ord_point_%03d.rds", i)))
}

centroid_distances <- do.call(
  rbind,
  lapply(seq_along(pre_ord_ace_list), function(i) {
    do.call(
      rbind,
      lapply(names(pre_ord_ace_list[[i]]), function(rate) {
        do.call(
          rbind,
          lapply(names(pre_ord_ace_list[[i]][[rate]]), function(fossil_level) {
            matrix <- pre_ord_ace_list[[i]][[rate]][[fossil_level]]

            data <- distance.from.centroid(matrix)
            data$transition_rate <- rate
            data$fossil_level <- fossil_level
            data$replicate <- i
            data
          })
        )
      })
    )
  })
)

centroid_distances <- subset(
  centroid_distances,
  type %in% c("tip", "node")
)

library(ggplot2)

centroid_distances$fossil_level <- factor(
  centroid_distances$fossil_level,
  levels = c("living", "fossil_low", "fossil_med", "fossil_high", "all")
)

centroid_distances$type <- factor(
  centroid_distances$type,
  levels = c("tip", "node")
)

tip_nodes_est_states <- ggplot(
  centroid_distances,
  aes(
    x = fossil_level,
    y = distance_from_centroid,
    fill = type
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  facet_wrap(~ transition_rate) +
  labs(
    x = "Fossil sampling",
    y = "Distance from centroid",
    fill = "Entity type"
  ) +
  scale_fill_manual(
    values = c(tip = "steelblue", node = "firebrick")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("tips_nodes_est_states.png", tip_nodes_est_states, )














distances_ace <- lapply(pre_ord_ace_list, lapply, lapply, char.diff, method = "hamming", by.col = FALSE)


distances_ace_slow <- lapply(distances_ace, function(x) x$slow)






ace_node_dist_slow <- lapply(distances_ace_slow, lapply, function(x){
  mat <- x[grepl("^n",rownames(x))]
  return(mean(mat))
})



ace_tip_dist_slow <- lapply(distances_ace_slow, lapply, function(x){
  mat <- x[grepl("^t",rownames(x))]
  return(mean(mat))
})


node_values <- unlist(ace_node_dist_slow)
tip_values <- unlist(ace_tip_dist_slow)

boxplot(
  c(node_values, tip_values) ~
    factor(rep(c("Node", "Tip"), each = length(node_values))),
  xlab = "",
  ylab = "Mean pairwise distance",
  col = c("firebrick", "steelblue")
)

disparity_vals_ace <- lapply(
  distances_ace_slow,lapply, mean
)



no_ace_list <- list()
for (i in 1:2){
  no_ace_list[[i]] <- readRDS(paste0(sprintf("../../Data/revisions/wagner2/discrete/matrices/11532243_fossil_matrices_%03d.rds", i)))
}


distances_no_ace <- lapply(no_ace_list, lapply, lapply, char.diff, method = "hamming", by.col = FALSE)




living_matrix <- no_ace_list[[1]]$slow$living$matrix
living_matrix_dist <- char.diff(living_matrix, method = "mord",  by.col = FALSE)
result <- cmdscale(living_matrix_dist, k = ncol(living_matrix_dist) - 1, add = TRUE)
result$ac   # the additive constant applied

corrected_D <- living_matrix_dist
corrected_D[upper.tri(corrected_D)] <- corrected_D[upper.tri(corrected_D)] + result$ac
sumsq_corrected <- sum(corrected_D[upper.tri(corrected_D)]^2)
sumsq_full <- sum(as.matrix(dist(living_matrix_ord))[upper.tri(as.matrix(dist(living_matrix_ord)))]^2)

# distances_no_ace_slow <- lapply(distances_no_ace, function(x) x$slow)

disparity_vals_no_ace <- lapply(
  distances_ace_slow,lapply, mean
)




## compare sum of variances across each
living_matrix_dist <- char.diff(living_matrix, method = "mord",  by.col = FALSE)
ordination <- cmdscale(living_matrix_dist, k = ncol(living_matrix_dist) - 1, add = TRUE)$points
sum(living_matrix_dist[upper.tri(living_matrix_dist)]^2) / (2 * (ncol(living_matrix_dist) ^2))
dispRity(ordination, metric = c(sum, variances))$disparity


### so we see a marked difference in the node vs tip disparity when it comes to the actual raw distance matrix, versus the ordination...



