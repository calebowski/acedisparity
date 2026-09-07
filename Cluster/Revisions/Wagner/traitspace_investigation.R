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
  true = ord_true$med,
  no_ace = ord_no_ace$med$fossil_med,
  ace = ord_point$med$fossil_med,
  post_ord_ace = post_ord_point_fossil_med

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

  sampled_nodes <- grepl("^f_n", row_names)
  ordinary_tips <- grepl("^t", row_names)
  ordinary_nodes <- grepl("^n", row_names)

  plot(
    x[, 2], x[, 3],
    type = "n",
    main = title,
    xlab = "PCoA axis 1",
    ylab = "PCoA axis 2"
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
    "topright",
    legend = c(
      "Tips",
      "Sampled ancestral nodes",
      "Ancestral nodes"
    ),
    pch = 19,
    col = c("steelblue", "darkorange", "firebrick"),
    bty = "n",
    cex = 1.5
  )
}

par(mfrow = c(2, 2))

plot_space(spaces$true, "True", crown_tree)
plot_space(spaces$no_ace, "No ACE: fossil_med")
plot_space(spaces$ace, "ACE: fossil_med")
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




par(mfrow = c(1, 1))







### TIP TO PARENT TEST
ord_true_list <- list()
trees <- list()
for (i in 1:100){
  ord_true_list[[i]] <- readRDS(paste0(sprintf("../../Data/revisions/wagner/discrete/ord/11429393_ord_true_%03d.rds", i)))
  trees[[i]] <- extract.crown.tree(read.tree(sprintf("../../Data/trees/tree_50t_%03d.tre", i)))
}

get.parent.child <- function(tree) {
  tree_names <- c(tree$tip.label, tree$node.label)

  parent_child <- data.frame(
    parent = tree_names[tree$edge[, 1]],
    child = tree_names[tree$edge[, 2]],
    branch_length = tree$edge.length,
    morphological_distance = NA_real_,
    stringsAsFactors = FALSE
  )

  return(parent_child)

}


parent_children <- lapply(trees, get.parent.child)



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