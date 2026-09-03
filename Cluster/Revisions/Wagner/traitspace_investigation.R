library(dispRity)


tree <- read.tree("../../Data/trees/tree_50t_001.tre")

ages <- tree.age(tree) # get tip ages
extant <- ages$element[ages$ages == 0]
living_tree <- keep.tip(tree, extant) ## root at mrca - alike to crown group analyses.

crown_tree  <- extract.clade(tree, living_tree$node.label[1]) ## this is the crown tree

fossil_matrices <- readRDS("..//Data/revisions/wagner/discrete/matrices/11429393_fossil_matrices_001.rds")
ord_no_ace <- readRDS("..//Data/revisions/wagner/discrete/ord/11429393_ord_no_ace_001.rds")

ord_true <- readRDS("..//Data/revisions/wagner/discrete/ord/11429393_ord_true_001.rds")

ord_point <- readRDS("..//Data/revisions/wagner/discrete/ord/11429393_ord_point_001.rds")


par(mfrow = c(2,2))
plot(ord_no_ace$slow$fossil_high)
plot(ord_point$slow$fossil_high)
plot(ord_true$slow)


ords <- list(true = ord_true, no_ace = ord_no_ace, ace = ord_point)


# get.disparity(dispRity(ord_no_ace$slow$fossil_high, metric = c(sum, variances)))
# get.disparity(dispRity(ord_point$slow$fossil_high, metric = c(sum, variances)))
# get.disparity(dispRity(ord_true$slow, metric = c(sum, variances)))


# lapply(ords, function(x) get.disparity(dispRity(x, metric = c(sum, variances))))

# Use the same two axes and limits for every panel
spaces <- list(
  true = ord_true$med,
  no_ace = ord_no_ace$med$fossil_med,
  ace = ord_point$med$fossil_med
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
    x[, 1], x[, 2],
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

par(mfrow = c(1, 3))

plot_space(spaces$true, "True", crown_tree)
plot_space(spaces$no_ace, "No ACE: fossil_med")
plot_space(spaces$ace, "ACE: fossil_med")

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





## show nodes that descended from another node, vs nodes that descended from tips. this can prove the point that derived nodes, clsoer to end of tree, are important to sample because they fill in the gaps.