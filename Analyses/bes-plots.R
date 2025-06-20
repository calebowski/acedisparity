library(treats)

trait <- make.traits(BM.process, trait.names = "Trait value")
bd_params <- make.bd.params(speciation = 1, extinction = 0.7)
stop_rule <- list(max.living = 20)

# tree <- rtree(n = )

trees <- treats(traits = trait, stop.rule = stop_rule, bd.params = bd_params, null.error = 100, replicates = 5)

## First plot is just tips and fossils, no edges or 
png(filename = "../besMacro/BMnoedge.png", width = 12, height = 10, units = "in", res = 300)
plot.treats(trees[[5]], cex = 3, cex.axis = 2.0, edge.width = 3, 
           legend.cex = 1.5, cex.lab = 3, col = c("nodes" = "NA"), 
           edges = NULL, legend = FALSE)

# Define colors for legend
# legend_colors <- c("red", "blue")  # or whatever colors you want
legend("topleft", col = c("nodes" = "orange", "fossils" = "lightblue", "livings" = "blue"), 
       legend = c("nodes", "fossils", "living"), 
       pch = 19, bg = "white", cex = 4, bty = "n")
dev.off()

png(filename = "../besMacro/BMnonode.png", width = 12, height = 10, units = "in", res = 300)
plot.treats(trees[[5]], legend = FALSE, cex = 3, cex.axis = 2.0, legend.cex = 1.5,edge.width = 15, cex.lab = 1.4, col = c("nodes" = "NA"), edges = "black", xlab = "", ylab = "")
dev.off()


png(filename = "../besMacro/BMfull.png", width = 12, height = 10, units = "in", res = 300)
plot.treats(trees[[5]], legend = FALSE, cex = 3, cex.axis = 2.0, legend.cex = 1.5,edge.width = 15, cex.lab = 1.4, edges = "black",xlab = "", ylab = "")
dev.off()


# plot.treats(trees[[2]], legend = TRUE, cex = 2, cex.axis = 1.2, legend.cex = 1.5, cex.lab = 1.4, col = c("nodes" = "NA"), edges = NULL)

discrete <- make.traits(discrete.process)
discrete_trees <- treats(traits = discrete, stop.rule = stop_rule, bd.params = bd_params, null.error = 100, replicates = 5)
plot.treats(discrete_trees[[1]], legend = TRUE, cex = 2, cex.axis = 1.2, legend.cex = 1.5, cex.lab = 1.4, tips.nodes = "blue")



###########################################################################################
# Make figures for different fossil sampling levels
bd_params <- make.bd.params(speciation = 1, extinction = 0.7)

stop_rule <- list(max.living = 8)
set.seed(123)
trees <- treats(stop.rule = stop_rule, bd.params = bd_params, null.error = 100, replicates = 3)
full_fossil_tree <- trees[[1]]
ages <- tree.age(full_fossil_tree)
living <- ages$element[ages$ages == 0]
living_tree <- keep.tip(full_fossil_tree, c(living))

original_root_height <- max(nodeHeights(full_fossil_tree))
# full_fossil_tree <- set.root.time(full_fossil_tree)
new_root_height <- max(nodeHeights(living_tree))

# Shift the pruned tree to match original height
height_diff <- original_root_height - new_root_height
living_tree$edge.length[living_tree$edge[,1] == Ntip(living_tree) + 1] <- 
living_tree$edge.length[living_tree$edge[,1] == Ntip(living_tree) + 1] + height_diff


png(filename = "../besMacro/fulltree.png", width = 10, height = 10, units = "in", res = 300, bg = "transparent")
plot.phylo(full_fossil_tree, edge.color = "darkorange", edge.width = 14, show.tip.label = FALSE)
dev.off()


png(filename = "../besMacro/livingtree.png", width = 10, height = 10, units = "in", res = 300, bg = "transparent")
plot.phylo(living_tree, edge.color = "darkorange", edge.width = 14, show.tip.label = FALSE)
dev.off()

## What if use a similar logic, but try and project the estimated values onto the real treats plot to compare how far they are.

living <- ages$element[ages$ages == 0]
# Or provide a manual vector of living taxa:
# living_tips <- c("tip1", "tip3", "tip7")

# Prune the tree
pruned_tree <- drop.tip(full_fossil_tree, setdiff(full_fossil_tree$tip.label, living))

# Now, to preserve the original root height:
# Shift the root node height back to match the full tree

# Compare root heights
original_root_height <- max(nodeHeights(full_fossil_tree))
# full_fossil_tree <- set.root.time(full_fossil_tree)
new_root_height <- max(nodeHeights(pruned_tree))

# Shift the pruned tree to match original height
height_diff <- original_root_height - new_root_height
pruned_tree$edge.length[pruned_tree$edge[,1] == Ntip(pruned_tree) + 1] <- 
pruned_tree$edge.length[pruned_tree$edge[,1] == Ntip(pruned_tree) + 1] + height_diff


library(ggplot2)
library(reshape2)

## Make discrete table

discrete_data <- data.frame(
  row.names = c("Taxon 1", "Taxon 2", "Taxon 3", "Taxon 4"),
  Trait1 = c(0, 1, 0, 1),
  Trait2 = c(1, 1, 0, 0),
  Trait3 = c(0, 0, 1, 1)
)

library(grid)

library(ggplot2)
library(gridExtra)

blue_orange_theme <- ttheme_default(
  core = list(
    fg_params = list(cex = 1.2, col = "navy"),                    
    bg_params = list(fill = c("lightblue", "peachpuff"),          
                     col = "white", lwd = 1)
  ),
  colhead = list(
    fg_params = list(cex = 1.4, fontface = "bold", col = "white"),
    bg_params = list(fill = "#4682B4", col = "white", lwd = 2)  
  ),
  rowhead = list(
    fg_params = list(cex = 1.2, fontface = "bold", col = "white"),
    bg_params = list(fill = "#FF8C00", col = "white", lwd = 2) 
  )
)

# Create the plot
p <- ggplot() + 
  annotation_custom(tableGrob(discrete_data, theme = blue_orange_theme)) + 
  theme_void()

ggsave("../besMacro/discrete_matrix.png", p, width = 8, height = 8, dpi = 100)


## Make continuous table

continuous_data <- data.frame(
  row.names = c("Taxon 1", "Taxon 2", "Taxon 3", "Taxon 4"),
  Trait1 = c(1.46, 2.31, 0.87, 1.23),
  Trait2 = c(3.12, 1.95, 2.67, 4.08),
  Trait3 = c(0.74, 1.89, 1.34, 0.52)
)

library(grid)

library(ggplot2)
library(gridExtra)

blue_orange_theme <- ttheme_default(
  core = list(
    fg_params = list(cex = 1.2, col = "#000080"),                    
    bg_params = list(fill = c("#ADD8E6", "#FFDAB9"),          
                     col = "white", lwd = 1)
  ),
  colhead = list(
    fg_params = list(cex = 1.4, fontface = "bold", col = "white"),
    bg_params = list(fill = "#4682B4", col = "white", lwd = 2)  
  ),
  rowhead = list(
    fg_params = list(cex = 1.2, fontface = "bold", col = "white"),
    bg_params = list(fill = "#FF8C00", col = "white", lwd = 2) 
  )
)

# Create the plot
p <- ggplot() + 
  annotation_custom(tableGrob(continuous_data, theme = blue_orange_theme)) + 
  theme_void()

ggsave("../besMacro/continuous_matrix.png", p, width = 8, height = 8, dpi = 100)

## Simulate BM process

png(filename = "../besMacro/bmprocess.png", width = 12, height = 9, units = "in", res = 600)
BM.process(x0 = 0)
plot(make.traits(process = BM.process), cex.axis = 1.4,cex.lab = 1.7)
dev.off()


## Simulate OU process

png(filename = "../besMacro/ouprocess.png", width = 12, height = 9, units = "in", res = 600)
OU.process(x0 = 0)
plot(make.traits(process = OU.process), cex.axis = 1.4,cex.lab = 1.7)
dev.off()



library(gridExtra)
library(ggplot2)

# Create your three matrices
high_transition <- matrix(c(0.1, 0.9, 0.9, 0.1), nrow = 2, byrow = TRUE)
medium_transition <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE)
low_transition <- matrix(c(0.99, 0.01, 0.01, 0.99), nrow = 2, byrow = TRUE)

# Set names
colnames(high_transition) <- rownames(high_transition) <- c("0", "1")
colnames(medium_transition) <- rownames(medium_transition) <- c("0", "1")
colnames(low_transition) <- rownames(low_transition) <- c("0", "1")

# Create ggplot objects
p1 <- ggplot() + 
  annotation_custom(tableGrob(high_transition, theme = blue_orange_theme)) + 
  theme_void() 

p2 <- ggplot() + 
  annotation_custom(tableGrob(medium_transition, theme = blue_orange_theme)) + 
  theme_void() 

p3 <- ggplot() + 
  annotation_custom(tableGrob(low_transition, theme = blue_orange_theme)) + 
  theme_void() 

# Arrange vertically
combined_plot <- arrangeGrob(p1, p2, p3, ncol = 1)
ggsave("../besMacro/transitions.png", combined_plot, width = 8, height = 15, dpi = 100)

medium_transition <- matrix(c(
  0.9, 0.1,
  0.1, 0.9
), nrow = 2, byrow = TRUE)

low_transition <- matrix(c(
  0.99, 0.01,
  0.01, 0.99
), nrow = 2, byrow = TRUE)


library(xtable)

# Your matrix
library(xtable)

# Your matrix
high_rate <- matrix(c(0.1, 0.9, 0.9, 0.1), nrow = 2, byrow = TRUE)
rownames(high_rate) <- colnames(high_rate) <- c("0", "1")
medium_rate <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE)
rownames(medium_rate) <- colnames(medium_rate) <- c("0", "1")
low_rate <- matrix(c(0.99, 0.01, 0.01, 0.99), nrow = 2, byrow = TRUE)


plot_matrix_brackets <- function(mat) {
  # Format matrix values with consistent spacing
  row1 <- paste(sprintf("%6.2f", mat[1,]), collapse = "   ")
  row2 <- paste(sprintf("%6.2f", mat[2,]), collapse = "   ")
  
  # Create bracket matrix string - exactly 2 lines
  bracket_text <- paste0("⎡ ", row1, " ⎤\n⎣ ", row2, " ⎦")
  
  ggplot() + 
    annotate("text", x = 0.5, y = 0.5, 
             label = bracket_text, 
             size = 8, family = "mono", hjust = 0.5, vjust = 0.5,
             color = "white") +  # White text
    theme_void() +
    theme(plot.margin = margin(30, 30, 30, 30),
          plot.background = element_rect(fill = "transparent", color = NA),  # Transparent background
          panel.background = element_rect(fill = "transparent", color = NA)) # Transparent panel
}

# Create plots for each
p1 <- plot_matrix_brackets(high_rate)
p2 <- plot_matrix_brackets(medium_rate)
p3 <- plot_matrix_brackets(low_rate)


# Arrange and save

ggsave("../besMacro/high_rate_matrix.png", p1, 
       width = 6, height = 4, dpi = 300, bg = "transparent")

ggsave("../besMacro/medium_rate_matrix.png", p2, 
       width = 6, height = 4, dpi = 300, bg = "transparent")

ggsave("../besMacro/low_rate_matrix.png", p3, 
       width = 6, height = 4, dpi = 300, bg = "transparent")


###################################################
## PIE Chart creation

library(ggplot2)
library(dplyr)

# Data
# data <- data.frame(
#   state = c(0, 1),
#   probability = c(0.4, 0.6)
# )

# data <- data %>%
#   arrange(desc(state)) %>%
#   mutate(
#     highlight = probability >= 0.5,
#     label = paste0("state ", state, " (", round(probability, 1), ")")
#   )
# ggplot(data, aes(x = factor(state), y = probability, fill = factor(state))) +
#   geom_col(width = 0.6) +
#   geom_hline(yintercept = relative_threshold, color = "black", size = 6, linetype = "dashed") +
#   scale_fill_manual(values = c("0" = "#FF8C00", "1" = "#FF8C00")) +
#   scale_y_continuous(limits = c(0, 0.8)) +
#   labs(x = "State", y = "") +
#   theme_void() +
#   theme(axis.text.x = element_text(size = 25, color = "white"),          # White x-axis text
#         axis.title.x = element_text(size = 25, color = "white", margin = margin(t = 10)), # White x-axis title
#         plot.background = element_rect(fill = "transparent", color = NA),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         legend.position = "none") +
#   # annotate("text", x = 2.3, y = relative_threshold + 0.03, 
#   #          label = paste0("Threshold = ", round(relative_threshold, 2)), 
#   #          color = "red", size = 5, fontface = "bold") +
#   geom_text(aes(label = probability), 
#             vjust = -0.3, size = 12, fontface = "bold", color = "white") 
# ggsave("../besMacro/strictmajoritypie.png", width = 6, height = 6, dpi = 300, bg = "transparent")


# # Your data
# # Your data
# # Your data
# data <- data.frame(
#   state = c(0, 1),
#   probability = c(0.4, 0.6)
# )

# Calculate relative majority threshold
k <- 2  
relative_threshold <- max(data$probability) - (1/k)

# Create super clean bar chart with no bar outlines
ggplot(data, aes(x = factor(state), y = probability, fill = factor(state))) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = relative_threshold, color = "black", size = 3, linetype = "dashed") +
  scale_fill_manual(values = c("0" = "#ADD8E6", "1" = "#FF8C00")) +
  scale_y_continuous(limits = c(0, 0.8)) +
  labs(x = "State", y = "") +
  theme_void() +
  theme(axis.text.x = element_text(size = 30, color = "white"),          # White x-axis text
        axis.title.x = element_text(size = 30, color = "white", margin = margin(t = 10)), # White x-axis title
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none") +
  # annotate("text", x = 2.3, y = relative_threshold + 0.03, 
  #          label = paste0("Threshold = ", round(relative_threshold, 2)), 
  #          color = "red", size = 5, fontface = "bold") +
  geom_text(aes(label = probability), 
            vjust = -0.3, size = 18, fontface = "bold", color = "white") 
# Save
ggsave("../besMacro/relativebarchart.png", width = 6, height = 6, bg = "transparent", dpi = 300)

# Save with transparent background

### strict threshold
strict_threshold <- 0.5


ggplot(data, aes(x = factor(state), y = probability, fill = factor(state))) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = strict_threshold, color = "black", size = 3, linetype = "dashed") +
  scale_fill_manual(values = c("0" = "#ADD8E6", "1" = "#FF8C00")) +
  scale_y_continuous(limits = c(0, 0.8)) +
  labs(x = "State", y = "") +
  theme_void() +
  theme(axis.text.x = element_text(size = 30, color = "white"),          # White x-axis text
        axis.title.x = element_text(size = 30, color = "white", margin = margin(t = 10)), # White x-axis title
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none") +
  # annotate("text", x = 2.3, y = relative_threshold + 0.03, 
  #          label = paste0("Threshold = ", round(relative_threshold, 2)), 
  #          color = "red", size = 5, fontface = "bold") +
  geom_text(aes(label = probability), 
            vjust = -0.3, size = 18, fontface = "bold", color = "white") 
# Save
ggsave("../besMacro/strictbarchart.png", width = 6, height = 6, bg = "transparent", dpi = 300)




prob_state_0 <- 0.4
prob_state_1 <- 0.6

# Create 10x10 grid data
grid_data <- expand.grid(x = 1:10, y = 1:10) %>%
  mutate(
    # Randomly assign states based on probabilities
    state = sample(c(0, 1), 100, replace = TRUE, prob = c(prob_state_0, prob_state_1)),
    # Or create exact proportions (40 state 0, 60 state 1)
    # state = c(rep(0, 40), rep(1, 60))
  )

# Create the grid plot
ggplot(grid_data, aes(x = x, y = y, fill = factor(state))) +
  geom_tile(color = "white", size = 0.5) +  # White gridlines
  scale_fill_manual(values = c("0" = "#ADD8E6", "1" = "#FF8C00")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA))

ggsave("../besMacro/samplechart.png", width = 6, height = 6, bg = "transparent", dpi = 300)
