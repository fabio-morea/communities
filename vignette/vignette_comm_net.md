---
title: "Building Community Networks with make_community_network()"
author: "communities package"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Building Community Networks}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)
```

## Introduction

The `make_community_network()` function creates a meta-network representation where each node represents a community from the original network. This is useful for:

- Understanding relationships between communities
- Visualizing high-level network structure
- Analyzing inter-community connectivity patterns
- Simplifying large networks for easier interpretation

## Example: Zachary's Karate Club

We'll use the famous Zachary's Karate Club network, which shows social ties between members of a university karate club. The network famously split into two groups after a dispute.

### Load Required Libraries

```{r libraries, message=FALSE}
library(igraph)
library(communities)  # Our package
library(dplyr)
```

### Load and Explore the Network

```{r load_data}
# Load Zachary's karate club network
g <- make_graph("Zachary")

# Basic network information
cat("Number of nodes:", vcount(g), "\n")
cat("Number of edges:", ecount(g), "\n")
cat("Network density:", edge_density(g), "\n")
```

### Detect Communities

First, we need to detect communities in the original network. We'll use the Louvain algorithm, which is popular for finding modular structure.

```{r detect_communities}
# Detect communities using Louvain algorithm
louvain_comms <- cluster_louvain(g)

# Assign community labels to vertices
V(g)$community <- membership(louvain_comms)

# Set edge weights (required by make_community_network)
E(g)$w <- 1.0  # Unweighted network

# Summary of detected communities
cat("Number of communities:", max(V(g)$community), "\n")
cat("Modularity:", modularity(louvain_comms), "\n")
cat("\nCommunity sizes:\n")
table(V(g)$community)
```

### Visualize Original Network with Communities

```{r plot_original, fig.width=8, fig.height=6}
# Set colors for communities
community_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
V(g)$color <- community_colors[V(g)$community]

# Plot original network
set.seed(42)
plot(g,
     vertex.label = NA,
     vertex.size = 8,
     vertex.color = V(g)$color,
     edge.color = "gray70",
     main = "Zachary's Karate Club Network\nColored by Community",
     layout = layout_with_fr(g))

# Add legend
legend("bottomleft",
       legend = paste("Community", 1:max(V(g)$community)),
       col = community_colors[1:max(V(g)$community)],
       pch = 19,
       bty = "n",
       cex = 0.8)
```

### Build the Community Network

Now we create the meta-network where each node is a community:

```{r build_comm_network}
# Create community network
gc <- make_community_network(g)

# Explore the community network
cat("Community network:\n")
cat("Number of nodes (communities):", vcount(gc), "\n")
cat("Number of edges (inter-community links):", ecount(gc), "\n")
cat("\nCommunity sizes:\n")
print(data.frame(
  community = V(gc)$name,
  size = V(gc)$size,
  percentage = round(100 * V(gc)$size / vcount(g), 1)
))

cat("\nEdge weights (connections between communities):\n")
print(as_data_frame(gc, what = "edges"))
```

### Visualize the Community Network

```{r plot_comm_network, fig.width=7, fig.height=6}
# Set node colors to match original
V(gc)$color <- community_colors[as.numeric(V(gc)$name)]

# Plot community network
# Node size proportional to community size
set.seed(123)
plot(gc,
     vertex.label = paste0("C", V(gc)$name, "\n(n=", V(gc)$size, ")"),
     vertex.size = sqrt(V(gc)$size) * 10,  # Scale by sqrt for better visualization
     vertex.color = V(gc)$color,
     vertex.label.color = "black",
     vertex.label.cex = 0.9,
     edge.width = E(gc)$weight / 5,  # Edge width proportional to connection strength
     edge.color = "gray50",
     edge.label = E(gc)$weight,
     edge.label.cex = 0.8,
     main = "Community Network\n(Node size = community size, Edge width = connection strength)",
     layout = layout_nicely(gc))
```

### Analyze Inter-Community Connectivity

Let's examine the connectivity patterns:

```{r analyze_connectivity}
# Create adjacency matrix for community network
adj_matrix <- as_adjacency_matrix(gc, attr = "weight", sparse = FALSE)
rownames(adj_matrix) <- paste("Community", rownames(adj_matrix))
colnames(adj_matrix) <- paste("Community", colnames(adj_matrix))

cat("Inter-community connectivity matrix:\n")
print(adj_matrix)

# Calculate connectivity metrics
edge_df <- as_data_frame(gc, what = "edges")

cat("\n\nInter-community connectivity summary:\n")
cat("Total inter-community edges:", sum(edge_df$weight), "\n")
cat("Average edge weight:", round(mean(edge_df$weight), 2), "\n")
cat("Max edge weight:", max(edge_df$weight), "\n")
cat("Min edge weight:", min(edge_df$weight), "\n")

# Identify strongest connection
strongest <- edge_df[which.max(edge_df$weight), ]
cat("\nStrongest inter-community link:\n")
cat(sprintf("  Community %s <-> Community %s (weight: %d)\n",
            strongest$from, strongest$to, strongest$weight))
```

### Compare with Ground Truth

The karate club historically split into two factions. Let's see if our communities capture this:

```{r ground_truth}
# The network has a known split (captured in vertex attributes if available)
# Community 1 and 2 are typically the two main factions

# Calculate mixing between communities
cat("Community separation analysis:\n\n")

for (i in 1:max(V(g)$community)) {
  # Internal edges
  subg <- induced_subgraph(g, V(g)[V(g)$community == i])
  internal_edges <- ecount(subg)
  
  # All edges involving this community
  all_edges <- sum(E(g)$w[which(V(g)$community[ends(g, E(g))[,1]] == i | 
                                 V(g)$community[ends(g, E(g))[,2]] == i)])
  
  # External edges
  external_edges <- all_edges - internal_edges
  
  cat(sprintf("Community %d:\n", i))
  cat(sprintf("  Internal edges: %d\n", internal_edges))
  cat(sprintf("  External edges: %d\n", external_edges))
  cat(sprintf("  Internal ratio: %.2f%%\n\n", 
              100 * internal_edges / (internal_edges + external_edges)))
}
```

### Practical Applications

The community network representation is useful for:

1. **Hierarchical Analysis**: Understanding multi-scale structure
2. **Boundary Analysis**: Identifying bridges between communities
3. **Coarse-Graining**: Simplifying large networks
4. **Communication Patterns**: Analyzing information flow between groups
5. **Intervention Planning**: Identifying key inter-community links

### Example: Finding Bridge Nodes

```{r bridge_nodes}
# Find nodes that connect different communities (high betweenness)
V(g)$betweenness <- betweenness(g)

# Identify top bridge nodes
bridge_nodes <- order(V(g)$betweenness, decreasing = TRUE)[1:5]

cat("Top 5 bridge nodes (high betweenness):\n")
for (i in bridge_nodes) {
  cat(sprintf("  Node %d: Community %d, Betweenness = %.1f\n",
              i, V(g)$community[i], V(g)$betweenness[i]))
}

# These nodes are important for inter-community communication
```

## Advanced: Weighted Networks

For weighted networks, edge weights are aggregated between communities:

```{r weighted_example, eval=FALSE}
# Example with weighted network
g_weighted <- g
E(g_weighted)$w <- runif(ecount(g_weighted), 0.1, 1.0)  # Random weights

gc_weighted <- make_community_network(g_weighted)

# Edge weights now represent total connection strength
print(as_data_frame(gc_weighted, what = "edges"))
```

## Summary

The `make_community_network()` function provides a powerful way to:

- Reduce complex networks to their community-level structure
- Visualize relationships between groups
- Analyze inter-community connectivity patterns
- Identify important boundary regions

This is particularly valuable for large networks where individual node-level analysis becomes unwieldy.

## References

- Zachary, W. W. (1977). An information flow model for conflict and fission in small groups. *Journal of Anthropological Research*, 33(4), 452-473.
- Blondel, V. D., et al. (2008). Fast unfolding of communities in large networks. *Journal of Statistical Mechanics: Theory and Experiment*, 2008(10), P10008.

## Session Info

```{r session_info}
sessionInfo()
```
