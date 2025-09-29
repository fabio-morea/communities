communities
================

<!-- badges: start -->

[![License: CC BY
4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13594209.svg)](https://doi.org/10.5281/zenodo.13594209)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The `communities` package addresses critical reproducibility issues in
network community detection by providing systematic tools to explore the
complete solution space of detection algorithms. Built on `igraph`, it
implements a Bayesian framework for assessing algorithm stability,
managing input ordering bias, and identifying outlier nodes through
repeated trials with controlled randomization.

For detailed examples and usage, see the vignettes: - [Introduction to
communities](https://github.com/fabio-morea/communities/blob/main/vignette/getting-started.html) -
[Exploring the Solution
Space](https://github.com/fabio-morea/communities/blob/main/vignette/vignette-karate.html)

## **Key Features**

### **Solution Space Exploration**

- **Systematic exploration**: Address multiplicity of solutions through
  repeated algorithm runs with vertex permutation
- **Bayesian convergence**: Model solution probabilities using
  Dirichlet-Multinomial framework with defensible stopping rules
- **Input ordering bias mitigation**: Systematic permutation of node
  orderings to ensure topology-dependent results
- **Solution space taxonomy**: Classify solution spaces as Single,
  Dominant, Multiple, Sparse, or Empty

### **Outlier Detection and Quality Assessment**

- **Node outlier identification**: Data-driven approach using pairwise
  agreement matrices (Γ) from solution space exploration
- **Quality metrics**: Modularity, empirical mixing parameter (μ),
  internal connectivity assessment
- **Validity criteria**: Systematic filtering of degenerate solutions
  (trivial partitions, disconnected communities)

### **Advanced Visualization**

- **Solution space plots**: Visualize probability distributions and
  convergence patterns
- **Community networks**: Generate meta-networks representing community
  relationships
- **Consensus visualization**: Display agreement patterns across
  multiple solutions
- **Custom layouts**: Community-aware graph layouts for better
  interpretation

### **Benchmark Networks**

- **Ring-of-cliques generator**: Parameterizable test networks (RC,
  RC+B, RC+C, RC+BC) with known ground truth
- **Controlled complexity**: Systematic testing of algorithm behavior
  across network structures
- **Outlier scenarios**: Networks with central nodes to test outlier
  detection capabilities

## **Installation**

### Development Version (Recommended)

``` r
# Install from GitHub using pak (recommended)
install.packages("pak")
pak::pak("fabio-morea/communities")

# Alternative: using devtools
install.packages("devtools")
devtools::install_github("fabio-morea/communities")
```

### Dependencies

The package requires: - **R** ≥ 4.0.0 - **igraph** ≥ 1.3.0 (network
analysis) - **aricode** (partition similarity metrics) - **dplyr** (data
manipulation) - **ggplot2** (visualization) - **tibble** (modern data
frames)

## **Quick Start**

### Basic Solution Space Analysis

``` r
library(communities)
library(igraph)

# Create a test network with known community structure
graph <- make_ring_of_cliques(n_cliques = 4, clique_size = 6)

# Explore the solution space using Infomap
solution_space <- solutions_space_DM(
  graph = graph,
  n_trials = 100,
  method = "IM",
  precision_threshold = 0.05,
  verbose = TRUE
)

# View solution probabilities
print(solution_space$probabilities)
#>    id count  phat  pmin  pmax
#> 1 s01    78 0.780 0.682 0.862
#> 2 s02    15 0.150 0.086 0.234
#> 3 s03     7 0.070 0.028 0.139

# Classify the solution space
space_type <- solution_space_type(solution_space)
print(space_type)  # "Dominant"
```

### Quality Assessment

``` r
# Assess solution quality with multiple metrics
quality_results <- assess_solution_quality(
  graph = graph,
  solution_space = solution_space,
  mu_threshold = 0.5,
  require_connected_communities = TRUE
)

print(quality_results)
#> # A tibble: 3 × 7
#>   solution_id n_communities modularity mixing_parameter internally_connected is_valid invalidity_reason
#>         <int>         <int>      <dbl>            <dbl>                <lgl>    <lgl> <chr>           
#> 1           1             4      0.652            0.125                 TRUE     TRUE ""              
#> 2           2             3      0.543            0.234                 TRUE     TRUE ""              
#> 3           3             5      0.423            0.456                 TRUE     TRUE ""
```

### Consensus Community Detection

``` r
# Build co-occurrence matrix
co_occurrence_matrix <- co_occurrence(solution_space)

# Find consensus communities
consensus_result <- consensus_communities(
  co_occurrence_matrix, 
  p = 0.5,                    # Threshold for community membership
  group_outliers = FALSE
)

print(consensus_result)
#>    name cons_comm_label gamma comm_size single
#> 1   c1_1               1 0.89         6  FALSE
#> 2   c1_2               1 0.92         6  FALSE
#> ...
```

### Visualization

``` r
# Plot solution space evolution and statistics
plots <- plot_sol_space(solution_space, quality_results)

# Display probability evolution over trials
plot_sol_space_evolution(solution_space, show_ci = TRUE)

# Visualize all solutions found
plot_solutions(graph, solution_space, 
               add_node_labels = TRUE,
               add_prob_labels = TRUE)
```

## **Advanced Usage**

### Custom Network Generation

``` r
# Create different ring-of-cliques variants
rc_basic <- make_ring_of_cliques(4, 5, variant = "RC")       # Basic ring
rc_bridged <- make_ring_of_cliques(4, 5, variant = "RC+B")   # With bridge nodes  
rc_central <- make_ring_of_cliques(4, 5, variant = "RC+C")   # With central hub
rc_both <- make_ring_of_cliques(4, 5, variant = "RC+BC")     # Bridges + center

# Examine network properties
empirical_mu(rc_basic, V(rc_basic)$gt_community)
internally_connected(rc_basic, V(rc_basic)$gt_community)
```

### Algorithm Comparison

``` r
# Compare multiple detection methods
methods <- c("IM", "LV", "LD", "WT")
results <- map(methods, ~solutions_space_DM(graph, method = .x, n_trials = 50))
names(results) <- methods

# Analyze method stability
stabilities <- map_dbl(results, ~{
  max_prob <- max(.$probabilities$phat)
  return(max_prob)
})
print(stabilities)
```

### Community-Based Network Analysis

``` r
# Create community-level network
V(graph)$community <- consensus_result$cons_comm_label
community_network <- make_community_network(graph)

# Analyze community relationships
plot(community_network, 
     vertex.size = V(community_network)$size * 5,
     vertex.label = V(community_network)$name,
     edge.width = E(community_network)$weight)
```

## **Core Functions**

| Function | Purpose | Key Parameters |
|----|----|----|
| `solutions_space_DM()` | Explore solution space with Bayesian framework | `method`, `n_trials`, `precision_threshold` |
| `assess_solution_quality()` | Comprehensive quality assessment | `mu_threshold`, `min_communities` |
| `consensus_communities()` | Build consensus from multiple solutions | `p` (threshold), `group_outliers` |
| `co_occurrence()` | Compute node co-occurrence probabilities | \- |
| `plot_solutions()` | Visualize all detected solutions | `add_node_labels`, `add_prob_labels` |
| `make_ring_of_cliques()` | Generate benchmark networks | `n_cliques`, `clique_size`, `variant` |
| `solution_space_type()` | Classify solution space structure | \- |

## **Method Support**

| Code   | Algorithm           | Reference                  |
|--------|---------------------|----------------------------|
| `"IM"` | Infomap             | Rosvall & Bergstrom (2008) |
| `"LV"` | Louvain             | Blondel et al. (2008)      |
| `"LD"` | Leiden              | Traag et al. (2019)        |
| `"WT"` | Walktrap            | Pons & Latapy (2005)       |
| `"LP"` | Label Propagation   | Raghavan et al. (2007)     |
| `"EV"` | Leading Eigenvector | Newman (2006)              |
| `"EB"` | Edge Betweenness    | Girvan & Newman (2002)     |

## **Theoretical Background**

The package implements several methodological innovations:

### **Solution Space Exploration**

Community detection returns a single partition ad each run, but the
result may vary across different runs. This reflects a feature of the
network (rather than a bug in the algorithm): networks may admit
multiple equally valid community structures. This package systematically
explores this solution space by:

1.  **Controlled randomization**: Using vertex permutation to generate
    algorithm diversity
2.  **Bayesian modeling**: Treating solutions as draws from a
    multinomial distribution
3.  **Convergence assessment**: Using credible intervals to determine
    when sufficient exploration has occurred

### **Quality Assessment Framework**

The package provides multiple complementary quality metrics:

- **Modularity** (Q): Measures strength of community structure
- **Mixing parameter** (μ): Proportion of inter-community edges  
- **Internal connectivity**: Whether communities form connected
  subgraphs
- **Stability**: Consistency across multiple algorithm runs

### **Consensus Detection**

When multiple solutions exist, consensus methods aggregate information:

1.  **Co-occurrence matrix**: Probability that node pairs belong to same
    community
2.  **Threshold-based clustering**: Group nodes with high co-occurrence
3.  **Uncertainty quantification**: Measure confidence in community
    assignments

## **Use Cases**

### **Research Applications**

- **Algorithm evaluation**: Systematic assessment of community detection
  method stability and reliability  
- **Reproducible analysis**: Address the “single run fallacy” in network
  studies through principled exploration
- **Outlier identification**: Data-driven detection of nodes that don’t
  clearly belong to any community
- **Input bias assessment**: Quantify sensitivity to node ordering and
  other implementation artifacts

### **Practical Applications**

- **Collaboration networks**: Identify stable research communities with
  confidence measures (as demonstrated on EU Horizon projects)
- **Social network analysis**: Robust community detection in social
  platforms with uncertainty quantification
- **Biological networks**: Find reliable functional modules while
  handling ambiguous node assignments  
- **Infrastructure networks**: Detect communities robust to data
  collection artifacts and processing variations

### **Reporting Issues**

Please report bugs and feature requests through [GitHub
Issues](https://github.com/fabio-morea/communities/issues). Include: -
Minimal reproducible example - Session information (`sessionInfo()`) -
Expected vs. actual behavior

## **License**

This project is licensed under the Creative Commons Attribution 4.0
International License (CC BY 4.0)

## **How to cite this work**

If you use this package in your research, please cite:

``` bibtex
@misc{communities2024,
  title = {communities: A Comprehensive Framework for Solution Space 
           Exploration in Community Detection},
  author = {Morea, Fabio and De Stefano, Domenico},
  year = {2024},
  url = {https://github.com/fabio-morea/communities}, 
  note = {R package version 0.1.0}
}
```

The theoretical framework underlying this package is described in:

``` bibtex
@article{morea2024comprehensive,
  title = {A Comprehensive Framework for Solution Space Exploration in 
           Community Detection},
  author = {Morea, Fabio and De Stefano, Domenico},
  year = {2024},
  journal = {Under Review},
  note = {Preprint available at: https://zenodo.org/records/13825471}
}
```

## **related Work**

For applications of solution space exploration and consensus methods in
network analysis, see:

``` bibtex
@article{morea2025mapping,
  author = {Morea, F. and Soraci, A. and De Stefano, D.},
  title = {Mapping leadership and communities in EU-funded research through 
           network analysis},
  journal = {Open Research Europe},
  volume = {4},
  number = {268},
  year = {2025},
  doi = {10.12688/openreseurope.18544.2},
  note = {Version 2; peer review: 1 approved, 2 approved with reservations}
}
```

## 🙏 **Acknowledgments**

This package builds upon the excellent work of: - The **igraph** for
foundational network analysis tools - The **aricode** for partition
similarity metrics  
- Community detection algorithm developers whose methods we implement

**Author**: Fabio Morea, Area Science Park, Trieste, Italy -
<fabio.morea@areasciencepark.it> this code was developed within the
framework of my PhD in Applied Data Science and Artificial Intelligence,
under the supervizion of Prof. Domenico De Stefano, University of
Trieste, Italy.
