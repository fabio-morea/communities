communities
================

<!-- badges: start -->

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13594209.svg)](https://doi.org/10.5281/zenodo.13594209)
[![Lifecycle:experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
![R](https://img.shields.io/badge/language-R-blue?logo=r)
<!-- badges: end -->

The `communities` package addresses critical reproducibility issues in network community detection by providing systematic tools to explore the
complete solution space of detection algorithms. Built on `igraph`, it implements a Bayesian framework for assessing algorithm stability,
managing input ordering bias, and identifying outliers through repeated trials.

A Python implementation of this package is available at [communities-py](https://github.com/fabio-morea/communities-py/)

# **Description**
## The Challenge
Community detection algorithms take a network as input and return a partition. This leaves researchers with the critical task of validation and interpretation—a task that becomes surprisingly difficult as networks grow in complexity.
Standard questions arise: Do the identified communities align with expectations? Are nodes genuinely more connected within their assigned community than across boundaries? And crucially, when running the algorithm multiple times yields different results, how should we interpret this variability?
The conventional approach treats such variability as algorithmic noise to be minimized or ignored—run once, select the most plausible result, and proceed.
We argue this approach discards valuable information.

## Our Approach
Solution Space Exploration: We conceptualize each algorithm execution as sampling a point within a broader solution space. Only through multiple runs can we map this space and understand the range of plausible community structures present in the network.
Statistical Convergence Criterion: The exploration process continues until reaching a statistically grounded stopping point. Each solution is evaluated with Bayesian confidence intervals, providing a principled basis for determining when sufficient information has been gathered.
Structural Taxonomy: The toolkit classifies solution spaces into five categories—Single, Dominant, Multiple, Sparse, or Empty—offering a formal characterization of community structure clarity and stability.
Mitigation of Input Ordering Effects: The network undergoes random permutation prior to each algorithm run, substantially reducing the dependence of results on arbitrary input data ordering.
Systematic Outlier Identification: A formal methodology identifies nodes that frequently transition between communities across multiple runs. Such nodes may represent structural bridges, peripheral members, or candidates for single-node community assignment.

## **Key Features**
- Solution space exploration with Bayesian convergence criteria
- Quality assessment using modularity, mixing parameters, and connectivity metrics
- Consensus community detection from multiple algorithm runs
- Outlier identification for ambiguous node assignments
- Benchmark network generation (ring-of-cliques variants) for algorithm testing
  
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

## **Community dertection algorithms supported**
| Code   | Algorithm           | Reference                  |
|--------|---------------------|----------------------------|
| `"IM"` | Infomap             | Rosvall & Bergstrom (2008) |
| `"LV"` | Louvain             | Blondel et al. (2008)      |
| `"LD"` | Leiden              | Traag et al. (2019)        |
| `"WT"` | Walktrap            | Pons & Latapy (2005)       |
| `"LP"` | Label Propagation   | Raghavan et al. (2007)     |
| `"EV"` | Leading Eigenvector | Newman (2006)              |
| `"EB"` | Edge Betweenness    | Girvan & Newman (2002)     |
 

# **Quick Start**
## **Installation**
``` r
install.packages("devtools")
devtools::install_github("fabio-morea/communities")
```
## Usage
``` r
rlibrary(communities)
library(igraph)

# Create test network with known structure
graph <- make_ring_of_cliques(n_cliques = 4, clique_size = 6)

# Explore solution space
solution_space <- solutions_space_DM(
  graph = graph,
  n_trials = 100,
  method = "IM",
  precision_threshold = 0.05
)

# Assess quality
quality <- assess_solution_quality(graph, solution_space)

# Build consensus communities
consensus <- consensus_communities(co_occurrence(solution_space), p = 0.5)

# Visualize results
plot_solutions(graph, solution_space, add_node_labels = TRUE)
```

## **Reporting Issues**

Please report bugs and feature requests through [GitHubIssues](https://github.com/fabio-morea/communities/issues). Include: -
Minimal reproducible example - Session information (`sessionInfo()`) -
Expected vs. actual behavior

# **License**
This project is licensed under the Creative Commons Attribution 4.0
International License (CC BY 4.0)

# **How to cite this work**

If you use this package in your research, please cite:

``` bibtex
@misc{communities2024,
  title = {A Comprehensive Framework for Solution Space Exploration in Community Detection},
  author = {Morea, Fabio and De Stefano, Domenico},
  year = {2024},
  url = {https://github.com/fabio-morea/communities}, 
  note = {R package version 0.1.0}
}
```

The theoretical framework underlying this package is described in:

``` bibtex
@article{morea2024comprehensive,
  title = {A Comprehensive Framework for Solution Space Exploration in Community Detection},
  author = {Morea, Fabio and De Stefano, Domenico},
  year = {2024},
  journal = {Under Review} 
}
```

For applications of solution space exploration and consensus methods in network analysis, see also:

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

## **Acknowledgments**

This package builds upon the excellent work of: - The **igraph** for foundational network analysis tools - The **aricode** for partition
similarity metrics - Community detection algorithm developers whose methods we implement

**Author**: Fabio Morea, Area Science Park, Trieste, Italy -
<fabio.morea@areasciencepark.it> this code was developed within the
framework of my PhD in Applied Data Science and Artificial Intelligence,
under the supervision of Prof. Domenico De Stefano, University of
Trieste, Italy.
