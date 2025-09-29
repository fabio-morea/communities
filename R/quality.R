#' Calculate the Empirical Mixing Parameter (Mu)
#'
#' Calculates the mixing parameter μ (mu) for a given network, which represents 
#' the proportion of edges that connect nodes from different communities. This 
#' parameter quantifies how mixed or modular a network is.
#'
#' @param graph An igraph object representing the network. Edge weights will be 
#'   set to 1.0 if not present.
#' @param community_labels A vector of community labels, ordered according to 
#'   the vertices in \code{graph}. These labels define the community membership 
#'   of each node.
#'
#' @return A numeric value in [0, 1] representing the empirical mixing parameter 
#'   (μ), which is the ratio of inter-community edge weights to total edge weights.
#'   Higher values indicate more mixing between communities.
#'
#' @details 
#' The function ensures all edges have weights (defaulting to 1.0 if missing), 
#' then calculates the proportion of edge weight that connects nodes from 
#' different communities relative to the total edge weight in the graph.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' graph <- make_ring_of_cliques(3, 5)
#' community_labels <- V(graph)$gt_community
#' mu <- empirical_mu(graph, community_labels)
#' print(mu)
#' }
#'
#' @export
empirical_mu <- function(graph, community_labels) {
    
    # Ensure all edges have weights (default to 1.0)
    igraph::E(graph)$weight <- 1.0
    
    # Assign community labels to vertices
    igraph::V(graph)$community <- community_labels
    
    # Convert graph to long-format edge data frame
    edge_data <- igraph::as_long_data_frame(graph)
    
    # Identify inter-community edges
    edge_data$is_inter_community <- (edge_data$from_community != edge_data$to_community)
    
    # Calculate total weight of inter-community edges
    inter_community_weight <- sum(edge_data$weight[edge_data$is_inter_community == TRUE])
    
    # Calculate mixing parameter (ratio of inter-community to total weight)
    mixing_parameter <- inter_community_weight / sum(edge_data$weight)
    
    return(mixing_parameter)
}

#' Calculate the Number of Connected Components within Each Community
#'
#' Calculates how many connected components exist within each community of a 
#' given network. For each community, it extracts the induced subgraph and 
#' computes the number of connected components.
#'
#' @param graph An igraph object representing the network to be analyzed.
#' @param community_labels A vector of community labels, where each label 
#'   corresponds to a node in the graph. The length of this vector must match 
#'   the number of vertices in the graph.
#'
#' @return A numeric vector where each element corresponds to the number of 
#'   connected components in one of the communities. The order matches the 
#'   order of unique community labels in the input.
#'
#' @details 
#' The function assigns community labels to vertices, then for each community, 
#' extracts the induced subgraph and counts connected components. A community 
#' is internally connected if it has exactly one connected component.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' graph <- make_ring_of_cliques(4, 5)
#' community_labels <- V(graph)$gt_community
#' internal_components <- internally_connected(graph, community_labels)
#' print(internal_components)
#' }
#'
#' @export
internally_connected <- function(graph, community_labels) {
    
    # Validate input lengths match
    stopifnot(length(community_labels) == igraph::vcount(graph))
    
    # Assign community labels to vertices
    igraph::V(graph)$community <- community_labels
    
    # Get unique community identifiers
    unique_communities <- unique(igraph::V(graph)$community)
    
    # Initialize result vector
    n_components_per_community <- c()
    
    # Process each community
    for (community_id in unique_communities) {
        
        # Extract induced subgraph for this community
        community_subgraph <- igraph::induced_subgraph(
            graph, 
            igraph::V(graph)$community == community_id
        )
        
        # Count connected components in this community
        n_components <- igraph::components(community_subgraph)$no
        
        # Append to results
        n_components_per_community <- c(n_components_per_community, n_components)
    }
    
    return(n_components_per_community)
}

#' Perform Quality Checks on Community Detection Solutions
#'
#' Evaluates the quality of community detection solutions by computing multiple 
#' metrics including number of communities, modularity, mixing parameter, and 
#' internal connectivity. Flags solutions as valid or invalid based on quality 
#' criteria.
#'
#' @param graph An igraph object representing the network to be analyzed.
#' @param solution_space_result A list object from \code{solutions_space_DM()} 
#'   containing \code{partitions} (matrix of community labels for each solution).
#' @param mu_max Numeric. Maximum acceptable mixing parameter value (default: 0.5). 
#'   Solutions with μ > mu_max are flagged as invalid.
#'
#' @return A tibble with one row per solution containing:
#' \describe{
#'   \item{\code{solution_id}}{Solution identifier (1 to n_solutions)}
#'   \item{\code{k}}{Number of unique communities}
#'   \item{\code{modularity}}{Newman's modularity score}
#'   \item{\code{mu}}{Empirical mixing parameter}
#'   \item{\code{int_conn}}{Whether all communities are internally connected}
#'   \item{\code{valid}}{Overall validity flag (TRUE/FALSE)}
#'   \item{\code{reason}}{Explanation when valid = FALSE}
#' }
#'
#' @details 
#' Quality assessment includes:
#' \itemize{
#'   \item{Number of communities (k)}
#'   \item{Modularity score (higher is better)}
#'   \item{Mixing parameter μ (should be ≤ mu_max)}
#'   \item{Internal connectivity (all communities should be connected)}
#' }
#'
#' Validity criteria (all must be satisfied):
#' \itemize{
#'   \item{μ ≤ mu_max}
#'   \item{k > 1 (more than one community)}
#'   \item{All communities are internally connected}
#'   \item{No missing metric values}
#' }
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' graph <- make_ring_of_cliques(4, 5)
#' solution_space_result <- solutions_space_DM(graph, n_trials = 50)
#' quality_results <- quality_check(graph, solution_space_result)
#' print(quality_results)
#' }
#'
#' @export
quality_check <- function(graph, solution_space_result, mu_max = 0.5) {
    
    # ---- Input validation ----
    stopifnot(inherits(graph, "igraph"))
    
    if (is.null(solution_space_result$partitions)) {
        stop("solution_space_result$partitions is missing", call. = FALSE)
    }
    
    # ---- Extract and standardize partition matrix ----
    partition_matrix <- solution_space_result$partitions
    
    # Handle single solution (vector) by converting to matrix
    if (is.null(ncol(partition_matrix))) {
        partition_matrix <- matrix(partition_matrix, ncol = 1)
    }
    
    n_solutions <- ncol(partition_matrix)
    n_vertices <- igraph::vcount(graph)
    
    # ---- Initialize result vectors ----
    n_communities <- numeric(n_solutions)
    modularity_scores <- numeric(n_solutions)
    mixing_parameters <- numeric(n_solutions)
    internal_connectivity <- logical(n_solutions)
    validity_flags <- rep(TRUE, n_solutions)
    invalidity_reasons <- character(n_solutions)
    
    # ---- Helper function to append invalidity reasons ----
    append_reason <- function(solution_idx, reason_text) {
        if (is.na(invalidity_reasons[solution_idx]) || 
            invalidity_reasons[solution_idx] == "") {
            invalidity_reasons[solution_idx] <<- reason_text
        } else {
            invalidity_reasons[solution_idx] <<- paste(
                invalidity_reasons[solution_idx], 
                reason_text, 
                sep = " | "
            )
        }
    }
    
    # ---- Assess each solution ----
    for (solution_idx in seq_len(n_solutions)) {
        
        # Extract raw community labels for this solution
        raw_labels <- partition_matrix[, solution_idx]
        
        # Normalize labels to consecutive integers 1..k
        # This avoids issues with sparse labeling (e.g., labels 0, 5, 7)
        normalized_labels <- as.integer(
            factor(raw_labels, levels = unique(raw_labels))
        )
        
        # Compute quality metrics
        
        # 1) Number of communities
        n_communities[solution_idx] <- dplyr::n_distinct(normalized_labels)
        
        # 2) Modularity score
        modularity_scores[solution_idx] <- igraph::modularity(graph, normalized_labels)
        
        # 3) Mixing parameter μ
        mixing_parameters[solution_idx] <- empirical_mu(graph, normalized_labels)
        
        # 4) Internal connectivity check
        connectivity_vector <- internally_connected(graph, normalized_labels)
        internal_connectivity[solution_idx] <- all(
            as.logical(connectivity_vector), 
            na.rm = TRUE
        )
        
        # Apply validity rules (with NA-safe checks)
        
        # Rule 1: Mixing parameter threshold
        if (!is.na(mixing_parameters[solution_idx]) && 
            mixing_parameters[solution_idx] > mu_max) {
            validity_flags[solution_idx] <- FALSE
            append_reason(
                solution_idx, 
                sprintf("mu=%.3f > %.3f", mixing_parameters[solution_idx], mu_max)
            )
        }
        
        # Rule 2: Trivial partition (single community)
        if (!is.na(n_communities[solution_idx]) && 
            n_communities[solution_idx] == 1) {
            validity_flags[solution_idx] <- FALSE
            append_reason(solution_idx, "k=1 (trivial partition)")
        }
        
        # Rule 3: Internal connectivity requirement
        if (!is.na(internal_connectivity[solution_idx]) && 
            !internal_connectivity[solution_idx]) {
            validity_flags[solution_idx] <- FALSE
            append_reason(solution_idx, "communities not internally connected")
        }
        
        # Rule 4: Missing metrics
        if (any(is.na(c(n_communities[solution_idx], 
                        modularity_scores[solution_idx], 
                        mixing_parameters[solution_idx], 
                        internal_connectivity[solution_idx])))) {
            validity_flags[solution_idx] <- FALSE
            append_reason(solution_idx, "NA in one or more metrics")
        }
    }
    
    # ---- Return results as tibble ----
    tibble::tibble(
        solution_id = seq_len(n_solutions),
        k = n_communities,
        modularity = modularity_scores,
        mu = mixing_parameters,
        int_conn = internal_connectivity,
        valid = validity_flags,
        reason = invalidity_reasons
    )
}

#' Compute Similarity Matrix Between Community Partitions Using NMI
#'
#' Computes a pairwise similarity matrix between community detection solutions
#' using Normalized Mutual Information (NMI) with square-root normalization.
#' The measure is invariant to community label permutations.
#'
#' @param partitions A tibble, data.frame, or matrix with n rows (nodes)
#'   and k columns (solutions). Each column contains community labels for
#'   one partitioning of the same graph.
#'
#' @return A k × k symmetric numeric matrix with values in [0, 1]:
#' \itemize{
#'   \item{\code{1}}{Identical partitions}
#'   \item{\code{0}}{Completely dissimilar partitions}
#' }
#' Row and column names correspond to the input column names (or
#' \code{"sol1"}, \code{"sol2"}, ... if missing).
#'
#' @details
#' This function uses \code{\link[aricode]{NMI}} from the aricode package
#' with \code{variant = "sqrt"}, which ensures values are in [0, 1].
#' The diagonal is always set to 1 (perfect self-similarity).
#'
#' @examples
#' \dontrun{
#' # Compare multiple partitions
#' partition_matrix <- cbind(
#'   sol1 = c(1, 1, 2, 2, 3, 3),
#'   sol2 = c(1, 1, 1, 2, 2, 2),
#'   sol3 = c(1, 2, 1, 2, 1, 2)
#' )
#' similarity <- similarity_matrix_nmi(partition_matrix)
#' print(similarity)
#' }
#'
#' @export
similarity_matrix_nmi <- function(partitions) {
    
    # Convert to matrix
    partition_matrix <- as.matrix(partitions)
    
    # Handle single partition (vector) case
    if (is.null(ncol(partition_matrix))) {
        partition_matrix <- matrix(partition_matrix, ncol = 1)
    }
    
    n_solutions <- ncol(partition_matrix)
    
    # Assign column names if missing
    if (is.null(colnames(partition_matrix))) {
        colnames(partition_matrix) <- paste0("sol", seq_len(n_solutions))
    }
    
    # Initialize k × k similarity matrix with diagonal = 1
    similarity_matrix <- diag(1, n_solutions)
    dimnames(similarity_matrix) <- list(
        colnames(partition_matrix), 
        colnames(partition_matrix)
    )
    
    # Compute pairwise similarities (upper triangle only, then mirror)
    if (n_solutions > 1) {
        for (i in 1:(n_solutions - 1)) {
            for (j in (i + 1):n_solutions) {
                
                # Compute NMI between partitions i and j
                nmi_value <- aricode::NMI(
                    partition_matrix[, i], 
                    partition_matrix[, j], 
                    variant = "sqrt"
                )
                
                # Fill upper and lower triangles (symmetric matrix)
                similarity_matrix[i, j] <- nmi_value
                similarity_matrix[j, i] <- nmi_value
            }
        }
    }
    
    return(similarity_matrix)
}