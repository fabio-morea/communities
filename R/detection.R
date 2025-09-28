#' Explore Community Detection Solution Space Using Dirichlet-Multinomial Framework
#'
#' Systematically explores the solution space of community detection algorithms by
#' running multiple trials under randomized conditions (vertex permutation) and
#' modeling solution probabilities using a Dirichlet-Multinomial framework with
#' Beta credible intervals for convergence assessment.
#'
#' @param graph An \code{igraph} object representing the network to analyze.
#'   Must have at least 2 vertices. Vertex names will be added if missing.
#' @param n_trials Integer. Maximum number of algorithm trials to execute (default: 200).
#' @param method Character. Community detection algorithm to use. One of:
#'   \itemize{
#'     \item{\code{"IM"}}{Infomap}
#'     \item{\code{"WT"}}{Walktrap}
#'     \item{\code{"LV"}}{Louvain}
#'     \item{\code{"LD"}}{Leiden}
#'     \item{\code{"LP"}}{Label Propagation}
#'     \item{\code{"EV"}}{Leading Eigenvector}
#'     \item{\code{"EB"}}{Edge Betweenness}
#'   }
#'   Default is \code{"IM"}.
#' @param precision_threshold Numeric. Maximum allowed width of Beta credible 
#'   intervals for early stopping (default: 0.05). Smaller values require more precision.
#' @param confidence_level Numeric. Confidence level for credible intervals (default: 0.95).
#'   Must be between 0 and 1.
#' @param prior_concentration Numeric. Symmetric Dirichlet prior concentration 
#'   parameter (default: 1). Higher values express stronger prior beliefs.
#' @param random_seed Optional integer. Random seed for reproducible results.
#' @param verbose Logical. Whether to print progress information (default: TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{partitions}}{Matrix (n_vertices × n_solutions) of unique community 
#'     assignments found. Each column represents one solution.}
#'   \item{\code{probabilities}}{Data frame with posterior statistics for each solution:
#'     \code{id} (solution identifier), \code{count} (frequency), \code{phat} (posterior mean),
#'     \code{pmin}, \code{pmax} (credible interval bounds).}
#'   \item{\code{log}}{List of convergence diagnostics including \code{stop_trial},
#'     \code{stop_reason}, trial-by-trial counts, and probability evolution matrices.}
#' }
#'
#' @details
#' The algorithm uses vertex permutation to explore different algorithm trajectories,
#' ensuring solutions are mapped back to the original vertex ordering using vertex names.
#' Posterior inference assumes a symmetric Dirichlet prior over solution probabilities,
#' with early stopping when all marginal credible interval widths fall below the
#' precision threshold.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' 
#' # Simple example with ring network
#' ring_graph <- make_ring(20)
#' result <- solutions_space_DM(ring_graph, n_trials = 50, method = "IM")
#' print(result$probabilities)
#' 
#' # More complex example with multiple communities
#' complex_graph <- make_ring_of_cliques(n_cliques = 4, clique_size = 5)
#' result <- solutions_space_DM(
#'   graph = complex_graph,
#'   n_trials = 100,
#'   precision_threshold = 0.02,
#'   random_seed = 42
#' )
#' }
#'
#'
#' @export
solutions_space_DM <- function(
        graph,
        n_trials = 200,
        method = "IM",
        precision_threshold = 0.05,
        confidence_level = 0.95,
        prior_concentration = 1,
        random_seed = NULL,
        verbose = TRUE
) {
    
    # ---- Input validation ----
    .validate_graph_input(graph)
    .validate_algorithm_parameters(n_trials, method, precision_threshold, 
                                   confidence_level, prior_concentration)
    
    # ---- Initialize algorithm state ----
    n_vertices <- igraph::vcount(graph)
    if (is.null(igraph::V(graph)$name)) {
        igraph::V(graph)$name <- as.character(seq_len(n_vertices))
    }
    if (!is.null(random_seed)) set.seed(random_seed)
    
    # Algorithm-specific parameters
    algorithm_params <- list(
        shuffle_vertices = TRUE,
        infomap_trials = 10,
        walktrap_steps = 3,
        leiden_resolution = 1.0
    )
    
    # Internal helper functions
    partitions_identical <- function(x, y) {
        isTRUE(all.equal(aricode::ARI(x, y), 1))
    }
    
    compute_beta_credible_intervals <- function(counts) {
        n_solutions <- length(counts)
        alpha_level <- (1 - confidence_level) / 2
        
        if (n_solutions == 1L) {
            # Single solution case: Beta distribution
            alpha_param <- counts + prior_concentration
            beta_param <- prior_concentration
            
            list(
                posterior_means = alpha_param / (alpha_param + beta_param),
                lower_bounds = stats::qbeta(alpha_level, alpha_param, beta_param),
                upper_bounds = stats::qbeta(1 - alpha_level, alpha_param, beta_param)
            )
        } else {
            # Multiple solutions: marginal Beta distributions from Dirichlet
            alpha_params <- counts + prior_concentration
            total_alpha <- sum(alpha_params)
            
            list(
                posterior_means = alpha_params / total_alpha,
                lower_bounds = stats::qbeta(alpha_level, alpha_params, total_alpha - alpha_params),
                upper_bounds = stats::qbeta(1 - alpha_level, alpha_params, total_alpha - alpha_params)
            )
        }
    }
    
    # ---- Initialize solution space tracking ----
    partition_matrix <- matrix(numeric(0), nrow = n_vertices, ncol = 0)
    solution_counts <- numeric(0)
    names(solution_counts) <- character(0)
    rownames(partition_matrix) <- igraph::V(graph)$name
    
    # Trial-by-trial tracking
    counts_by_trial <- vector("list", n_trials)
    max_ci_width <- numeric(n_trials)
    stop_trial <- NA_integer_
    stop_reason <- NA_character_
    
    # ---- Main exploration loop ----
    for (trial in seq_len(n_trials)) {
        
        # Generate permuted graph for algorithmic diversity
        if (algorithm_params$shuffle_vertices) {
            vertex_permutation <- sample(n_vertices)
            permuted_graph <- igraph::permute(graph, vertex_permutation)
        } else {
            permuted_graph <- graph
        }
        
        # Apply community detection method
        membership_permuted <- .apply_detection_method(permuted_graph, method, algorithm_params)
        
        # Map results back to original vertex order using names
        original_names <- igraph::V(graph)$name
        permuted_names <- igraph::V(permuted_graph)$name
        membership_original <- membership_permuted[match(original_names, permuted_names)]
        
        # Check if this solution already exists
        solution_matched <- FALSE
        if (ncol(partition_matrix) > 0) {
            for (solution_idx in seq_len(ncol(partition_matrix))) {
                existing_membership <- partition_matrix[, solution_idx]
                if (partitions_identical(membership_original, existing_membership)) {
                    solution_counts[solution_idx] <- solution_counts[solution_idx] + 1
                    solution_matched <- TRUE
                    break
                }
            }
        }
        
        # Add new solution if not matched
        if (!solution_matched) {
            solution_id <- sprintf("s%02d", ncol(partition_matrix) + 1)
            partition_matrix <- cbind(partition_matrix, membership_original)
            colnames(partition_matrix)[ncol(partition_matrix)] <- solution_id
            solution_counts <- c(solution_counts, 1)
            names(solution_counts) <- colnames(partition_matrix)
        }
        
        # Store counts for this trial
        counts_by_trial[[trial]] <- solution_counts
        
        # Assess convergence using credible intervals
        credible_intervals <- compute_beta_credible_intervals(solution_counts)
        ci_widths <- credible_intervals$upper_bounds - credible_intervals$lower_bounds
        max_ci_width[trial] <- if (length(solution_counts) == 1L) {
            ci_widths
        } else {
            max(ci_widths)
        }
        
        # Progress reporting
        if (verbose) {
            message(sprintf("[Trial %d] Found %d unique solutions, max CI width = %.4f", 
                            trial, length(solution_counts), max_ci_width[trial]))
        }
        
        # Early stopping check
        if (max_ci_width[trial] <= precision_threshold) {
            stop_trial <- trial
            stop_reason <- "precision_achieved"
            break
        }
    }
    
    # Set final stop information if not already set
    if (is.na(stop_trial)) {
        stop_trial <- n_trials
        stop_reason <- "trial_budget_exhausted"
    }
    
    # ---- Finalize results ----
    
    # Compute final posterior statistics
    final_intervals <- compute_beta_credible_intervals(solution_counts)
    
    probabilities_df <- data.frame(
        id = names(solution_counts),
        count = as.integer(solution_counts),
        phat = final_intervals$posterior_means,
        pmin = final_intervals$lower_bounds,
        pmax = final_intervals$upper_bounds,
        row.names = NULL
    )
    
    # Sort by posterior probability (descending)
    sort_order <- order(probabilities_df$phat, decreasing = TRUE)
    probabilities_df <- probabilities_df[sort_order, , drop = FALSE]
    partition_matrix <- partition_matrix[, match(probabilities_df$id, 
                                                 colnames(partition_matrix)), 
                                         drop = FALSE]
    
    # Build trial-wise probability evolution matrices
    all_solution_ids <- colnames(partition_matrix)
    n_solutions <- length(all_solution_ids)
    
    # Only create matrices if we have solutions
    if (n_solutions > 0) {
        posterior_means_matrix <- matrix(NA_real_, nrow = stop_trial, ncol = n_solutions,
                                         dimnames = list(trial = seq_len(stop_trial), 
                                                         solution = all_solution_ids))
        lower_bounds_matrix <- posterior_means_matrix
        upper_bounds_matrix <- posterior_means_matrix
        
        for (trial in seq_len(stop_trial)) {
            trial_counts <- counts_by_trial[[trial]]
            if (length(trial_counts) > 0) {
                # Compute credible intervals for solutions known at this trial
                trial_intervals <- compute_beta_credible_intervals(trial_counts)
                
                # Align with solution names
                posterior_means <- as.numeric(trial_intervals$posterior_means)
                names(posterior_means) <- names(trial_counts)
                lower_bounds <- as.numeric(trial_intervals$lower_bounds)
                names(lower_bounds) <- names(trial_counts)
                upper_bounds <- as.numeric(trial_intervals$upper_bounds)
                names(upper_bounds) <- names(trial_counts)
                
                # Fill matrices for existing solutions at this trial
                existing_solutions <- intersect(names(trial_counts), all_solution_ids)
                posterior_means_matrix[trial, existing_solutions] <- posterior_means[existing_solutions]
                lower_bounds_matrix[trial, existing_solutions] <- lower_bounds[existing_solutions]
                upper_bounds_matrix[trial, existing_solutions] <- upper_bounds[existing_solutions]
            }
        }
        
        # Build counts matrix for reference
        counts_matrix <- matrix(0L, nrow = stop_trial, ncol = n_solutions,
                                dimnames = list(trial = seq_len(stop_trial), 
                                                solution = all_solution_ids))
        for (trial in seq_len(stop_trial)) {
            trial_counts <- counts_by_trial[[trial]]
            if (length(trial_counts) > 0) {
                counts_matrix[trial, names(trial_counts)] <- trial_counts
            }
        }
        
        # Create long-format probability evolution log (for plot_sol_space_evolution)
        probability_evolution_log <- do.call(
            rbind,
            lapply(seq_len(stop_trial), function(trial) {
                data.frame(
                    trial = trial,
                    solution = all_solution_ids,
                    phat = posterior_means_matrix[trial, all_solution_ids],
                    pmin = lower_bounds_matrix[trial, all_solution_ids],
                    pmax = upper_bounds_matrix[trial, all_solution_ids],
                    row.names = NULL
                )
            })
        )
        
    } else {
        # No solutions found - create empty structures
        posterior_means_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
        lower_bounds_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
        upper_bounds_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
        counts_matrix <- matrix(integer(0), nrow = 0, ncol = 0)
        probability_evolution_log <- data.frame(
            trial = integer(0),
            solution = character(0),
            phat = numeric(0),
            pmin = numeric(0),
            pmax = numeric(0)
        )
    }
    
    # ---- Return structured results ----
    list(
        partitions = partition_matrix,
        probabilities = probabilities_df,
        log = list(
            stop_trial = stop_trial,
            stop_reason = stop_reason,
            precision_threshold = precision_threshold,
            confidence_level = confidence_level,
            prior_concentration = prior_concentration,
            max_ci_width = head(max_ci_width, stop_trial),
            counts_matrix = counts_matrix,
            posterior_means_matrix = posterior_means_matrix,
            lower_bounds_matrix = lower_bounds_matrix,
            upper_bounds_matrix = upper_bounds_matrix,
            probability_evolution_log = probability_evolution_log,
            prob_long = probability_evolution_log  # Alias for backward compatibility
        )
    )
}

# ---- Helper Functions ----

#' Validate graph input parameters
#' @noRd
.validate_graph_input <- function(graph) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required but not available", call. = FALSE)
    }
    if (!requireNamespace("aricode", quietly = TRUE)) {
        stop("Package 'aricode' is required but not available", call. = FALSE)
    }
    
    if (!inherits(graph, "igraph")) {
        stop("`graph` must be an igraph object", call. = FALSE)
    }
    
    n_vertices <- igraph::vcount(graph)
    if (n_vertices < 2L) {
        stop("Graph must have at least 2 vertices", call. = FALSE)
    }
}

#' Validate algorithm parameters
#' @noRd
.validate_algorithm_parameters <- function(n_trials, method, precision_threshold, 
                                           confidence_level, prior_concentration) {
    if (!is.numeric(n_trials) || length(n_trials) != 1L || n_trials < 1L) {
        stop("`n_trials` must be a positive integer", call. = FALSE)
    }
    
    valid_methods <- c("IM", "WT", "LV", "LD", "LP", "EV", "EB")
    if (!method %in% valid_methods) {
        stop(sprintf("`method` must be one of: %s", 
                     paste(valid_methods, collapse = ", ")), call. = FALSE)
    }
    
    if (!is.numeric(precision_threshold) || precision_threshold <= 0) {
        stop("`precision_threshold` must be a positive number", call. = FALSE)
    }
    
    if (!is.numeric(confidence_level) || confidence_level <= 0 || confidence_level >= 1) {
        stop("`confidence_level` must be between 0 and 1", call. = FALSE)
    }
    
    if (!is.numeric(prior_concentration) || prior_concentration <= 0) {
        stop("`prior_concentration` must be positive", call. = FALSE)
    }
}

#' Apply the specified community detection method
#' @noRd
.apply_detection_method <- function(graph, method, algorithm_params) {
    switch(method,
           "IM" = igraph::infomap.community(graph, nb.trials = algorithm_params$infomap_trials)$membership,
           "WT" = igraph::walktrap.community(graph, steps = algorithm_params$walktrap_steps)$membership,
           "LD" = igraph::cluster_leiden(graph, resolution_parameter = algorithm_params$leiden_resolution)$membership,
           "LV" = igraph::cluster_louvain(graph)$membership,
           "LP" = igraph::label.propagation.community(graph)$membership,
           "EV" = igraph::cluster_leading_eigen(graph)$membership,
           "EB" = igraph::cluster_edge_betweenness(graph)$membership,
           stop("Unsupported method: ", method, call. = FALSE)
    )
}


#' Classify Solution Space Taxonomy
#'
#' Given the result of \code{solutions_space_DM()}, this function assigns a
#' high-level label to the solution space, based on the number and posterior
#' credibility of distinct partitions found.
#'
#' @param solution_space_result A list output from \code{solutions_space_DM()}. 
#'   Must contain at least \code{partitions}, \code{probabilities}, and \code{log}.
#'
#' @return A single character string indicating the type of solution space:
#' \describe{
#'   \item{\code{"Empty"}}{No valid partitions found (only degenerate ones).}
#'   \item{\code{"Single"}}{Only one solution is credible and clearly identified.}
#'   \item{\code{"Dominant"}}{One solution has a strictly dominant posterior.}
#'   \item{\code{"Multiple"}}{Multiple solutions with no clear dominance.}
#'   \item{\code{"Sparse"}}{Too many solutions relative to trials (unstable space).}
#'   \item{\code{"Uncertain"}}{Atypical configuration not falling into other types.}
#' }
#'
#' @details
#' The function checks the number of credible solutions, posterior intervals, and 
#' convergence log data to categorize the structure of the solution space. Degenerate 
#' solutions (e.g., all nodes in one group or all in separate groups) are excluded 
#' from classification.
#'
#' Classification logic:
#' \itemize{
#'   \item{\strong{Empty}}: No solutions or all solutions are degenerate
#'   \item{\strong{Single}}: Exactly one valid solution
#'   \item{\strong{Dominant}}: One solution's lower CI bound exceeds all others' upper bounds
#'   \item{\strong{Sparse}}: Many solutions relative to trials (ratio ≥ 0.20) with budget exhausted
#'   \item{\strong{Multiple}}: Multiple solutions, no clear dominance, max lower bound < 0.5
#'   \item{\strong{Uncertain}}: All other cases
#' }
#'
#' @seealso \code{\link{solutions_space_DM}} for generating the input object.
#'
#' @examples
#' \dontrun{
#' graph <- igraph::make_ring(10)
#' solution_space_result <- solutions_space_DM(graph, n_trials = 100, method = "IM")
#' solution_space_type(solution_space_result)
#' }
#'
#' @export
solution_space_type <- function(solution_space_result) {
    
    # ---- Input validation ----
    if (!is.list(solution_space_result)) {
        stop("`solution_space_result` must be a list", call. = FALSE)
    }
    
    required_elements <- c("probabilities", "partitions", "log")
    missing_elements <- setdiff(required_elements, names(solution_space_result))
    if (length(missing_elements) > 0) {
        stop(sprintf("`solution_space_result` must contain elements: %s. Missing: %s",
                     paste(required_elements, collapse = ", "),
                     paste(missing_elements, collapse = ", ")), call. = FALSE)
    }
    
    # ---- Extract components ----
    probabilities <- solution_space_result$probabilities
    partition_matrix <- solution_space_result$partitions
    log_info <- solution_space_result$log
    
    # ---- Handle empty solution space ----
    if (is.null(partition_matrix) || ncol(partition_matrix) == 0L || nrow(probabilities) == 0L) {
        return("Empty")
    }
    
    n_solutions <- ncol(partition_matrix)
    n_vertices <- nrow(partition_matrix)
    
    # ---- Check for degenerate partitions ----
    # Degenerate = all nodes in one community OR all nodes in separate communities
    valid_solutions <- logical(n_solutions)
    
    for (solution_idx in seq_len(n_solutions)) {
        community_labels <- partition_matrix[, solution_idx]
        
        # Skip solutions with missing values
        if (anyNA(community_labels)) {
            valid_solutions[solution_idx] <- FALSE
            next
        }
        
        n_communities <- length(unique(community_labels))
        
        # Valid if not degenerate (not all same, not all different)
        valid_solutions[solution_idx] <- (n_communities > 1L && n_communities < n_vertices)
    }
    
    # If no valid solutions, return Empty
    if (!any(valid_solutions)) {
        return("Empty")
    }
    
    # Filter to only valid solutions for analysis
    n_valid_solutions <- sum(valid_solutions)
    
    # ---- Handle single valid solution ----
    if (n_valid_solutions == 1L) {
        return("Single")
    }
    
    # ---- Validate probability data ----
    if (any(is.na(probabilities$phat)) || 
        any(is.na(probabilities$pmin)) || 
        any(is.na(probabilities$pmax))) {
        return("Uncertain")
    }
    
    # ---- Extract probability statistics ----
    posterior_means <- probabilities$phat
    lower_bounds <- probabilities$pmin
    upper_bounds <- probabilities$pmax
    
    if (length(posterior_means) == 0 || all(is.na(posterior_means))) {
        return("Empty")
    }
    
    # ---- Check for dominant solution ----
    # Dominant: top solution's lower bound > all other solutions' upper bounds
    if (length(posterior_means) >= 2L) {
        top_solution_idx <- which.max(posterior_means)
        other_indices <- setdiff(seq_along(posterior_means), top_solution_idx)
        
        top_lower_bound <- lower_bounds[top_solution_idx]
        max_other_upper_bound <- max(upper_bounds[other_indices], na.rm = TRUE)
        
        if (!is.na(top_lower_bound) && !is.na(max_other_upper_bound) && 
            top_lower_bound > max_other_upper_bound) {
            return("Dominant")
        }
    }
    
    # ---- Check for sparse solution space ----
    # Sparse: many solutions relative to trials AND budget exhausted
    if (!is.null(log_info$stop_reason) && !is.null(log_info$stop_trial)) {
        if (log_info$stop_reason == "trial_budget_exhausted") {
            solution_trial_ratio <- n_solutions / log_info$stop_trial
            if (solution_trial_ratio >= 0.20) {
                return("Sparse")
            }
        }
    }
    
    # ---- Check for multiple solutions ----
    # Multiple: no single solution dominates (max lower bound < 0.5)
    max_lower_bound <- max(lower_bounds, na.rm = TRUE)
    if (!is.na(max_lower_bound) && max_lower_bound < 0.5) {
        return("Multiple")
    }
    
    # ---- Default case ----
    return("Uncertain")
}


#' Build Pairwise Agreement Matrix (Gamma) from Solution Space
#'
#' Computes the pairwise agreement matrix Gamma where entry (u,v) represents
#' the probability that nodes u and v are assigned to the same community,
#' weighted by solution probabilities.
#'
#' @param partitions Matrix (n_vertices × n_solutions) of community assignments,
#'   OR a vector (length n_vertices) for single-partition case.
#' @param solution_probabilities Optional numeric vector of length n_solutions 
#'   with non-negative weights that sum to 1. Defaults to equal weights.
#' @param tolerance Numeric tolerance for probability sum validation (default: 1e-9).
#'
#' @return An n_vertices × n_vertices symmetric matrix with entries in [0, 1] 
#'   and unit diagonal.
#'
#' @details
#' The agreement matrix is computed as: Gamma[u,v] = sum_i p[i] * I(u and v 
#' co-occur in partition i), where I is the indicator function and p[i] is 
#' the probability of partition i.
#'
#' @examples
#' # Multi-partition example
#' partition_matrix <- cbind(
#'   c(1, 1, 2, 2, 2),  # partition 1
#'   c(1, 2, 2, 2, 2)   # partition 2
#' )
#' gamma_matrix <- build_gamma_matrix(partition_matrix, solution_probabilities = c(0.6, 0.4))
#'
#' @export
build_gamma_matrix <- function(partitions, solution_probabilities = NULL, tolerance = 1e-9) {
    
    # Convert vector to matrix for unified handling
    if (!is.matrix(partitions)) {
        if (!is.atomic(partitions) && !is.factor(partitions)) {
            stop("`partitions` must be a matrix or an atomic/factor vector", call. = FALSE)
        }
        partitions <- matrix(partitions, ncol = 1L)
    }
    
    # Validate dimensions
    n_vertices <- nrow(partitions)
    n_solutions <- ncol(partitions)
    if (n_vertices < 1L || n_solutions < 1L) {
        stop("`partitions` must have at least one row (node) and one column (partition)", 
             call. = FALSE)
    }
    if (anyNA(partitions)) {
        stop("`partitions` contains NA values; please remove them", call. = FALSE)
    }
    
    # Handle solution probabilities
    if (is.null(solution_probabilities)) {
        if (n_solutions == 1L) {
            solution_probabilities <- 1.0
        } else {
            stop("`solution_probabilities` must be provided for multiple partitions", 
                 call. = FALSE)
        }
    }
    
    # Validate probabilities
    if (!is.numeric(solution_probabilities) || length(solution_probabilities) != n_solutions) {
        stop("`solution_probabilities` must be numeric vector of length ncol(partitions)", 
             call. = FALSE)
    }
    if (any(!is.finite(solution_probabilities)) || any(solution_probabilities < 0)) {
        stop("`solution_probabilities` must be finite and non-negative", call. = FALSE)
    }
    
    probability_sum <- sum(solution_probabilities)
    if (!is.finite(probability_sum) || probability_sum <= 0) {
        stop("Sum of `solution_probabilities` must be positive and finite", call. = FALSE)
    }
    if (abs(probability_sum - 1) > max(tolerance, 1e-12)) {
        warning("`solution_probabilities` do not sum to 1; normalizing", call. = FALSE)
        solution_probabilities <- solution_probabilities / probability_sum
    }
    
    # Initialize agreement matrix
    gamma_matrix <- matrix(0.0, nrow = n_vertices, ncol = n_vertices)
    
    # Compute pairwise co-occurrence probabilities
    for (solution_idx in seq_len(n_solutions)) {
        weight <- solution_probabilities[solution_idx]
        community_labels <- partitions[, solution_idx]
        
        for (vertex_u in seq_len(n_vertices)) {
            for (vertex_v in seq_len(n_vertices)) {
                if (community_labels[vertex_u] == community_labels[vertex_v]) {
                    gamma_matrix[vertex_u, vertex_v] <- gamma_matrix[vertex_u, vertex_v] + weight
                }
            }
        }
    }
    
    # Ensure numerical stability and invariants
    gamma_matrix[gamma_matrix < 0] <- 0
    gamma_matrix[gamma_matrix > 1] <- 1
    diag(gamma_matrix) <- 1  # Each node always co-occurs with itself
    
    return(gamma_matrix)
}

#' Calculate Co-occurrence Matrix from Solution Space
#'
#' Computes a weighted co-occurrence matrix indicating how often pairs of nodes
#' appear together in the same community across different solutions, weighted
#' by solution probabilities.
#'
#' @param solution_space_result A list from \code{solutions_space_DM()} containing
#'   \code{partitions} (matrix of community assignments) and \code{probabilities}
#'   (data frame with solution weights).
#'
#' @return A symmetric n_vertices × n_vertices matrix with entries in [0, 1] 
#'   and unit diagonal, showing co-occurrence probabilities.
#'
#' @details
#' The co-occurrence matrix D[u,v] represents the probability that nodes u and v
#' are assigned to the same community, weighted by the posterior probabilities
#' of each solution. The diagonal is always 1 (each node co-occurs with itself).
#'
#' @examples
#' \dontrun{
#' graph <- make_ring_of_cliques(4, 5)
#' solution_space_result <- solutions_space_DM(graph, n_trials = 50)
#' co_occurrence_matrix <- co_occurrence(solution_space_result)
#' }
#'
#' @export
co_occurrence <- function(solution_space_result) {
    
    # Extract and validate partition matrix
    partition_matrix <- solution_space_result$partitions
    if (is.data.frame(partition_matrix)) partition_matrix <- as.matrix(partition_matrix)
    if (!is.numeric(partition_matrix)) {
        stop("solution_space_result$partitions must be a numeric matrix of community labels", 
             call. = FALSE)
    }
    
    n_vertices <- nrow(partition_matrix)
    n_solutions <- ncol(partition_matrix)
    if (is.null(n_vertices) || is.null(n_solutions) || n_vertices < 1 || n_solutions < 1) {
        stop("Partition matrix must have at least one row (node) and one column (solution)", 
             call. = FALSE)
    }
    
    # Extract solution weights with fallback to equal weights
    solution_weights <- solution_space_result$probabilities$phat
    if (is.null(solution_weights) || length(solution_weights) != n_solutions || 
        any(!is.finite(solution_weights))) {
        solution_weights <- rep(1 / n_solutions, n_solutions)
    } else {
        weight_sum <- sum(solution_weights)
        if (!isTRUE(all.equal(weight_sum, 1))) {
            solution_weights <- solution_weights / weight_sum
        }
    }
    
    # Initialize co-occurrence matrix with proper dimnames
    if (is.null(rownames(partition_matrix))) {
        rownames(partition_matrix) <- as.character(seq_len(n_vertices))
    }
    
    co_occurrence_matrix <- matrix(0, nrow = n_vertices, ncol = n_vertices,
                                   dimnames = list(rownames(partition_matrix), 
                                                   rownames(partition_matrix)))
    
    # Accumulate weighted co-occurrences
    for (solution_idx in seq_len(n_solutions)) {
        community_labels <- partition_matrix[, solution_idx]
        
        # Skip solutions with missing labels
        if (anyNA(community_labels)) next
        
        # Process each community in this solution
        unique_communities <- sort(unique(community_labels))
        for (community_id in unique_communities) {
            community_members <- which(community_labels == community_id)
            n_members <- length(community_members)
            
            # Update pairwise co-occurrences within this community
            if (n_members >= 2) {
                for (member_i in 1:(n_members - 1)) {
                    for (member_j in (member_i + 1):n_members) {
                        vertex_i <- community_members[member_i]
                        vertex_j <- community_members[member_j]
                        
                        # Add weighted co-occurrence (symmetric)
                        weight <- solution_weights[solution_idx]
                        co_occurrence_matrix[vertex_i, vertex_j] <- co_occurrence_matrix[vertex_i, vertex_j] + weight
                        co_occurrence_matrix[vertex_j, vertex_i] <- co_occurrence_matrix[vertex_i, vertex_j]
                    }
                }
            }
        }
    }
    
    # Ensure diagonal entries are 1 (each node co-occurs with itself)
    diag(co_occurrence_matrix) <- 1
    
    return(co_occurrence_matrix)
}

#' Identify Consensus Communities from Co-occurrence Matrix
#'
#' Identifies consensus communities from a co-occurrence matrix by grouping nodes
#' based on a threshold for pairwise co-occurrence. Computes uncertainty coefficients
#' reflecting confidence in community assignments.
#'
#' @param co_occurrence_matrix A symmetric co-occurrence matrix where entries 
#'   represent pairwise co-occurrence probabilities across multiple trials.
#' @param co_occurrence_threshold Numeric threshold for defining communities. 
#'   Nodes are considered in the same community if their pairwise co-occurrence 
#'   exceeds this value.
#' @param group_outliers Logical indicating whether single-node communities 
#'   (outliers) should be grouped together (default: FALSE).
#' @param verbose Logical indicating whether to print progress information 
#'   (default: FALSE).
#'
#' @return A data frame containing:
#' \describe{
#'   \item{\code{node_name}}{Name of each node}
#'   \item{\code{consensus_community_label}}{Consensus community assignment}
#'   \item{\code{uncertainty_coefficient}}{Confidence measure (gamma) for assignment}
#'   \item{\code{community_size}}{Size of assigned community}
#'   \item{\code{is_singleton}}{Whether node forms single-node community}
#' }
#'
#' @details
#' The algorithm iteratively identifies communities as connected components in
#' the thresholded co-occurrence matrix. The uncertainty coefficient gamma is
#' computed as 1 - mean(co-occurrence values) for nodes that co-occur at least
#' once with the focal node.
#'
#' @examples
#' \dontrun{
#' # Assume co_occurrence_matrix is a co-occurrence matrix
#' consensus_result <- consensus_communities(
#'   co_occurrence_matrix, 
#'   co_occurrence_threshold = 0.5, 
#'   group_outliers = TRUE
#' )
#' print(consensus_result)
#' }
#'
#' @export
consensus_communities <- function(
        co_occurrence_matrix, 
        co_occurrence_threshold, 
        group_outliers = FALSE, 
        verbose = FALSE
) {
    
    # Validate input matrix
    if (!is.matrix(co_occurrence_matrix) || !isSymmetric(co_occurrence_matrix)) {
        stop("co_occurrence_matrix must be a symmetric matrix", call. = FALSE)
    }
    if (is.null(colnames(co_occurrence_matrix))) {
        colnames(co_occurrence_matrix) <- rownames(co_occurrence_matrix) <- 
            as.character(seq_len(nrow(co_occurrence_matrix)))
    }
    
    # Initialize results data frame
    node_names <- colnames(co_occurrence_matrix)
    n_nodes <- length(node_names)
    
    results <- data.frame(
        node_name = node_names,
        is_processed = FALSE,
        temporary_community_label = NA_integer_,
        uncertainty_coefficient = NA_real_,
        community_size = NA_integer_,
        is_singleton = FALSE,
        stringsAsFactors = FALSE
    )
    
    community_counter <- 0
    nodes_remaining <- n_nodes
    
    # Main community identification loop
    while (nodes_remaining > 0) {
        community_counter <- community_counter + 1
        
        # Find first unprocessed node
        first_unprocessed <- which.max(results$is_processed == FALSE)
        
        # Identify nodes connected to this node above threshold
        connected_nodes <- (co_occurrence_matrix[first_unprocessed, ] > co_occurrence_threshold)
        
        # Calculate uncertainty coefficients for nodes in this community
        community_co_occurrence <- co_occurrence_matrix[connected_nodes, , drop = FALSE]
        
        # Set zero co-occurrences to NA (nodes never in same community)
        community_co_occurrence[community_co_occurrence == 0] <- NA
        
        if (sum(connected_nodes) > 1) {
            # Multi-node community: compute row-wise means
            uncertainty_values <- 1 - apply(community_co_occurrence, 1, mean, na.rm = TRUE)
        } else {
            # Single-node community: compute overall mean
            uncertainty_values <- 1 - mean(community_co_occurrence, na.rm = TRUE)
        }
        
        # Update results for nodes in this community
        results$temporary_community_label[connected_nodes] <- community_counter
        results$community_size[connected_nodes] <- sum(connected_nodes)
        results$uncertainty_coefficient[connected_nodes] <- uncertainty_values
        results$is_processed[connected_nodes] <- TRUE
        
        # Update remaining nodes count
        nodes_remaining <- sum(!results$is_processed)
        
        if (verbose) {
            message(sprintf("Community %d: %d nodes, %d remaining", 
                            community_counter, sum(connected_nodes), nodes_remaining))
        }
    }
    
    # Handle missing uncertainty coefficients
    results$uncertainty_coefficient[is.na(results$uncertainty_coefficient)] <- 0.0
    
    # Identify singleton communities
    results$is_singleton <- (results$community_size == 1)
    
    # Optionally group outliers
    if (group_outliers) {
        results$temporary_community_label[results$is_singleton] <- 0
    }
    
    # Create final community labels sorted by size (largest first)
    community_size_summary <- results %>%
        dplyr::group_by(temporary_community_label) %>%
        dplyr::summarise(n_members = dplyr::n(), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(n_members)) %>%
        dplyr::mutate(consensus_community_label = dplyr::row_number())
    
    # Join with original results and clean up
    final_results <- results %>%
        dplyr::inner_join(community_size_summary, by = "temporary_community_label") %>%
        dplyr::select(
            node_name, 
            consensus_community_label, 
            uncertainty_coefficient, 
            community_size, 
            is_singleton
        ) %>%
        dplyr::arrange(consensus_community_label, node_name)
    
    return(final_results)
}