#' Calculate the Empirical Mixing Parameter (Mu)
#'
#' The `empirical_mu` function calculates the mixing parameter "mu" for a given network, which represents the proportion of edges that exist between different communities. This parameter is often used in community detection algorithms to quantify how mixed or modular a network is.
#'
#' @param g An iGraph object representing the network. The edges of the graph should ideally be weighted, but if not, all edge weights are set to 1.0 by default.
#' @param community_labels A vector of community labels, ordered according to the vertices in `g`. These labels define the community membership of each node in the graph.
#'
#' @return A numeric value representing the empirical mixing parameter (mu), which is the ratio of inter-community edges to the total number of edges.
#'
#' @details The function works by first ensuring that all edges have weights (defaulting to 1.0 if no weights are provided). It then calculates the number of edges that connect nodes from different communities (inter-community edges) and computes the proportion of such edges relative to the total number of edges in the graph.
#'
#' @examples
#' # Create a simple graph and calculate the mixing parameter
#' g <- make_ring_of_cliques(3, 5)
#' community_labels <- V(g)$community
#' mu <- empirical_mu(g, community_labels)
#' print(mu)
#'
#' @export
empirical_mu <- function(g, community_labels) {
    # Calculates the value of mixing parameter "mu"
    E(g)$weight <- 1.0
    V(g)$community <- community_labels
    gdf<-as_long_data_frame(g)
    gdf$inter_comm <- (gdf$from_community != gdf$to_community)
    inter_community_links <- sum(gdf$weight[ gdf$inter_comm == TRUE])
    mu = sum(inter_community_links) / sum(gdf$weight)
    return(mu)
}

#' Calculate the Number of Connected Components within Each Community
#'
#' The `internally_connected` function calculates how many connected components exist within each community of a given network. For each community, it extracts the subgraph containing only the nodes from that community and computes the number of connected components within that subgraph.
#'
#' @param g An iGraph object representing the network to be analyzed.
#' @param community_labels A vector of community labels, where each label corresponds to a node in the graph `g`. The length of this vector must match the number of vertices in the graph.
#'
#' @return A numeric vector where each element corresponds to the number of connected components in one of the communities. The order of the communities matches the order of the unique community labels in the input.
#'
#' @details The function first checks that the length of the `community_labels` vector matches the number of nodes in the graph. It then assigns the community labels to the vertices in the graph. For each community, the function extracts the subgraph corresponding to that community and counts the number of connected components in that subgraph. This can be useful for analyzing the internal structure of communities within a network.
#'
#' @examples
#' # Create a simple graph and assign community labels
#' g <- make_ring_of_cliques(4, 5)
#' community_labels <- V(g)$community
#' internal_components <- internally_connected(g, community_labels)
#' print(internal_components)
#' @export
internally_connected <- function(g, community_labels) {
    stopifnot(length(community_labels)==vcount(g))
    V(g)$community<-community_labels
    communities <- unique(V(g)$community)
    int_conn<- c()
    for (community in communities) {
        subgraph <- induced_subgraph(g, V(g)$community == community)
        int_conn <- c(int_conn, components(subgraph)$no) 
    }
    return(int_conn)  
}

#' Perform Quality Checks on a Set of Community Detection Solutions
#'
#' The `quality_check` function evaluates the quality of a set of community detection solutions 
#' stored in the `sol_space` object. The function computes various metrics, including 
#' the number of communities, modularity, mixing parameter (mu), and internal connectivity. 
#' It also flags solutions as valid or invalid based on certain thresholds.
#'
#' @param g An iGraph object representing the network to be analyzed.
#' @param sol_space A list-link object containing the solution space for community detection, 
#' produced by the 'solution_space()' function that is composed of two elements:
#'   - `M`: A matrix of community labels for each solution.
#'   - `data`: A dataframe where each row represents a solution, and columns will be added for 
#'     computed metrics like `k`, `mod`, `mu`, `int_conn`, and `valid`.
#'
#' @return The `sol_space` object with additional columns in `data`:
#'   - `k`: The number of unique communities in each solution.
#'   - `mod`: The modularity score for each solution.
#'   - `mu`: The empirical mixing parameter for each solution.
#'   - `int_conn`: The number of connected components within each community.
#'   - `valid`: A boolean value indicating whether the solution is valid (`TRUE`) or invalid (`FALSE`), based on a set of criteria.
#'
#' @details This function assesses the quality of each community detection solution in the solution space by calculating:
#'   - The number of unique communities (`k`).
#'   - The modularity score (`mod`), which measures the strength of community structure.
#'   - The mixing parameter (`mu`), indicating the proportion of inter-community edges.
#'   - A flag to highlight Whether communities are internally connected (`int_conn`).
#'   - A flag to signal the validity of each solution. A solution is not valid if any of the following condition applies: 
#'     - `mu > 0.5`: If the mixing parameter exceeds 0.5, the solution is marked as invalid.
#'     - `k == 1`: If there is only one community, the solution is marked as invalid.
#'     - `int_conn > 1`: If any community has more than one connected component, the solution is marked as invalid.
#'
#' @examples
#' # Assuming 'g' is a graph and 'sol_space' contains community detection results
#' sol_space_checked <- quality_check(g, sol_space)
#' print(sol_space_checked$data)
#'
#' @export
quality_check <- function(g, ssp, mu_max = 0.5) {
    stopifnot(inherits(g, "igraph"))
    if (is.null(ssp$partitions)) stop("ssp$partitions is missing.")
    
    # Partitions: accept either a vector (single solution) or a matrix (multiple solutions)
    M <- ssp$partitions
    if (is.null(ncol(M))) M <- matrix(M, ncol = 1)
    
    n_solutions <- ncol(M)
    n_nodes     <- igraph::vcount(g)
    
    # Preallocate result vectors
    k         <- numeric(n_solutions)
    mod       <- numeric(n_solutions)
    mu        <- numeric(n_solutions)
    int_conn  <- logical(n_solutions)
    valid     <- rep(TRUE, n_solutions)
    reasons   <- character(n_solutions)
    
    # Helper to append reasons for invalidity
    add_reason <- function(idx, txt) {
        if (is.na(reasons[idx]) || reasons[idx] == "") reasons[idx] <<- txt
        else reasons[idx] <<- paste(reasons[idx], txt, sep = " | ")
    }
    
    for (i in seq_len(n_solutions)) {
        labels_raw <- M[, i]
        
        # Normalize community labels to 1..K (avoids issues with 0 or sparse labels)
        labels <- as.integer(factor(labels_raw, levels = unique(labels_raw)))
        
        # 1) Number of communities
        k[i] <- dplyr::n_distinct(labels)
        
        # 2) Modularity
        mod[i] <- igraph::modularity(g, labels)
        
        # 3) Mixing parameter mu
        mu[i] <- communities::empirical_mu(g, labels)
        
        # 4) Internal connectivity: TRUE if all communities are internally connected
        ic_vec <- communities::internally_connected(g, labels) # usually logical per community
        int_conn[i] <- all(as.logical(ic_vec), na.rm = TRUE)
        
        # 5) Validity rules (all NA-safe)
        if (!is.na(mu[i]) && mu[i] > mu_max) {
            valid[i] <- FALSE; add_reason(i, sprintf("mu=%.3f > %.3f", mu[i], mu_max))
        }
        if (!is.na(k[i]) && k[i] == 1) {
            valid[i] <- FALSE; add_reason(i, "k=1 (trivial partition)")
        }
        if (!is.na(int_conn[i]) && !int_conn[i]) {
            valid[i] <- FALSE; add_reason(i, "communities not internally connected")
        }
        
        # If any metric is NA, mark solution as invalid and explain
        if (any(is.na(c(k[i], mod[i], mu[i], int_conn[i])))) {
            valid[i] <- FALSE
            add_reason(i, "NA in one or more metrics")
        }
    }
    
    tibble::tibble(
        solution_id = seq_len(n_solutions),
        k           = k,
        modularity  = mod,
        mu          = mu,
        int_conn    = int_conn,
        valid       = valid,
        reason      = reasons
    )
}


#' Similarity Matrix Between Community Partitions (NMI)
#'
#' Computes a pairwise similarity matrix between community-detection solutions
#' using Normalized Mutual Information (NMI, square-root normalization).
#' The measure is invariant to community label permutations.
#'
#' @param partitions A tibble, data.frame, or matrix with \eqn{n} rows (nodes)
#'   and \eqn{k} columns (solutions). Each column contains community labels for
#'   one partitioning of the same graph.
#'
#' @return A \eqn{k \times k} symmetric numeric matrix with values in \[0, 1\]:
#'   \itemize{
#'     \item \code{1} = identical partitions,
#'     \item \code{0} = completely dissimilar partitions.
#'   }
#'   Row and column names correspond to the input column names (or
#'   \code{"sol1"}, \code{"sol2"}, … if missing).
#'
#' @details
#' This function calls \code{\link[aricode]{NMI}} from the \pkg{aricode} package
#' with \code{variant = "sqrt"}, which ensures values are in \[0, 1\].
#' The diagonal is always set to \code{1}.
 
#' @export
similarity_matrix_nmi <- function(partitions) {
    M <- as.matrix(partitions)
    if (is.null(ncol(M))) M <- matrix(M, ncol = 1)
    
    k <- ncol(M)
    if (is.null(colnames(M))) colnames(M) <- paste0("sol", seq_len(k))
    
    # initialize k x k matrix, diagonals = 1
    S <- diag(1, k)
    dimnames(S) <- list(colnames(M), colnames(M))
    
    if (k > 1) {
        for (i in 1:(k - 1)) {
            for (j in (i + 1):k) {
                sim <- aricode::NMI(M[, i], M[, j], variant = "sqrt")
                S[i, j] <- sim   # upper triangle
                S[j, i] <- sim   # lower triangle (mirror)
            }
        }
    }
    return(S)
}


     
