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
#'     computed metrics like `k`, `k_n`, `mod`, `mu`, `int_conn`, and `valid`.
#'
#' @return The `sol_space` object with additional columns in `data`:
#'   - `k`: The number of unique communities in each solution.
#'   - `k_n`: The number of communities normalized by the total number of nodes.
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
#'     - `k_n == 1`: If the number of communities is equal to the total number of nodes, the solution is marked as invalid.
#'     - `int_conn > 1`: If any community has more than one connected component, the solution is marked as invalid.
#'
#' @examples
#' # Assuming 'g' is a graph and 'sol_space' contains community detection results
#' sol_space_checked <- quality_check(g, sol_space)
#' print(sol_space_checked$data)
#'
#' @export
quality_check <- function(g, sol_space){
    
    M <- sol_space$M 

    sol_space$data$valid <- TRUE
    
    n_solutions = nrow(sol_space$data)
    
    for (i in 1:n_solutions){
        if (n_solutions == 1){
            labels <- M    
        } else {
            labels <- M[,i]
        } 
        sol_space$data$k[i] <- length(labels %>% unique())
        sol_space$data$k_n[i] <- sol_space$data$k[i] / vcount(g)
        
        # modularity: need to use c_membership + 1 to handle community label 0
        sol_space$data$mod[i] <- igraph::modularity (g,  labels + 1)
        
        # mixing parameter mu
        sol_space$data$mu[i] <- communities::empirical_mu(g, labels)
        
        # communities are internally connected
        sol_space$data$int_conn[i]  <- max(communities::internally_connected(g, labels))
        
        # validity
        if (sol_space$data$mu[i] > 0.5){sol_space$data$valid[i] <- FALSE}
        if (sol_space$data$k[i] == 1){sol_space$data$valid[i] <- FALSE}
        if (sol_space$data$k_n[i] == 1){sol_space$data$valid[i] <- FALSE}
        if (sol_space$data$int_conn[i] > 1 ){sol_space$data$valid[i] <- FALSE}
    }

    return(sol_space)
}