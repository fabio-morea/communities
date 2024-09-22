#' Create a Clique and Assign a Community Label
#'
#' The `make_clique` function generates a fully connected subgraph (i.e. a clique) of a specified size 
#' and assigns a community label to all nodes in the clique. This function is used mainly by `make_ring_of_cliques()`.
#'
#' @param clique_size An integer representing the number of nodes in the clique. All nodes will 
#' be fully connected to each other.
#' @param comm_label A numeric or character value used to assign a community label to the nodes 
#' within the clique. The label will be stored in the `community` node attribute.
#'
#' @return An undirected iGraph object representing a clique, where all nodes are connected. 
#' Each node will have a `community` attribute assigned according to the `comm_label` parameter.
#'
#' @details A clique is a subset of a graph where every node is connected to every other node. 
#' This function constructs such a graph using a given number of nodes and allows the assignment 
#' of a common community label to those nodes. The returned graph is undirected and can be combined 
#' with other cliques or networks.
#'
#' @examples
#' # Create a clique of size 5 with community label 1
#' clique <- make_clique(clique_size = 5, comm_label = 1)
#' 
#' @export
make_clique <- function(clique_size, comm_label) {
    require(igraph)
    G <- graph.empty(n = clique_size)
    edges <- t(combn(1:clique_size, 2))
    for (e in 1:nrow(edges)) {
        G <- add_edges(G, edges[e, ])
    }
    V(G)$community <- comm_label
    
    return(as.undirected(G))
}


#' Create a Ring of Cliques Graph with Optional Bridge and Central Node
#'
#' The `make_ring_of_cliques` function generates a graph composed of multiple fully connected 
#' subgraphs (cliques) arranged in a ring. Optionally, "bridge nodes" can be added between cliques, and 
#' a central node can be added that connects to all cliques.
#' 
#' @param num_cliques An integer specifying the number of cliques o include in the ring.
#' @param clique_size An integer specifying the number of nodes in each clique.
#' @param add_center A logical value indicating whether to add a central node that connects to all cliques. 
#' Default is `TRUE`.
#' @param add_bridges A logical value indicating whether to add bridge nodes between adjacent cliques in the ring. 
#' Default is `TRUE`.
#'
#' @return An iGraph object representing the ring of cliques. Each node is assigned a community label, 
#' with cliques labeled as `C1, C2, ...`, bridge nodes labeled as `B1, B2, ...`, and the central node labeled as `A`. The graph is undirected. All edges have weight = 1.
#'
#' @details This function creates a network structure where cliques are connected in a ring formation. 
#' Rings of cliques are useful for testing community detection algorithms: varying the parameters one can create several toy examples of controlled complexity, ranging from trivially simple to particularly challenging.  
#' - If `add_bridges` is `TRUE`, additional bridge nodes are inserted between each adjacent pair of cliques, 
#' with edges connecting the cliques.
#' - If `add_center` is `TRUE`, a central node is added, connected to a node in each clique, creating a star-like structure.
#'
#' @examples
#' # Create a ring of 4 cliques, each with 5 nodes, with both bridges and a central node
#' ring_of_cliques <- make_ring_of_cliques(num_cliques = 4, clique_size = 5, add_center = TRUE, add_bridges = TRUE)
#'
#' # Create a ring of 3 cliques without bridges or a central node
#' ring_of_cliques_no_center <- make_ring_of_cliques(num_cliques = 3, clique_size = 4, add_center = FALSE, add_bridges = FALSE)
#'
#' @export 
make_ring_of_cliques <- function(num_cliques,
                                 clique_size,
                                 add_center = TRUE,
                                 add_bridges = TRUE ) {
    require(igraph)
    gg <- as.undirected(graph.empty())
    
    for (i in 1:num_cliques) {
        next_clique <- make_clique(clique_size, comm_label = paste0("C", i))
        gg <- gg + next_clique
    }
    
    b <- vcount(gg)
    
    
    if (add_bridges) {
        gg <- add_vertices(gg, num_cliques)
    }
    
    
    for (j in 1:(num_cliques)) {
        b <- b + 1
        b_start <- (j-1) * clique_size +1
        b_end <- b_start + clique_size +1
        if (b_end > (clique_size * num_cliques)) {b_end <- 2}
        if (add_bridges) {
            gg <- add_edges(gg, c(b_start, b))
            gg <- add_edges(gg, c(b, b_end))
            V(gg)$community[b] <- paste0("B", j)
        } else {
            gg <- add_edges(gg, c(b_start, b_end))
        }
        
    }
    
    if (add_center) {
        gg <- add_vertices(gg, 1)
        id_center <- vcount(gg)
        V(gg)$community[id_center] <- "A"
        for (j in 1:(num_cliques)) {
            c_start <- (j-1) * clique_size +3
            gg <- add_edges(gg, c(c_start , id_center))
        }
        
    }
    
    E(gg)$weight <- 1.0
    
    V(gg)$name <- 1:vcount(gg)
    V(gg)$id <-   1:vcount(gg)
    
    
    return(gg)
    
}

 