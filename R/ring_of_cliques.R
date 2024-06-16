
#' Functions to make a "ring of cliques"
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

 