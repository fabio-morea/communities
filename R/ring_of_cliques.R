#' Ring of Cliques (RC, RC_B, RC_C, RC_BC)
#' - RC:    ring edges from c(k,1) to c(k+1,2)
#' - RC_B:  RC + bridge nodes: c(k,2) -- b_k -- c(k+1,3)
#' - RC_C:  RC + center connected to c(k,1)
#' - RC_BC: RC + bridges + center
#' @param n_cliques integer >= 2
#' @param clique_size integer >= 3 (needs nodes 1,2,3)
#' @param variant one of "RC","RC_B","RC_C","RC_BC"
#' @param set_names logical, add readable names
#' @return igraph
#' 
#' 
#' @export
make_ring_of_cliques <- function(n_cliques,
                                 clique_size,
                                 variant = c("RC", "RC+C", "RC+B", "RC+BC",
                                             "RC_C", "RC_B", "RC_BC"),
                                 set_names = TRUE) {
    stopifnot(is.numeric(n_cliques), length(n_cliques) == 1L, n_cliques >= 2)
    stopifnot(is.numeric(clique_size), length(clique_size) == 1L, clique_size >= 3)
    
    # normalize variant to the new naming
    variant <- match.arg(variant)
    variant <- switch(variant,
                      "RC_C"  = "RC+C",
                      "RC_B"  = "RC+B",
                      "RC_BC" = "RC+BC",
                      variant
    )
    
    # Build disjoint cliques
    cliques <- vector("list", n_cliques)
    for (i in seq_len(n_cliques)) cliques[[i]] <- igraph::make_full_graph(clique_size, directed = FALSE)
    g <- cliques[[1]]
    if (n_cliques > 1) for (i in 2:n_cliques) g <- igraph::disjoint_union(g, cliques[[i]])
    
    # Vertex attributes
    clique_id <- rep(seq_len(n_cliques), each = clique_size)
    within_id <- rep(seq_len(clique_size), times = n_cliques)
    igraph::V(g)$clique_id <- clique_id
    igraph::V(g)$within_id <- within_id
    igraph::V(g)$gt_community <- clique_id
    igraph::V(g)$gt_label <- paste0("C", clique_id)
    igraph::V(g)$role <- "clique"
    
    # Existing edges are intra-clique
    igraph::E(g)$edge_type <- "intra_clique"
    
    # Helpers
    vid    <- function(k, j) (k - 1L) * clique_size + j
    next_k <- function(k) if (k < n_cliques) k + 1L else 1L
    
    # --- RC (direct ring): c(k,1) -> c(k+1,2) ---
    if (variant %in% c("RC", "RC+C")) {
        ring_edges <- integer(0L)
        for (k in seq_len(n_cliques)) {
            ring_edges <- c(ring_edges, vid(k, 1), vid(next_k(k), 2))
        }
        e0 <- igraph::ecount(g)
        g  <- igraph::add_edges(g, ring_edges)
        igraph::E(g)$edge_type[seq(e0 + 1L, igraph::ecount(g))] <- "ring_edge"
    }
    
    # --- RC+B: add one bridge per clique and wire c(k,1) -- b_k -- c(k+1,2) ---
    bridge_ids <- integer(0L)
    if (variant %in% c("RC+B", "RC+BC")) {
        g <- igraph::add_vertices(g, n_cliques, role = "bridge")
        bridge_ids <- (igraph::vcount(g) - n_cliques + 1L):igraph::vcount(g)
        igraph::V(g)$clique_id[bridge_ids] <- NA_integer_
        igraph::V(g)$within_id[bridge_ids] <- NA_integer_
        igraph::V(g)$gt_community[bridge_ids] <- 0L
        igraph::V(g)$gt_label[bridge_ids] <- paste0("B", seq_len(n_cliques))
        
        bridge_edges <- integer(0L)
        for (k in seq_len(n_cliques)) {
            b <- bridge_ids[k]
            bridge_edges <- c(bridge_edges, vid(k, 1), b,  b, vid(next_k(k), 2))
        }
        e1 <- igraph::ecount(g)
        g  <- igraph::add_edges(g, bridge_edges)
        igraph::E(g)$edge_type[seq(e1 + 1L, igraph::ecount(g))] <- "bridge_edge"
    }
    
    # --- +C: center with spokes to c(k,3) ---
    center_vid <- integer(0L)
    if (variant %in% c("RC+C", "RC+BC")) {
        g <- igraph::add_vertices(g, 1L, role = "central")
        center_vid <- igraph::vcount(g)
        igraph::V(g)$clique_id[center_vid] <- NA_integer_
        igraph::V(g)$within_id[center_vid] <- NA_integer_
        igraph::V(g)$gt_community[center_vid] <- 0L
        igraph::V(g)$gt_label[center_vid] <- "CENTRAL"
        
        spokes <- integer(0L)
        for (k in seq_len(n_cliques)) spokes <- c(spokes, center_vid, vid(k, 3))
        e2 <- igraph::ecount(g)
        g  <- igraph::add_edges(g, spokes)
        igraph::E(g)$edge_type[seq(e2 + 1L, igraph::ecount(g))] <- "center_spoke"
    }
    
    # Names
    if (isTRUE(set_names)) {
        v_names <- character(igraph::vcount(g))
        idx_clique <- seq_len(n_cliques * clique_size)
        v_names[idx_clique] <- paste0("c", igraph::V(g)$clique_id[idx_clique], "_", igraph::V(g)$within_id[idx_clique])
        if (length(bridge_ids)) v_names[bridge_ids] <- paste0("bridge", seq_along(bridge_ids))
        if (length(center_vid)) v_names[center_vid] <- "center"
        igraph::V(g)$name <- v_names
    }
    
    g
}

 
