#' Generate Ring of Cliques Benchmark Networks
#'
#' Creates synthetic benchmark networks consisting of cliques arranged in a ring 
#' topology, with optional bridge nodes and/or a central hub. These networks are 
#' useful for testing community detection algorithms as they have known ground 
#' truth community structure.
#'
#' @param n_cliques Integer. Number of cliques in the ring (must be ≥ 2).
#' @param clique_size Integer. Number of nodes in each clique (must be ≥ 3).
#' @param variant Character. Network variant to generate:
#'   \itemize{
#'     \item{\code{"RC"}}{Basic ring: cliques connected in a ring}
#'     \item{\code{"RC+B"}}{Ring + bridge nodes between cliques}
#'     \item{\code{"RC+C"}}{Ring + central hub connected to all cliques}
#'     \item{\code{"RC+BC"}}{Ring + bridges + central hub}
#'   }
#'   Legacy names (\code{"RC_B"}, \code{"RC_C"}, \code{"RC_BC"}) are also accepted.
#' @param set_names Logical. If TRUE, assigns readable names to vertices 
#'   (default: TRUE).
#'
#' @return An igraph object with vertex and edge attributes:
#' \describe{
#'   \item{Vertex attributes:}{
#'     \itemize{
#'       \item{\code{clique_id}}{Clique membership (NA for bridges/center)}
#'       \item{\code{within_id}}{Position within clique (1 to clique_size)}
#'       \item{\code{gt_community}}{Ground truth community label}
#'       \item{\code{gt_label}}{Human-readable community label}
#'       \item{\code{role}}{Node type: "clique", "bridge", or "central"}
#'       \item{\code{name}}{Vertex name (if set_names = TRUE)}
#'     }
#'   }
#'   \item{Edge attributes:}{
#'     \itemize{
#'       \item{\code{edge_type}}{Type: "intra_clique", "ring_edge", "bridge_edge", 
#'         or "center_spoke"}
#'     }
#'   }
#' }
#'
#' @details
#' Network variants:
#' \itemize{
#'   \item{\strong{RC}}: Basic ring where clique k node 1 connects to clique k+1 node 2
#'   \item{\strong{RC+B}}: Adds bridge nodes between cliques: k→bridge_k→(k+1)
#'   \item{\strong{RC+C}}: Adds central hub connected to node 3 of each clique
#'   \item{\strong{RC+BC}}: Combines both bridges and central hub
#' }
#'
#' Ground truth communities: Each clique forms one community. Bridge nodes and 
#' the central node have gt_community = 0 (outliers).
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' 
#' # Basic ring of 4 cliques, each with 5 nodes
#' graph_basic <- make_ring_of_cliques(
#'   n_cliques = 4, 
#'   clique_size = 5, 
#'   variant = "RC"
#' )
#' 
#' # Ring with central hub (creates outlier node)
#' graph_with_center <- make_ring_of_cliques(
#'   n_cliques = 4, 
#'   clique_size = 5, 
#'   variant = "RC+C"
#' )
#' 
#' # Visualize ground truth communities
#' plot(graph_with_center, 
#'      vertex.color = V(graph_with_center)$gt_community,
#'      vertex.label = NA)
#' }
#'
#' @export
make_ring_of_cliques <- function(
        n_cliques,
        clique_size,
        variant = c("RC", "RC+C", "RC+B", "RC+BC", "RC_C", "RC_B", "RC_BC"),
        set_names = TRUE
) {
    
    # ---- Input validation ----
    stopifnot(
        is.numeric(n_cliques), 
        length(n_cliques) == 1L, 
        n_cliques >= 2
    )
    stopifnot(
        is.numeric(clique_size), 
        length(clique_size) == 1L, 
        clique_size >= 3
    )
    
    # ---- Normalize variant naming (support legacy names) ----
    variant <- match.arg(variant)
    variant <- switch(
        variant,
        "RC_C"  = "RC+C",
        "RC_B"  = "RC+B",
        "RC_BC" = "RC+BC",
        variant  # Keep as-is for new naming
    )
    
    # ---- Build disjoint cliques ----
    clique_list <- vector("list", n_cliques)
    for (clique_idx in seq_len(n_cliques)) {
        clique_list[[clique_idx]] <- igraph::make_full_graph(
            clique_size, 
            directed = FALSE
        )
    }
    
    # Combine cliques into single graph via disjoint union
    graph <- clique_list[[1]]
    if (n_cliques > 1) {
        for (clique_idx in 2:n_cliques) {
            graph <- igraph::disjoint_union(graph, clique_list[[clique_idx]])
        }
    }
    
    # ---- Assign vertex attributes ----
    
    # Clique membership for each node
    clique_membership <- rep(seq_len(n_cliques), each = clique_size)
    
    # Position within clique (1 to clique_size)
    position_within_clique <- rep(seq_len(clique_size), times = n_cliques)
    
    igraph::V(graph)$clique_id <- clique_membership
    igraph::V(graph)$within_id <- position_within_clique
    igraph::V(graph)$gt_community <- clique_membership  # Ground truth
    igraph::V(graph)$gt_label <- paste0("C", clique_membership)
    igraph::V(graph)$role <- "clique"
    
    # ---- Assign edge attributes ----
    # All existing edges are intra-clique
    igraph::E(graph)$edge_type <- "intra_clique"
    
    # ---- Helper functions for vertex indexing ----
    
    # Get vertex ID for node j in clique k
    get_vertex_id <- function(clique_k, node_j) {
        (clique_k - 1L) * clique_size + node_j
    }
    
    # Get next clique ID (wraps around for ring topology)
    get_next_clique <- function(clique_k) {
        if (clique_k < n_cliques) clique_k + 1L else 1L
    }
    
    # ---- Add ring edges (for RC and RC+C variants) ----
    if (variant %in% c("RC", "RC+C")) {
        
        ring_edge_list <- integer(0L)
        
        # Connect node 1 of clique k to node 2 of clique k+1
        for (clique_k in seq_len(n_cliques)) {
            next_clique <- get_next_clique(clique_k)
            ring_edge_list <- c(
                ring_edge_list,
                get_vertex_id(clique_k, 1),
                get_vertex_id(next_clique, 2)
            )
        }
        
        # Add ring edges to graph
        n_edges_before <- igraph::ecount(graph)
        graph <- igraph::add_edges(graph, ring_edge_list)
        
        # Mark new edges as ring edges
        new_edge_indices <- seq(n_edges_before + 1L, igraph::ecount(graph))
        igraph::E(graph)$edge_type[new_edge_indices] <- "ring_edge"
    }
    
    # ---- Add bridge nodes (for RC+B and RC+BC variants) ----
    bridge_vertex_ids <- integer(0L)
    
    if (variant %in% c("RC+B", "RC+BC")) {
        
        # Add n_cliques bridge nodes
        graph <- igraph::add_vertices(graph, n_cliques, role = "bridge")
        
        # Get IDs of newly added bridge nodes
        bridge_vertex_ids <- (igraph::vcount(graph) - n_cliques + 1L):igraph::vcount(graph)
        
        # Set bridge node attributes
        igraph::V(graph)$clique_id[bridge_vertex_ids] <- NA_integer_
        igraph::V(graph)$within_id[bridge_vertex_ids] <- NA_integer_
        igraph::V(graph)$gt_community[bridge_vertex_ids] <- 0L  # Outliers
        igraph::V(graph)$gt_label[bridge_vertex_ids] <- paste0("B", seq_len(n_cliques))
        
        # Build bridge edge list: clique k → bridge k → clique k+1
        bridge_edge_list <- integer(0L)
        
        for (clique_k in seq_len(n_cliques)) {
            bridge_id <- bridge_vertex_ids[clique_k]
            next_clique <- get_next_clique(clique_k)
            
            bridge_edge_list <- c(
                bridge_edge_list,
                get_vertex_id(clique_k, 1), bridge_id,  # k → bridge
                bridge_id, get_vertex_id(next_clique, 2)  # bridge → k+1
            )
        }
        
        # Add bridge edges to graph
        n_edges_before <- igraph::ecount(graph)
        graph <- igraph::add_edges(graph, bridge_edge_list)
        
        # Mark new edges as bridge edges
        new_edge_indices <- seq(n_edges_before + 1L, igraph::ecount(graph))
        igraph::E(graph)$edge_type[new_edge_indices] <- "bridge_edge"
    }
    
    # ---- Add central hub (for RC+C and RC+BC variants) ----
    center_vertex_id <- integer(0L)
    
    if (variant %in% c("RC+C", "RC+BC")) {
        
        # Add single central node
        graph <- igraph::add_vertices(graph, 1L, role = "central")
        center_vertex_id <- igraph::vcount(graph)
        
        # Set central node attributes
        igraph::V(graph)$clique_id[center_vertex_id] <- NA_integer_
        igraph::V(graph)$within_id[center_vertex_id] <- NA_integer_
        igraph::V(graph)$gt_community[center_vertex_id] <- 0L  # Outlier
        igraph::V(graph)$gt_label[center_vertex_id] <- "CENTRAL"
        
        # Build spoke edges: center → node 3 of each clique
        spoke_edge_list <- integer(0L)
        
        for (clique_k in seq_len(n_cliques)) {
            spoke_edge_list <- c(
                spoke_edge_list,
                center_vertex_id,
                get_vertex_id(clique_k, 3)
            )
        }
        
        # Add spoke edges to graph
        n_edges_before <- igraph::ecount(graph)
        graph <- igraph::add_edges(graph, spoke_edge_list)
        
        # Mark new edges as center spokes
        new_edge_indices <- seq(n_edges_before + 1L, igraph::ecount(graph))
        igraph::E(graph)$edge_type[new_edge_indices] <- "center_spoke"
    }
    
    # ---- Assign vertex names (optional) ----
    if (isTRUE(set_names)) {
        
        vertex_names <- character(igraph::vcount(graph))
        
        # Name clique nodes as "c{clique_id}_{position}"
        clique_node_indices <- seq_len(n_cliques * clique_size)
        vertex_names[clique_node_indices] <- paste0(
            "c", 
            igraph::V(graph)$clique_id[clique_node_indices], 
            "_", 
            igraph::V(graph)$within_id[clique_node_indices]
        )
        
        # Name bridge nodes as "bridge1", "bridge2", ...
        if (length(bridge_vertex_ids) > 0) {
            vertex_names[bridge_vertex_ids] <- paste0(
                "bridge", 
                seq_along(bridge_vertex_ids)
            )
        }
        
        # Name central node as "center"
        if (length(center_vertex_id) > 0) {
            vertex_names[center_vertex_id] <- "center"
        }
        
        igraph::V(graph)$name <- vertex_names
    }
    
    return(graph)
}