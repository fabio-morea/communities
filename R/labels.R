
#' @export
comm_label_as_strongest <- function(g, comms) {
    g1 <- g
    node_strength <- strength(g1)
    V(g1)$comm_labels <- "--"
    for (i in 1:max(membership(comms))) {
        community_nodes <- which(membership(comms) == i)
        strongest_node_within_community <- names(community_nodes)[which.max(node_strength[community_nodes])]
        
        print(strongest_node_within_community)
        community_label <- paste("C_", strongest_node_within_community, sep = "")
        V(g1)$comm_labels[community_nodes] <- community_label
    }
    return(V(g1)$comm_labels)
}