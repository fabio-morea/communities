#' Assign Community Labels Based on the Strongest Node in Each Community

#' Typically, community labels generated from community detection algorithms are arbitrary and unrelated to node identifiers.
#' In some cases, it can be useful to assign more meaningful labels to communities that are easier to interpret. A practical choice is to label the community after its "strongest" node, i.e., the node with the highest strength (the sum of the weights of its connections). 
#' This approach provides a label that reflects the most central or influential node within each community, making the labels more interpretable.
#'
#' With the `comm_label_as_strongest`, each community is labeled with the identifier of its strongest node.
#'
#' @param g An iGraph object representing the network. The graph should have weights for the edges if node strength is to be calculated based on weighted degrees.
#' @param comms A community detection object (such as from `cluster_*` functions in the iGraph package) that contains membership information, mapping each node to its community.
#'
#' @return A vector of community labels, where each node is assigned the label of the strongest node in its community. The labels take the form `"C_<node_name>"`, where `<node_name>` is the name of the strongest node in that community.
#'
#' @details The function works by calculating the node strength (the sum of weights of adjacent edges) for all nodes in the graph. For each community, the node with the highest strength is identified, and that node's name is used as the label for the entire community. This labeling helps to identify and label communities based on the most influential node within them.
#'
#' @examples
#' # Assuming 'g' is a graph and 'comms' is a community detection object
#' community_labels <- comm_label_as_strongest(g, comms)
#' print(community_labels)
#'
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