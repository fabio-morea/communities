#' Assign Community Labels Based on the Strongest Node in Each Community
#'
#' Typically, community labels generated from community detection algorithms are 
#' arbitrary and unrelated to node identifiers. In some cases, it can be useful 
#' to assign more meaningful labels to communities that are easier to interpret. 
#' A practical choice is to label the community after its "strongest" node, i.e., 
#' the node with the highest strength (the sum of the weights of its connections). 
#' This approach provides a label that reflects the most central or influential 
#' node within each community, making the labels more interpretable.
#'
#' With `comm_label_as_strongest`, each community is labeled with the identifier 
#' of its strongest node.
#'
#' @param graph An igraph object representing the network. The graph should have 
#'   weights for the edges if node strength is to be calculated based on weighted 
#'   degrees.
#' @param community_membership A community detection object (such as from 
#'   \code{cluster_*} functions in the igraph package) that contains membership 
#'   information, mapping each node to its community.
#'
#' @return A vector of community labels, where each node is assigned the label 
#'   of the strongest node in its community. The labels take the form 
#'   \code{"C_<node_name>"}, where \code{<node_name>} is the name of the 
#'   strongest node in that community.
#'
#' @details 
#' The function works by calculating the node strength (the sum of weights of 
#' adjacent edges) for all nodes in the graph. For each community, the node with 
#' the highest strength is identified, and that node's name is used as the label 
#' for the entire community. This labeling helps to identify and label communities 
#' based on the most influential node within them.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' 
#' # Create a simple graph
#' graph <- make_graph("Zachary")
#' 
#' # Detect communities
#' community_result <- cluster_louvain(graph)
#' 
#' # Assign labels based on strongest nodes
#' community_labels <- comm_label_as_strongest(graph, community_result)
#' }
#'
#' @export
comm_label_as_strongest <- function(graph, community_membership) {
    
    # Create working copy of graph
    graph_working <- graph
    
    # Calculate node strength (sum of edge weights for each node)
    node_strength <- igraph::strength(graph_working)
    
    # Initialize community labels for all nodes
    igraph::V(graph_working)$comm_labels <- "--"
    
    # Get the number of communities
    n_communities <- max(igraph::membership(community_membership))
    
    # Process each community
    for (community_id in seq_len(n_communities)) {
        
        # Find all nodes belonging to this community
        community_node_indices <- which(igraph::membership(community_membership) == community_id)
        
        # Get strength values for nodes in this community
        community_node_strengths <- node_strength[community_node_indices]
        
        # Find the strongest node within this community
        strongest_node_local_index <- which.max(community_node_strengths)
        strongest_node_name <- names(community_node_indices)[strongest_node_local_index]

        # Create community label based on strongest node
        community_label <- paste0("C_", strongest_node_name)
        
        # Assign this label to all nodes in the community
        igraph::V(graph_working)$comm_labels[community_node_indices] <- community_label
    }
    
    # Return the vector of community labels
    return(igraph::V(graph_working)$comm_labels)
}
