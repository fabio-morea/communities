#' Build a Community Network from an iGraph Object
#'
#' This function constructs a new network (`Gc`), where each node represents a community
#' within the original network (`G`). The resulting community network can be used to 
#' analyze the relationships between communities and explore their structure and interactions.
#'
#' @param g An iGraph object representing the original network to be analyzed. The network 
#' must include the following attributes:
#'   - `V(g)$community`: An integer vector representing the community assignment for each node.
#'   - `E(g)$w`: A numeric vector representing edge weights. If the network is unweighted, 
#' set `E(g)$w` to 1.0 for all edges.
#'
#' @return An iGraph object representing the community network (`Gc`) with the following attributes:
#'   - `V(Gc)$membership`: Community labels corresponding to the ones in the original network.
#'   - `E(Gc)$w`: Edge weights representing the sum of edge weights between communities in the original network.
#'   - `V(Gc)$size`: The number of nodes from the original network `G` that belong to each community in `Gc`.
#'
#' @details This function aggregates the nodes and edges of the original network based on 
#' their community membership, creating a higher-level network where each node represents 
#' a community, and edges between nodes represent the aggregated connections between those communities.
#'
#' @examples
#' # Assuming 'g' is a pre-existing iGraph object with community and edge weight attributes
#' community_network <- make_community_network(g)
#'
#'   
#' @export
make_community_network <- function (g) {
    edges_list <- g %>%
        as_long_data_frame() %>%
        select(from_community, to_community, w) %>%
        group_by(from_community, to_community) %>%
        summarize(weight = sum(w))
    
    gc <- graph_from_data_frame(edges_list, directed = FALSE)
    V(gc)$id <- (1:vcount(gc))
    
    comms <- data.frame(label = V(gc)$name)
    comms$size <- 0
    for (i in 1:length(comms$label)) {
        comms$size[i] <-
            length(V(g)$community[V(g)$community == comms$label[i]])
    }
    V(gc)$size <- comms$size
    
    return(gc)
}



#' Calculate Distances for Community-Based Network Layout
#'
#' The `layout_distance_comm` function computes distances between pairs of nodes in a network 
#' to enhance visualization. The function clusters nodes belonging to the same community more 
#' closely while pushing nodes from different communities further apart.
#'
#' @param g An iGraph object representing the network to be analyzed. The network must have edge weights 
#' stored in `E(g)$weight`. If the network is unweighted, set `E(g)$weight` to 1.0 for all edges.
#' @param membership A numeric vector representing the community membership of each node, where each 
#' entry corresponds to a node in the graph `g` (ordered as `V(g)`).
#' @param eps A numeric value representing the minimum distance between nodes within the same community. 
#' The default value is 0.02, and typical values range from 0.1 to 0.2.
#'
#' @return A numeric vector of distances between pairs of nodes, ordered according to the edges in `g` (i.e., `E(g)`).
#' Nodes in the same community will have smaller distances, while nodes in different communities will be spread further apart.
#'
#' @details The function calculates distances based on edge weights and community membership. The parameter 
#' `eps` sets the minimum distance between nodes within the same community, while the actual distance depends 
#' on the edge weight between the nodes.
#'
#' @examples
#' # Assuming 'g' is a pre-existing iGraph object and 'membership' is a vector of community memberships
#' distances <- layout_distance_comm(g, membership, eps = 0.15)
#'
#' @export
layout_distance_comm <- function(g, membership, eps = .02) {
    if (!is.igraph(g)) {
        stop("g is not an igraph graph.")
    }
    if (any(is.na(E(g)$weight))) {
        stop("Edge attribute 'weight' contains NA values.")
    }
    
    df <- as_long_data_frame(g) %>% select(from, to, weight)
    df$dist <- 0
    for (i in 1:nrow(df)) {
        same_comm = (membership[df$from[i]] == membership[df$to[i]])
        df$dist[i] <- eps + (1 - eps) * df$weight[i] * same_comm
    }
    return(df$dist)
}


#' Plot the Solution Space of Community Detection Results
#'
#' The `plot_sol_space` function generates a series of plots to visually explore the solution space
#' derived from community detection analysis. It provides insights into the distribution of solutions, 
#' their internal consistency, and similarities between different partitions.
#'
#' @param sol_space A dataframe produced by the `explore_solution_space` function. The dataframe should contain
#' data on community partitions and relevant statistics such as cumulative sums, medians, and upper/lower bounds 
#' for each solution.
#'
#' @return A list of up to 4 plots:
#'   - `pl1`: A line plot showing the cumulative sum and median of solutions, with error bars indicating variability.
#'   - `pl2`: A visualization of the solution space, showing the range of valid and non-valid solutions.
#'   - `pl3`: A scatter plot showing the distribution of community sizes for each solution.
#'   - `pl4`: A plot representing the similarity between the solutions (empty if solution space contains only one solution).
#'
#' @examples
#' # Assuming 'sol_space' is a dataframe produced by 'explore_solution_space'
#' plot_list <- plot_sol_space(sol_space)
#'
#' @export
plot_sol_space <- function(sol_space) {
 
   if(nrow(sol_space$data) == 0){
        print("Solution space is epmty")
        return(0)
    }
    
    # 1 ######################### 
    pl1 <- sol_space$data %>%
        ggplot(aes(x = y)) +
        geom_line(aes(y = cumsum), color = "black") +
        geom_point(aes(y = cumsum,  color = group), size = 3) +
        geom_col(aes(y = median, fill = group)) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
        geom_text(aes(label = median %>% round(3), y = median), vjust = 0) +
        geom_hline(yintercept = 0.5, color = "red") +
        theme_light()
    
    # 2 ######################### 
    pl2 <- sol_space$data %>%
        ggplot(aes(y = id)) +
        geom_rect(
            fill = "gray",
            aes(
                xmin = lower,
                xmax = upper,
                ymin = y - 0.4,
                ymax = y + 0.4
            ) ,
            alpha = 0.3
        ) +
        geom_segment(aes(
            x = lower,
            xend = upper,
            y = y,
            yend = y
        ), linewidth = 1) +
        geom_point(aes(x = median , y = y, 
                       shape = if_else(valid == TRUE, "v", "NV"), 
                       color = if_else(valid == TRUE, "v", "NV"), 
                   size = 3)) +
        scale_shape_manual(values = c("v" = 16, "NV" = 124)) +  # 16 = circle, 3 = cross
        scale_color_manual(values = c("v" = "black", "NV" = "red")) +
        geom_vline(xintercept = 0.5,
                   color = "red",
                   linetype = "dashed") +
        geom_vline(xintercept = 0.0, color = "black") +
        geom_vline(xintercept = 1.0, color = "black") +
        
        labs(x = "frequency of solutions",
             y = "solution")   +
        theme_minimal() 
    
    # 3 ######################### 
    nn <- nrow(sol_space$data)
    df <- data.frame()
    for (j in 1:nn) {
        if (nn == 1) {
            comm_labels <- sol_space$M
        } else {
            comm_labels <- sol_space$M[, j]
        }
        community_size_dist <- table(comm_labels) |>
            sort(decreasing = TRUE)  |>
            unname() |> as.integer()
        
        df <- rbind(df,
                    data.frame(
                        cs = community_size_dist %>% round(0),
                        x = 1:length(community_size_dist),
                        y = rep(j, length(community_size_dist))
                    ))
    }
    
    pl3 <- df %>%
        ggplot( aes(x = x, y = y)) +
        geom_point(aes(
            size = cs,
            color = if_else(cs == 1, "single", "comm"),
            shape = if_else(cs == 1, "single", "comm")
        )) +
        scale_color_manual(values = c("black", "blue")) +
        scale_shape_manual(values = c(1, 5)) +
        ylim(0.4, nrow(sol_space$data)+0.4)+
        
        theme_minimal() +
        labs(x = "community size", y = "", solution = "") 
    
    #heatmap
    if (nrow(sol_space$data)>=2) {
        
    
    simil_df <- as.data.frame(sol_space$simil)
    simil_df$Partition1 <- rownames(simil_df)
    simil_long <- simil_df %>%
        pivot_longer(cols = -Partition1, names_to = "Partition2", values_to = "similarity")
    
    simil_long <- simil_long %>%
        mutate(Partition1 = factor(Partition1, levels = unique(Partition1)),
               Partition2 = factor(Partition2, levels = unique(Partition2)))
    
    pl4 <-  ggplot(simil_long, aes(x = Partition2, y = Partition1, fill = similarity)) +
        geom_tile() +
        #scale_fill_gradient(low = "white",  high = "darkgreen", limits = c(0, 1)) +
        scale_fill_distiller(palette = "PiYG", direction = 1, limits = c(0, 1), 
                             na.value = "white", 
                             guide = "colorbar") +
        labs(x = "solutions", y = "solutions", title = "Heatmap of Similarity Matrix") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              axis.text.y = element_text(hjust = 1))+ theme(aspect.ratio = 1.0)
    
    } else {pl4 = NA}
    
    return(list(pl1 = pl1, pl2 = pl2, pl3 = pl3, pl4 = pl4))
}


#' Plot Multiple Solutions from a Solution Space
#'
#' The function plot_all_solutions plots all the the solutions from a given solution space on a graph using a grid layout. 
#' Each plot highlights a different solution with communities highlighted in red. 
#'
#' @param g An `igraph` object representing the network graph.
#' @param sol_space A list containing the solution space. It should include:
#'   \describe{
#'     \item{data}{A matrix with rows corresponding to different solutions.}
#'     \item{M}{A matrix where each column corresponds to the membership of nodes for a given solution.}
#'   }
#' @param device A string specifying the output device. Possible values are "screen" (default) and "png".
#' @param filename A string specifying the file name for saving the plot if `device = "png"`. Default is `NULL`.
#' @param width Width of the PNG image in pixels. Default is 1600.
#' @param height Height of the PNG image in pixels. Default is 1600.
#' @param res Resolution of the PNG image in DPI. Default is 300.
#' 
#' @return NULL. This function is called for creating plots.
#' 
#' @examples
#' # Example usage (assume `g` is an igraph object and `sol_space` is defined):
#' # plot_solutions(g, sol_space)
#'
#' @export


plot_solutions <- function(g, ssp, 
                           add_node_labels = TRUE,
                           add_prob_labels = TRUE,
                           add_title = TRUE,
                           device = "screen", 
                           filename = NULL, 
                           width = 1600, height = 1600, res = 300) {
    
    
    # If device is PNG, open the PNG file
    if (device == "png" && !is.null(filename)) {
        png(filename, width = width, height = height, res = res)  
    }

    ns <- nrow(ssp$data)
    
    # Set up the plotting area
    par(mfrow = c(ceiling(sqrt(ns)), ceiling(sqrt(ns))),  # Grid layout
        mar = c(1.0, 0.1, 3.0, 0.1),  # margins: bottom, left, top, right
        oma = c(0, 0, 0, 0))        # No outer margins
    
    # Plot each solution
    node_positions <- igraph::layout.fruchterman.reingold(g)
    p_median = ssp$data$median %>% round(2)
    
    for (i in 1:ns) {
        # Extract membership information
        if (ns == 1){
            membership_i <- ssp$M
        } else {
            membership_i <- ssp$M[, i]  
        }
        
        
        
        # Plot the graph highlighting the i-th solution  

        plot(g, 
             vertex.label = if(add_node_labels) {V(g)$name} else {NA},
             layout = node_positions,
             vertex.size = if(add_node_labels) {30} else {15},
             vertex.color = if(add_node_labels) {"white"} else {"lightblue"},        
             mark.groups = split(1:vcount(g), membership_i),          
             mark.col = rgb(0.5, 0.5, 0.5, alpha = 0.1),   
             mark.border = 'red',              
             
        )
        # Add title if `add_title` is TRUE
        if (add_title) {
            mtext(paste("Solution", i), side = 3, line = 2, cex = 0.8)
        }
        
        # Add probability labels if `add_prob_labels` is TRUE
        if (add_prob_labels) {
            mtext(paste("p = ", p_median[i]), side = 3, line = 1, cex = 0.6)
        }
    }
    
    # If PNG device was opened, close it
    if (device == "png") {
        dev.off()
    }
    
    # Reset to default plotting layout (if plotting to the screen)
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}

#' @export
plot_sol_space_evolution <- function(ssp, confidence){
    x = (1 - confidence) / 2
    df <- ssp$log %>% 
        mutate(p = a/(a+b)) %>%
        mutate(up = qbeta(1 - x, a, b)) %>%
        mutate(lo = qbeta(x, a, b)) %>%
        mutate(s = as.factor(s))
    
    pl <- df %>% ggplot() +
        geom_line(aes(x = t, y = p, group = s, color = s)) +
        geom_ribbon(aes(x = t, ymin = lo, ymax = up, fill = s, ), linewidth = 3, alpha = 0.1) +
        ylim(0.0, 1.0) +
        theme_minimal()
    
    return(pl)
}
