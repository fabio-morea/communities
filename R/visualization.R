#' This function builds a network (Gc) whose nodes represent communities within the original network (G). 
#' Gc can be utilized for in-depth analysis and visualization of the relationships among these communities.
#' 
#' @param g: The network to be analyzed. It should be an iGraph object with specific attributes:
#'   - V(g)$community: An integer value representing the community assignment for each node.
#'   - E(g)$w: A numeric vector containing edge weights. If the network is unweighted, set E(g)$w to 1.0.
#' 
#' @returns A network object (Gc) with the following attributes:
#'   - $membership: Stores the community labels.
#'   - E(Gc)$w: The sum of corresponding edge weights in the original network (G).
#'   - V(Gc)$size: the number of nodes of the original network G that are represented in a node of Gc.
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

#' The function layout_distance_gamma calculates a distance between pairs of nodes, that enhances the visualization,
#' grouping nodes of the same community, and spreading out different communities.
#'
#' @param g: The network to be analyzed. It should be an iGraph object with weights E(g)$weight (if the network is unweighted, set E(g)$weight to 1.0)
#' @param D: A matrix of co-occurrence as defined by library CCD.
#' @param eps: the distance between community members. typical value in the range 0.1 to 0.2.
#'
#' @returns An array of distances, ordered as E(g).
#'
#' @export
layout_distance_gamma <- function(g, D, eps) {
    if (!is.igraph(g)) {
        stop("g is not an igraph graph.")
    }
    if (any(is.na(E(g)$weight))) {
        stop("Edge attribute 'weight' contains NA values.")
    }
    if (any(is.na(D))) {
        stop("matrix D contains NA values.")
    }
    
    df <- as_long_data_frame(g) %>% select(from, to, weight)
    df$dist <- 0
    for (i in 1:nrow(df)) {
        df$dist[i] <- eps + (1 - eps) * df$weight[i] * D[df$from[i], df$to[i]]
    }
    return(df$dist)
}


#' The function layout_distance_gamma calculates a distance between pairs of nodes, that enhances the visualization,
#' grouping nodes of the same community, and spreading out different communities.
#'
#' @param g: The network to be analyzed. It should be an iGraph object with weights E(g)$weight (if the network is unweighted, set E(g)$weight to 1.0)
#' @param membership: The membership vector, ordered as V(g).
#' @param eps: the distance between community members. typical value in the range 0.1 to 0.2.
#'
#' @returns An array of distances, ordered as E(g).
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


#' The function plot_solution_space produces plots of the solution space
#' @param sol_space: a dataframe produced by explore_solution_space function
#' @returns A list of 2 plots
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
                        cs = community_size_dist,
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
        
        theme_minimal() +
        labs(x = "community", solution = "")
    
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
        scale_fill_gradient(low = "white",  high = "darkgreen", limits = c(0, 1)) +
        labs(x = "solutions", y = "solutions", title = "Heatmap of Similarity Matrix") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              axis.text.y = element_text(hjust = 1))+ theme(aspect.ratio = 1.0)
    
    } else {pl4 = NA}
    
    return(list(pl1 = pl1, pl2 = pl2, pl3 = pl3, pl4 = pl4))
}
