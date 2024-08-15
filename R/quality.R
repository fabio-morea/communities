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


#' @export
quality_check <- function(g, sol_space){
    M <- sol_space$M 
    
    sol_space$data$valid <- TRUE
    
    n_solutions = ncol(M)
    
    for (i in 1:n_solutions){
        labels <- M[,i] %>% unique() %>% sort()
        sol_space$data$k[i] <- length(labels)
        sol_space$data$k_n[i] <- sol_space$data$k[i] / vcount(g)
        
        # modularity: need to use c_membership + 1 to handle community label 0
        sol_space$data$mod[i] <- igraph::modularity (g,  M[,i] + 1)
        
        # mixing parameter mu
        sol_space$data$mu[i] <- communities::empirical_mu(g, M[,i])
        
        # communities are internally connected
        sol_space$data$int_conn[i]  <- max(communities::internally_connected(g, M[,i]))
        
        # validity
        if (sol_space$data$mu[i] > 0.5){sol_space$data$valid[i] <- FALSE}
        if (sol_space$data$k[i] == 1){sol_space$data$valid[i] <- FALSE}
        if (sol_space$data$k_n[i] == 1){sol_space$data$valid[i] <- FALSE}
        if (sol_space$data$int_conn[i] > 1 ){sol_space$data$valid[i] <- FALSE}
    }

    return(sol_space)
}