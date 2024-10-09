
bayesian_update_model <- function(prior, selected){
    # a, b: Parameters of the prior Beta(a, b) distribution
    # i: index of  solution that will be updates
    posterior <- prior
    for (i in 1:nrow(prior)){
        if (i == selected){
            # increment success count for the solution found
            posterior$a[i] <- prior$a[i] + 1 
        } else {
            # increment failure count for all other solutions
            posterior$b[i] <- prior$b[i]  + 1
        }
    }
    return(posterior)
}

bayesian_add_new_solution <- function(prior){
    # updates the prior distribution adding a new solution
    # the number of trials is embedded in the prior
    t = prior$a[1] + prior$b[1] - 1
    # add a new beta for the new solution
    prior <- rbind(prior, data.frame(a = 1, b = t))
    return(prior)
}



#' Explore the Solution Space of Community Detection Algorithms
#'
#' The `solutions_space` function generates the solution space of a given community detection algorithm, i.e. the set of all independent solutions identified after a number of independent trials. By running multiple trials of a selected community detection method, it identifies one or more unique solutions, and calculates similarity between them. The function also computes the probability associated with each unique solution, using a Bayesian approach (to estimate the confidence of each solution and organizes the results based on confidence intervals)Beta Binomial model).
#'
#' @param g An iGraph object representing the network on which community detection will be performed.
#' @param n_trials An integer specifying the number of trials to run for community detection.
#' #' @param tau: A value between 0.0 and 1.0 specifying the threshold of stability of the solution space. When p_stable >= tau, the solution spaces is considered stable, i.e. the probability of a nwe solution to be found at the next trial is below tau. and thethe explorarion is terminated. Default value is 0.95.
#' @param met A string specifying the community detection method to use. Available options include:
#'   - `"IM"`: Infomap
#'   - `"WT"`: Walktrap
#'   - `"LV"`: Louvain
#'   - `"LD"`: Leiden
#'   - `"LP"`: Label Propagation
#'   - `"EV"`: Leading Eigenvector
#'   - `"EB"`: Edge Betweenness
#' @param shuffle A logical value indicating whether to shuffle the network (i.e. permute the order of nodes and edges) order before running the community detection algorithm. Default is `TRUE`.
#' @param comp_method A string specifying the method for comparing solutions. Available options:
#'   - `"ami"`: Adjusted Mutual Information (default)
#'   - `"ari"`: Adjusted Rand Index
#' @param confidence A numeric value representing the confidence level for estimating the confidence intervals of the solutions. Default is 0.95.
#' @param resolution A numeric value for the resolution parameter used by some community detection algorithms, such as Louvain and Leiden. Default is 1.0.
#' @param IM.nb.trials An integer specifying the number of trials for the Infomap algorithm. Default is 10.
#' @param WT.steps An integer specifying the number of steps for the Walktrap algorithm. Default is 3.
#'
#' @return A list-like object containing 3 items:
#'   - `M`: A matrix where each column represents a unique community detection solution. Rows are ordered as the nodes in the input graph. Columns are ordered by decreasing probability (median of confidence intervals).
#'   - `data`: A dataframe summarizing the features of each solution, including confidence intervals, similarity scores, and grouping based on confidence levels.
#'   - `simil`: A similarity matrix (using ARI) between the different solutions found.
#'
#' @details The function performs multiple trials of the chosen community detection method and compares the solutions using the selected comparison method (AMI or ARI). It checks if new solutions are unique by comparing them to previously found solutions. Unique solutions are added to the result set, and a Bayesian updating procedure is applied to compute confidence intervals for the posterior distribution of each solution. A similarity matrix between solutions is also computed using ARI.
#'
#' The function supports several popular community detection algorithms, and it allows for shuffling the node order between trials to introduce randomness and explore different solutions.
#'
#' @examples
#' # Run 10 trials of the Louvain method and explore the solution space
#' solution_space <- solutions_space(g, n_trials = 10, tau = 0.9, met = 'LV', comp_method = 'ami')
#' print(solution_space$data)
#'
#' @export
solutions_space <-
    function(g,
             n_trials = 100,
             tau = 0.95, 
             met = 'IM',
             shuffle = TRUE,
             comp_method = 'ami',#ari
             confidence = .95,
             resolution = 1.0,
             IM.nb.trials = 10,
             WT.steps = 3, 
             verbose = TRUE) {
        
        M <- matrix(NA, nrow = vcount(g), ncol = 1)
        rownames(M) <- V(g)$name
        S <- matrix(0.0,  nrow = n_trials, ncol = n_trials)
        ns <- 0
        
        prior <- data.frame(a = 1, b = 1) # uninformative prior: no trials, no info
        log = data.frame()
        
        for (t in 1:n_trials) {
            if (shuffle == TRUE) {
                gs <- igraph::permute(g, sample(vcount(g)))
            } else {
                gs <- g
            }
            
            comms <- switch(
                met,
                "IM" = igraph::infomap.community(gs, nb.trials = IM.nb.trials),
                "WT" = igraph::walktrap.community(gs, steps = WT.steps),
                "LV" = igraph::cluster_louvain(gs, resolution = resolution),
                "LD" = igraph::cluster_leiden(gs, resolution_parameter = resolution),
                "LP" = igraph::label.propagation.community(gs),
                "EV" = igraph::cluster_leading_eigen(gs),
                "EB" = igraph::cluster_edge_betweenness(gs)
            )
            
            membership <- comms$membership
            
            if (t == 1) {
                # first solution found
                M[, 1] <- membership[match(V(g)$name, V(gs)$name)]
                ns <- 1
                posterior <- bayesian_update_model(prior, 1)

            } else {
                # further trials: check if the solution is new
                for (i in 1:ns) {
                    sim_score <- switch(
                        comp_method,
                        "ari" = aricode::ARI(membership[match(V(g)$name, V(gs)$name)], M[, i]),
                        "ami" = aricode::AMI(membership[match(V(g)$name, V(gs)$name)], M[, i]),
                        "Invalid selection"
                    )
                    
                    if (sim_score == 1) {
                        # we already have this solution
                        posterior <- bayesian_update_model(prior, i)
 
                        break # no need to explore further
                    }
                }#end for
                
                if (sim_score < 1) {
                    # it's a new solution
                    ns <- ns + 1
                    M <- cbind(M, membership[match(V(g)$name, V(gs)$name)])
                    prior <- bayesian_add_new_solution(prior)
                    posterior <- bayesian_update_model(prior, ns)
 
                    
                }
                # bayesian update
                results <- posterior %>%
                    mutate(lower = NA,
                           upper = NA,
                           median = NA)
                
                s = nrow(results)
                
                # confidence intervals
                x = (1 - confidence) / 2
                for (i in 1:s) {
                    results$upper[i] <- qbeta(1 - x, posterior$a[i], posterior$b[i])
                    results$lower[i] <-
                        qbeta(x  , posterior$a[i], posterior$b[i])
                    results$median[i] <-
                        qbeta(0.5, posterior$a[i], posterior$b[i])
                }
                results <- results %>% arrange(-median)
                
            }#end if
            
            prior <- posterior
            
            # save the results at each iteration in the log dataframe
            log <- rbind(log, posterior %>% select(a, b) %>%
                             mutate(s = 1:nrow(posterior)))
            
            # exit the for loop before t_max in case the solution is stable
            p_stable = (t-ns)/(t+1)
            
            if (p_stable > tau){
                if (verbose) {print(paste("Exit at t = ", t))}
                 break
                }
        }#end for
        
        #results <- results %>% filter(a > 0) #remove empty lines
        results$y <- 1:nrow(results)
        if (shuffle == TRUE) {
            results$id <- sprintf("s%02d", results$y)
        } else {
            results$id <- sprintf("u%02d", results$y)
        }
        
        results$group <- NA
        grp <- 1
        results$group[1] <- grp
        if (ns > 1) {
            for (i in 2:ns) {
                gap = results$lower[i - 1] > results$upper[i]
                if (gap) {
                    grp <- grp + 1
                }
                results$group[i] <- grp
            }
        }
        
        results <- results %>%
            mutate(group = factor(group)) %>%
            mutate(cumsum = cumsum(median))  %>%
            filter(a > 0) %>%
            arrange(y)
        
        #calculate similarity matrix
        similarity_matrix <- matrix(NA, nrow = ns, ncol = ns)
        for (i in 1:ns) {
            for (j in i:ns) {
                if (i == j) {
                    next
                }
                similarity_score <- aricode::ARI(M[, i], M[, j])
                similarity_matrix[i, j] <- similarity_score
                similarity_matrix[j, i] <- similarity_score
            }
        }
        
        rownames(similarity_matrix) <- results$id
        colnames(similarity_matrix) <- results$id
        
        if (verbose) {print(paste("pstable = ", p_stable))}
        
        
        
        return(list(
            M = M[, order(-posterior$a)],
            data = results,
            simil = similarity_matrix, 
            log = log
        ))
    }




#' Calculate Normalized Co-occurrence Matrix
#'
#' The `co_occurrence` function computes a normalized co-occurrence matrix based on the 
#' solution space provided by the input object `ssp`. 
#' The output matrix indicates the extent to which pairs of nodes appear together in the same community 
#' across different trials, with weights normalized by the median values of the solutions.
#' 
#' The function filters out invalid trials based on the `valid` column in the 
#' `ssp$data` dataframe. It then computes the co-occurrence matrix `D` by iterating 
#' through each trial and identifying pairs of nodes that belong to the same 
#' community. The contribution of each trial to the co-occurrence matrix is 
#' weighted by the `median` values in the `ssp$data` dataframe, normalized by the 
#' total sum of medians.
#'
#' @param ssp A list containing the solution space, including community membership 
#' matrix `M` and associated data. The list should have at least the following elements (any further elements are ignored):
#'   \itemize{
#'     \item `M`: A matrix where each column represents a solution (trial), and 
#'     each row corresponds to a node. The matrix contains integers representing 
#'     community assignments.
#'     \item `data`: A data frame where rows correspond to the same trials as the 
#'     columns of `M`. It should include a `valid` logical vector indicating 
#'     whether a trial is valid, and a `median` column representing the weight 
#'     for each solution.
#'   }
#'   
#'
#' @return A symmetric matrix `D` where each entry represents the weighted co-occurrence count of node pairs across all trials. The matrix shows how often nodes were grouped together in the same community, scaled by the weights in `alpha`.
#'
#' @details The function uses the results of multiple community detection trials stored in `M`. For each trial, nodes that belong to the same community are identified, and their co-occurrence count is incremented by the corresponding value from `alpha`. The co-occurrence matrix is symmetric, reflecting the fact that if node A co-occurs with node B, then node B also co-occurs with node A.
#'
#' @examples
#' D <- co_occurrence(ssp)
#'
#' @export
co_occurrence <- function(ssp) {
    # calculates normalized co-occurrence matrix from solution space
    # 
    
    #keep only valid results
    keep_valid_results = ssp$data$valid  
    M <- ssp$M[, keep_valid_results]
    sspdata<-ssp$data[ keep_valid_results,]
    
    # weight of each solution 
    alpha <- sspdata$median / sum(sspdata$median) 
    
    # initialize co-occurrence matrix D
    n_nodes <- nrow(M)
    D <- matrix(0, nrow = n_nodes, ncol = n_nodes)
    rownames(D) <- rownames(M)
    colnames(D) <- rownames(M)
    
    # calculate co-occurrence
    n_solutions <- nrow(sspdata)
    for (t in (1:n_solutions)) {
        n_comms <- max(M[, t])
        for (k in 1:n_comms) {
            comm_members <- which(M[, t] == k)
            nc <- length(comm_members)
            for (i in 1:nc) {
                for (j in (i+1):nc) {
                    D[comm_members[j], comm_members[i]] <- D[comm_members[j], comm_members[i]] + alpha[t]
                    D[comm_members[i], comm_members[j]] <- D[comm_members[j], comm_members[i]]
                }
            }
        }
    }
    
    diag(D)<-1.0

    return (D)
}



#' Identify Consensus Communities from Co-occurrence Matrix
#'
#' The `consensus_communities` function identifies consensus communities from a co-occurrence matrix, grouping nodes into communities based on a threshold `p` for pairwise co-occurrence. It also calculates an uncertainty coefficient (`gamma`) for each node, reflecting how confidently each node belongs to its community. The function can handle outliers (single-node communities).
#'
#' @param D A symmetric co-occurrence matrix where each entry represents the pairwise co-occurrence of nodes across multiple trials.
#' @param p A numeric threshold for defining communities. Nodes are considered to be in the same community if their pairwise co-occurrence value is greater than `p`.
#' @param group_outliers A logical value indicating whether single-node communities (outliers) should be grouped together. Default is `FALSE`.
#' @param verbose A logical value indicating whether to print detailed progress information. Default is `FALSE`.
#'
#' @return A dataframe containing the following columns:
#'   - `name`: The name of each node.
#'   - `cons_comm_label`: The consensus community label assigned to the node.
#'   - `gamma`: The uncertainty coefficient for each node, calculated as `1 - mean(di)` over all nodes that co-occur at least once in the same community.
#'   - `comm_size`: The size of the community to which the node belongs.
#'   - `single`: A boolean indicating whether the node is part of a single-node community (outlier).
#'
#' @details The function processes the co-occurrence matrix by iteratively grouping nodes into communities based on the threshold `p`. For each identified community, it computes the uncertainty coefficient `gamma` for each node, which quantifies how strongly a node is tied to its community. The function can optionally group outliers.
#'
#' @examples
#' # Assume D is a co-occurrence matrix
#' results <- consensus_communities(D, p = 0.5, group_outliers = TRUE)
#' print(results)
#'
#' @export

consensus_communities <- function(D, p, group_outliers = FALSE, verbose = FALSE) {
    
    # definition of community: block within D in which dij > p 
    # this definition includes single node communities (outliers)
    
    # definition of uncertainty coefficient gamma: 
    #     (1-MEAN of di) over all nodes that are at least once in the same community
    
    results <- data.frame(name = colnames(D))
    results$done <- FALSE
    results$tmp_comm_label <- NA
    results$gamma <- NA
    results$comm_size <- NA
    results$single <- FALSE
    community_label <- 0
    nodes_to_process = nrow(results)
    
    # definition of community: block within D in which dij > p
    # this definition includes single node communities (outliers)
    
    # definition of uncertainty coefficient gamma:
    # (1-MEAN of di) over all nodes that are at least once in the same community
    
    results <- data.frame(name = colnames(D))
    results$done <- FALSE
    results$tmp_comm_label <- NA
    results$gamma <- NA
    results$comm_size <- NA
    results$single
    community_label <- 0
    nodes_to_process = nrow(results)
    
    while (nodes_to_process > 0)  {
        community_label <- community_label + 1
        
        #select a block with respect to threshold p, first row not done
        nodes_internal <- (D[which.max(results$done == FALSE),] > p)
        
        # calculate gamma for eachnode in the block
        gammas <- D[nodes_internal,]
        # ignore nodes that are never in the same community
        gammas[gammas == 0] <- NA
        
        if (sum(nodes_internal) > 1) {
            # a proper block
            results$gamma[nodes_internal] <- 1 - apply(gammas, 1, mean, na.rm = T)
        } else {
            # a single node
            results$gamma[nodes_internal] <- 1 - mean(gammas,  na.rm = T)
        }
        
        results$tmp_comm_label[nodes_internal] <- community_label
        results$comm_size[nodes_internal] <- sum(nodes_internal)
        results$done[nodes_internal] <-  TRUE
        nodes_to_process <- sum(results$done == FALSE)
    }
    
    results$gamma[is.na(results$gamma)] <- 0.0
    results$single[results$comm_size == 1] <- TRUE
    
    if (group_outliers) { results$tmp_comm_label[results$single] <- 0 }
    
    x <- results %>%
        group_by(tmp_comm_label) %>%
        summarize(n = n()) %>% arrange(-n) %>%
        mutate(cons_comm_label = row_number())
    
    results <- results %>%
        inner_join (x, by = 'tmp_comm_label') %>%
        select(name, cons_comm_label, gamma, comm_size, single)
    

    return(results)
}


