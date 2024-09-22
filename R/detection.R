
bayesian_update_1 <- function(prior, selected){
    # a, b: Parameters of the prior Beta(a, b) distribution
    # i: index of  solution that will be updates
    posterior <- prior
    for (i in 1:nrow(prior)){
        if (i == selected){
            posterior$a[i] <- prior$a[i] + 1
        } else {
            posterior$b[i] <- prior$b[i]  + 1
        }
    }
    return(posterior)
}

bayesian_new_1 <- function(prior){
    posterior <- prior
    n = posterior$a[1] + posterior$b[1] - 2
    posterior <- rbind(posterior, data.frame(a = 2, b = n))
    posterior$b <- posterior$b + 1 
    return(posterior)
}



#' Explore the Solution Space of Community Detection Algorithms
#'
#' The `solutions_space` function generates the solution space of a given community detection algorithm, i.e. the set of all independent solutions identified after a number of independent trials. By running multiple trials of a selected community detection method, it identifies one or more unique solutions, and calculates similarity between them. The function also computes the probability associated with each unique solution, using a Bayesian approach (to estimate the confidence of each solution and organizes the results based on confidence intervals)Beta Binomial model).
#'
#' @param g An iGraph object representing the network on which community detection will be performed.
#' @param n_trials An integer specifying the number of trials to run for community detection.
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
#' solution_space <- solutions_space(g, n_trials = 10, met = 'LV', comp_method = 'ami')
#' print(solution_space$data)
#'
#' @export
solutions_space <-
    function(g,
             n_trials,
             met = 'IM',
             shuffle = TRUE,
             comp_method = 'ami',
             #ari
             confidence = .95,
             resolution = 1.0,
             IM.nb.trials = 10,
             WT.steps = 3) {
        M <- matrix(NA, nrow = vcount(g), ncol = 1)
        S <- matrix(0.0,  nrow = n_trials, ncol = n_trials)
        ns <- 0
        
        prior <- data.frame(a = 1, b = 1) # no trials, no info
        
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
                posterior <- bayesian_update_1(prior, 1)
                
            } else {
                # check if already exists
                for (i in 1:ns) {
                    sim_score <- switch(
                        comp_method,
                        "ari" = aricode::ARI(membership[match(V(g)$name, V(gs)$name)], M[, i]),
                        "ami" = aricode::AMI(membership[match(V(g)$name, V(gs)$name)], M[, i]),
                        "Invalid selection"
                    )
                    
                    if (sim_score == 1) {
                        #we already have this solution
                        posterior <- bayesian_update_1(posterior, i)
                        break # no need to explore further
                    }
                }#end for
                
                if (sim_score < 1) {
                    #it's a new solution
                    ns <- ns + 1
                    M <-
                        cbind(M, membership[match(V(g)$name, V(gs)$name)])
                    posterior <- bayesian_new_1(posterior)
                }
                
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
        
        return(list(
            M = M[, order(-posterior$a)],
            data = results,
            simil = similarity_matrix
        ))
    }
