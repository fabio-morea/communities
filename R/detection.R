
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


#' @export
solutions_space <-
    function(g,
             n_trials,
             met='IM',
             shuffle = TRUE,
             comp_method='ami',#ari
             confidence = .95,
             resolution = 1.0, IM.nb.trials = 10, WT.steps=3) {
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
                "IM" = igraph::infomap.community(gs, nb.trials=IM.nb.trials),  
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
                        "ari" = aricode::ARI(membership[match(V(g)$name, V(gs)$name)],M[, i]),
                        "ami" = aricode::AMI(membership[match(V(g)$name, V(gs)$name)],M[, i]),
                        "Invalid selection")  
                    
                    if (sim_score == 1) {
                        #we already have this solution
                        posterior <- bayesian_update_1(posterior, i)
                        break # no need to explore further
                    }
                }#end for

                if (sim_score < 1) {
                    #it's a new solution
                    ns <- ns + 1
                    M <- cbind(M, membership[match(V(g)$name, V(gs)$name)])
                    posterior <- bayesian_new_1(posterior)
                }
                  
                results <- posterior %>%  
                    mutate(lower = NA, upper = NA, median = NA)
                
                s = nrow(results)

                # confidence intervals
                x = (1-confidence)/2
                for (i in 1:s) {
                    results$upper[i] <- qbeta(1-x, posterior$a[i], posterior$b[i])
                    results$lower[i] <- qbeta(x  , posterior$a[i], posterior$b[i]) 
                    results$median[i] <-qbeta(0.5, posterior$a[i], posterior$b[i])
                }
                results <- results %>% arrange(-median)
                
            }#end if 
            
            prior <- posterior 

        }#end for
        
        #results <- results %>% filter(a > 0) #remove empty lines
        results$y <- 1:nrow(results)
        if(shuffle==TRUE){
            results$id <- sprintf("s%02d", results$y)
        } else {
            results$id <- sprintf("u%02d", results$y)
        }

        results$group <- NA
        grp<-1
        results$group[1] <- grp
        if (ns > 1){
            for (i in 2:ns) {
                gap = results$lower[i - 1] > results$upper[i]
                if (gap) {grp <- grp + 1}
                results$group[i] <- grp
            }
        }

        results <- results %>% 
            mutate(group = factor(group))%>%
            mutate(cumsum = cumsum(median))  %>%
            filter(a>0)%>%
            arrange(y)
        
        #calculate similarity matrix
        similarity_matrix <- matrix(NA, nrow = ns, ncol = ns)
        for (i in 1:ns) {
            for (j in i:ns) {
                if (i==j){next}
                similarity_score <- aricode::ARI(M[, i], M[, j])
                similarity_matrix[i, j] <- similarity_score
                similarity_matrix[j, i] <- similarity_score
            }
        }
        
        rownames(similarity_matrix) <- results$id
        colnames(similarity_matrix) <- results$id
        
        return(list(M = M[, order(-posterior$a)], data = results, simil = similarity_matrix))
    }
