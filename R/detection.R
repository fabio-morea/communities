
bayesian_update_1 <- function(prior, i){
    posterior <- prior
    # n increases by 1, hence b increases by 1 for all
    posterior$b <- prior$b + 1
    # except for the i-th element 
    posterior$b[i] <- prior$b[i]  
    
    # k increases by 1 only for the ith element
    posterior$a <- prior$a
    posterior$a[i] <- posterior$a[i] + 1  
    
    return(posterior)
}

bayesian_new_1 <- function(prior){
    n = prior$a[1] + prior$b[1]-2
    posterior <- rbind(prior, data.frame(a = 2, b = n))
    return(posterior)
}


#' @export
solutions_space <-
    function(g,
             tmax,
             met='IM',
             comp_method='ami',#ari
             confidence = .95,
             resolution = 1.0, IM.nb.trials = 10, WT.steps=3) {
        M <- matrix(NA, nrow = vcount(g), ncol = 1)
        S <- matrix(0.0,  nrow = tmax, ncol = tmax)
        ns <- 0   

        prior <- data.frame(a = 1, b = 1) # no trials, no info
        
        for (t in 1:tmax) {
            gs <- igraph::permute(g, sample(vcount(g)))
            comms <- switch(
                met,
                "IM" = igraph::infomap.community(gs, nb.trials=IM.nb.trials),  
                "WT" = igraph::walktrap.community(gs, steps = WT.steps),
                "LV" = igraph::cluster_louvain(gs, resolution = resolution),
                "LD" = igraph::cluster_leiden(gs, resolution_parameter = resolution),
                "LP" = igraph::label.propagation.community(gs)
            )
            
            membership <- comms$membership
            
            
            
            if (t == 1) {
                # first solution found
                M[, 1] <- membership[match(V(g)$name, V(gs)$name)]
                ns <- 1
                posterior <- bayesian_update_1(prior, i=1)

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
                        break
                    }
                }#end for

                if (sim_score < 1) {
                    
                    #it's a new solution
                    ns <- ns + 1
                    M <- cbind(M, membership[match(V(g)$name, V(gs)$name)])

                    # add a new solution
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

        }#end for
        
        
        
        results <- results %>% filter(a > 0) #remove empty line
        results$y <- 1:nrow(results)
        results$id <- paste0("s", results$y)
        
        # to add: check for non-valid communities k=1, k=n or disconnected
        # 
        
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
