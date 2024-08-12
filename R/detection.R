


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
        nn <- c()
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
                nn <- c(1)
                bayes_post <- data.frame(a = 2, b = 2)

            } else {
                
                # check if already esists
                for (i in 1:ns) {
                    sim_score <- switch(
                        comp_method,
                        "ari" = aricode::ARI(membership[match(V(g)$name, V(gs)$name)],M[, i]),
                        "ami" = aricode::AMI(membership[match(V(g)$name, V(gs)$name)],M[, i]),
                        "Invalid selection")  
                    
                    if (sim_score == 1) {
                        #we already have this solution
                        nn[i] <- nn[i] + 1
                        bayes_post$a[i] <- bayes_post$a[i] + 1
                        bayes_post$b <- t - bayes_post$a + 2
                        break
                    }
                }#end for

                if (sim_score < 1) {
                    
                    #it's a new solution
                    ns <- ns + 1
                    M <- cbind(M, membership[match(V(g)$name, V(gs)$name)])
                    nn <- c(nn, 1)
                    
                    # add a new empty row 
                    bayes_post <- rbind(bayes_post, data.frame(a = NA, b = NA))
                    #update
                    bayes_post$a[ns] <- 2
                    bayes_post$b <- t - bayes_post$a + 2
                 
                }
                  
                results <- bayes_post %>%    mutate(lower = NA, upper = NA, median = NA)
                s = nrow(bayes_post)

                # confidence intervals
                x = (1-confidence)/2
                for (i in 1:s) {
                    results$upper[i] <- qbeta(1-x,  bayes_post$a[i], bayes_post$b[i])
                    results$lower[i] <- qbeta(x  ,    bayes_post$a[i], bayes_post$b[i]) 
                    results$median[i] <-qbeta(0.5, bayes_post$a[i], bayes_post$b[i])
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
        
        return(list(M = M[, order(-nn)], data = results))
    }
