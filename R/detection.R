


#' @export
solutions_space <-
    function(g,
             tmax,
             met='IM',
             comp_method='adjusted.rand'
             confidence = .95,
             param = NA) {
        M <- matrix(NA, nrow = vcount(g), ncol = 1)
        S <- matrix(0.0,  nrow = tmax, ncol = tmax)
        ns <- 0
        nn <- c()
        for (t in 1:tmax) {
            gs <- igraph::permute(g, sample(vcount(g)))
            comms <- switch(
                met,
                "IM" = igraph::infomap.community(gs, nb.trials = param),
                "WT" = igraph::walktrap.community(gs, steps = param),
                "LV" = igraph::cluster_louvain(gs, resolution = param),
                "LD" = igraph::cluster_leiden(gs, resolution = param),
                "LP" = igraph::label.propagation.community(gs)
            )
            
            membership <- comms$membership
            
            if (t == 1) {
                # first solution found
                M[, 1] <- membership[match(V(g)$name, V(gs)$name)]
                ns <- 1
                nn <- c(1)
                bayes_post <- data.frame(a = rep(2, 2), b = rep(2, 2))

            } else {
                # check if already esists
                for (i in 1:ns) {
                    sim_score <- igraph::compare(membership[match(V(g)$name, V(gs)$name)],
                                        M[, i],
                                        method = comp_method)
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
                    M <-
                        cbind(M, membership[match(V(g)$name, V(gs)$name)])
                    nn <- c(nn, 1)
                    #update
                    bayes_post$a[ns] <- 2
                    bayes_post$b <- t - bayes_post$a + 2
                    #add a new distribution to test unseen solutions
                    bayes_post <- rbind(bayes_post, data.frame(a = 1, b = t + 1))
                }
                bayes_post$b <- t - bayes_post$a + 2
                s = nrow(bayes_post)
                results <- bayes_post %>%    mutate(lower = NA, upper = NA, median = NA)
                

                
                for (i in 1:s) {
                    results$upper[i] <- qbeta(confidence, bayes_post$a[i], bayes_post$b[i])
                    results$lower[i] <-qbeta(1 - confidence, bayes_post$a[i], bayes_post$b[i])
                    results$median[i] <-qbeta(0.5, bayes_post$a[i], bayes_post$b[i])
                }
                results <- results %>% arrange(-median)
                results$group[1] <- 1
                for (i in 2:s) {
                    gap = results$lower[i] < results$upper[i - 1]
                    results$group[i] <- if_else(gap,
                                                results$group[i - 1] + 1,
                                                results$group[i - 1])
                }
                
            }#end if
            
            ## 
            ## TO DO exit for loop before tmax if results are (almost) unchanged 
            ## 
        }#end for
        
        results$y <- 1:nrow(results)
        results$id <- paste0("s", results$y)
        
        results$id[results$y  == nrow(results)] <- "New"
        results$group[results$id == "New"] <-results$group[max(results$group)] + 1
        
        results$group <- factor(results$group)
        
        
        results <- results %>% mutate(cumsum = cumsum(median))  %>%
            arrange(y)
        
        return(list(M = M[, order(-nn)], data = results))
    }
