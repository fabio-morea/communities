#' Build a Community Network from an iGraph Object
#'
#' Constructs a meta-network where each node represents a community from the 
#' original network. Edges between communities represent aggregated connections
#' between nodes in those communities.
#'
#' @param g An igraph object representing the original network. Required attributes:
#'   \itemize{
#'     \item \code{V(g)$community}: Integer vector of community assignments
#'     \item \code{E(g)$weight}: Numeric vector of edge weights (use 1.0 for unweighted)
#'   }
#'
#' @return An igraph object representing the community network with attributes:
#'   \itemize{
#'     \item \code{V(gc)$id}: Community identifiers
#'     \item \code{V(gc)$size}: Number of nodes in each community
#'     \item \code{E(gc)$weight}: Aggregated edge weights between communities
#'   }
#'
#' @details Aggregates nodes and edges by community membership to create a 
#'   higher-level representation of network structure.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- make_ring(10)
#' V(g)$community <- rep(1:2, each = 5)
#' E(g)$weight <- 1.0
#' gc <- make_community_network(g)
#' }
#'
#' @export
make_community_network <- function(g) {
    # Aggregate edges by community membership
    edges_list <- g %>%
        as_long_data_frame() %>%
        select(from_community, to_community, weight) %>%
        group_by(from_community, to_community) %>%
        summarize(weight = sum(weight), .groups = "drop")
    
    # Create community-level graph
    gc <- graph_from_data_frame(edges_list, directed = FALSE)
    V(gc)$id <- seq_len(vcount(gc))
    
    # Calculate community sizes
    comms <- data.frame(label = V(gc)$name)
    comms$size <- vapply(comms$label, function(label) {
        sum(V(g)$community == label)
    }, FUN.VALUE = integer(1))
    
    V(gc)$size <- comms$size
    
    return(gc)
}


#' Calculate Distances for Community-Based Network Layout
#'
#' Computes edge distances that pull nodes in the same community closer together
#' while pushing nodes in different communities farther apart. Useful for layout
#' algorithms that accept distance matrices.
#'
#' @param g An igraph object with edge weights stored in \code{E(g)$weight}.
#'   For unweighted graphs, set all weights to 1.0.
#' @param membership Numeric vector of community assignments, one per node,
#'   ordered as \code{V(g)}.
#' @param eps Numeric. Minimum distance between nodes in the same community.
#'   Default is 0.02. Typical values: 0.1 to 0.2.
#'
#' @return Numeric vector of distances, one per edge in \code{E(g)}, ordered
#'   according to the edge list.
#'
#' @details Distance formula for edge \eqn{(u,v)}:
#'   \deqn{d_{uv} = \epsilon + (1 - \epsilon) \cdot w_{uv} \cdot I(c_u = c_v)}
#'   where \eqn{I(\cdot)} is the indicator function and \eqn{w_{uv}} is the
#'   edge weight.
#'
#' @examples
#' \dontrun{
#' g <- make_ring(10)
#' E(g)$weight <- runif(ecount(g))
#' membership <- rep(1:2, each = 5)
#' distances <- layout_distance_comm(g, membership, eps = 0.15)
#' }
#'
#' @export
layout_distance_comm <- function(g, membership, eps = 0.02) {
    # Input validation
    if (!is.igraph(g)) {
        stop("g must be an igraph object.")
    }
    if (any(is.na(E(g)$weight))) {
        stop("Edge attribute 'weight' contains NA values.")
    }
    
    # Extract edge information
    df <- as_long_data_frame(g) %>% 
        select(from, to, weight)
    
    # Calculate distance for each edge
    df$dist <- vapply(seq_len(nrow(df)), function(i) {
        same_comm <- membership[df$from[i]] == membership[df$to[i]]
        eps + (1 - eps) * df$weight[i] * same_comm
    }, FUN.VALUE = numeric(1))
    
    return(df$dist)
}


#' Plot Solution Space Diagnostics
#'
#' Generates diagnostic plots for exploring community detection solution spaces.
#' Visualizes posterior probabilities, credible intervals, community size 
#' distributions, and inter-solution similarity.
#'
#' @param sol_space List output from \code{solutions_space_DM()}, containing:
#'   \itemize{
#'     \item \code{partitions}: Matrix of community assignments
#'     \item \code{probabilities}: Data frame with posterior estimates
#'   }
#' @param qc Data frame from \code{quality_check()} with \code{valid} column
#'   indicating solution validity.
#'
#' @return List containing up to 4 ggplot objects:
#'   \itemize{
#'     \item \code{pl1}: Posterior probabilities with credible intervals
#'     \item \code{pl2}: Solution space visualization with validity markers
#'     \item \code{pl3}: Community size distributions per solution
#'     \item \code{pl4}: Pairwise solution similarity heatmap (NA if only 1 solution)
#'   }
#'
#' @examples
#' \dontrun{
#' ssp <- solutions_space_DM(g, n_trials = 100)
#' qc <- quality_check(g, ssp)
#' plots <- plot_sol_space(ssp, qc)
#' print(plots$pl1)
#' }
#'
#' @export
plot_sol_space <- function(sol_space, qc) {
    # Check for empty solution space
    if (ncol(sol_space$partitions) == 0) {
        message("Solution space is empty.")
        return(NULL)
    }
    
    # Plot 1: Posterior probabilities with error bars
    pl1 <- sol_space$probabilities %>%
        ggplot(aes(x = id)) +
        geom_line(aes(y = phat), color = "black") +
        geom_point(aes(y = phat), size = 3) +
        geom_errorbar(aes(ymin = pmin, ymax = pmax), width = 0.2) +
        geom_hline(yintercept = 0.5, color = "red", linetype = "dashed") +
        labs(
            x = "Solution",
            y = "Posterior probability",
            title = "Solution probabilities with credible intervals"
        ) +
        theme_light()
    
    # Plot 2: Solution space with validity indicators
    pl2 <- sol_space$probabilities %>%
        mutate(y = row_number()) %>%
        ggplot(aes(y = id)) +
        geom_rect(
            aes(xmin = pmin, xmax = pmax, ymin = y - 0.4, ymax = y + 0.4),
            fill = "gray", alpha = 0.3
        ) +
        geom_segment(
            aes(x = pmin, xend = pmax, y = y, yend = y),
            linewidth = 1
        ) +
        geom_point(
            aes(
                x = phat, y = y,
                shape = if_else(qc$valid, "Valid", "Invalid"),
                color = if_else(qc$valid, "Valid", "Invalid")
            ),
            size = 3
        ) +
        scale_shape_manual(
            name = "Status",
            values = c("Valid" = 16, "Invalid" = 18)
        ) +
        scale_color_manual(
            name = "Status",
            values = c("Valid" = "black", "Invalid" = "red")
        ) +
        geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
        geom_vline(xintercept = c(0, 1), color = "black") +
        labs(
            x = "Frequency of solutions",
            y = "Solution",
            title = "Solution space with validity markers"
        ) +
        theme_minimal()
    
    # Plot 3: Community size distributions
    n_solutions <- nrow(sol_space$probabilities)
    
    # Ensure partitions is a matrix
    M <- sol_space$partitions
    if (is.null(ncol(M))) M <- matrix(M, ncol = 1)
    
    # Build data frame for community size plotting
    df <- map_dfr(seq_len(n_solutions), function(j) {
        comm_labels <- if (ncol(M) == 1) M[, 1] else M[, j]
        
        community_size_dist <- table(comm_labels) %>%
            sort(decreasing = TRUE) %>%
            unname() %>%
            as.integer()
        
        tibble(
            comm_size = as.integer(community_size_dist),
            comm_index = seq_along(community_size_dist),
            solution_id = j,
            is_singleton = comm_size == 1L
        )
    })
    
    pl3 <- ggplot(df, aes(x = comm_index, y = solution_id)) +
        geom_point(aes(size = comm_size), shape = 1, color = "black", stroke = 0.8) +
        scale_size_continuous(
            range = c(min(df$comm_size), max(df$comm_size)),
            breaks = function(lims) {
                by <- max(1, floor((lims[2] - 1) / 5))
                seq(1, lims[2], by = by)
            },
            labels = scales::number_format(accuracy = 1),
            name = "Community size"
        ) +
        ylim(0.5, n_solutions + 0.5) +
        labs(
            x = "Community index (sorted by size)",
            y = "Solution",
            title = "Community size distribution per solution"
        ) +
        theme_minimal()
    
    # Plot 4: Similarity heatmap (only if multiple solutions exist)
    pl4 <- NA
    if (ncol(sol_space$partitions) >= 2) {
        # Calculate similarity matrix
        simil <- similarity_matrix_nmi(sol_space$partitions)
        
        # Convert to long format
        simil_long <- simil %>%
            as.data.frame() %>%
            tibble::rownames_to_column("partition_1") %>%
            pivot_longer(
                -partition_1,
                names_to = "partition_2",
                values_to = "similarity"
            ) %>%
            filter(similarity >= 0)
        
        # Ensure consistent factor ordering
        all_parts <- colnames(simil)
        simil_long <- simil_long %>%
            mutate(
                partition_1 = factor(partition_1, levels = all_parts),
                partition_2 = factor(partition_2, levels = all_parts)
            )
        
        # Create heatmap
        pl4 <- ggplot(simil_long, aes(x = partition_2, y = partition_1, fill = similarity)) +
            geom_tile(color = "grey70", linewidth = 0.5) +
            geom_text(
                aes(label = sprintf("%.2f", similarity)),
                size = 3, color = "black"
            ) +
            scale_fill_gradient(
                low = "red", high = "green",
                limits = c(0, 1),
                na.value = "white",
                guide = "colorbar"
            ) +
            labs(
                x = "Solution",
                y = "Solution",
                title = "Partition Similarity (NMI)"
            ) +
            theme_minimal(base_size = 14) +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(hjust = 1),
                plot.title = element_text(hjust = 0.5, face = "bold"),
                axis.title.x = element_text(margin = margin(t = 10)),
                axis.title.y = element_text(margin = margin(r = 10)),
                legend.position = "right",
                legend.title = element_blank(),
                aspect.ratio = 1
            )
    }
    
    return(list(pl1 = pl1, pl2 = pl2, pl3 = pl3, pl4 = pl4))
}


#' Plot Multiple Solutions from Solution Space
#'
#' Creates a grid of network visualizations showing different community detection
#' solutions. Each panel displays the same network with communities highlighted
#' according to one solution.
#'
#' @param g An igraph object representing the network.
#' @param ssp List output from \code{solutions_space_DM()}, containing:
#'   \itemize{
#'     \item \code{partitions}: Matrix of community assignments
#'     \item \code{probabilities}: Data frame with posterior probabilities
#'   }
#' @param add_node_labels Logical. Display node names? Default: TRUE.
#' @param add_prob_labels Logical. Display solution probabilities? Default: TRUE.
#' @param add_title Logical. Display solution numbers as titles? Default: TRUE.
#' @param device Character. Output device: "screen" or "png". Default: "screen".
#' @param filename Character. File path if device = "png". Default: NULL.
#' @param width Numeric. PNG width in pixels. Default: 1600.
#' @param height Numeric. PNG height in pixels. Default: 1600.
#' @param res Numeric. PNG resolution in DPI. Default: 300.
#'
#' @return NULL (called for side effects: plotting).
#'
#' @details Uses Fruchterman-Reingold layout with consistent node positions
#'   across all panels. Communities are highlighted with semi-transparent
#'   colored regions.
#'
#' @examples
#' \dontrun{
#' g <- make_ring(10)
#' ssp <- solutions_space_DM(g, n_trials = 50)
#' plot_solutions(g, ssp)
#' # Save to file:
#' plot_solutions(g, ssp, device = "png", filename = "solutions.png")
#' }
#'
#' @export
plot_solutions <- function(g, ssp,
                           add_node_labels = TRUE,
                           add_prob_labels = TRUE,
                           add_title = TRUE,
                           device = "screen",
                           filename = NULL,
                           width = 1600,
                           height = 1600,
                           res = 300) {
    
    # Open PNG device if requested
    if (device == "png" && !is.null(filename)) {
        png(filename, width = width, height = height, res = res)
    }
    
    # Determine number of solutions
    n_solutions <- ncol(ssp$partitions)
    
    # Set up grid layout
    n_cols <- ceiling(sqrt(n_solutions))
    n_rows <- ceiling(n_solutions / n_cols)
    par(
        mfrow = c(n_rows, n_cols),
        mar = c(1.0, 0.1, 3.0, 0.1),
        oma = c(0, 0, 0, 0)
    )
    
    # Compute consistent layout across all solutions
    node_positions <- igraph::layout.fruchterman.reingold(g)
    posterior_probs <- round(ssp$probabilities$phat, 3)
    
    # Plot each solution
    for (i in seq_len(n_solutions)) {
        # Extract membership for this solution
        membership_i <- if (n_solutions == 1) {
            ssp$partitions
        } else {
            ssp$partitions[, i]
        }
        
        # Create plot
        plot(
            g,
            vertex.label = if (add_node_labels) V(g)$name else NA,
            layout = node_positions,
            vertex.size = if (add_node_labels) 30 else 15,
            vertex.color = if (add_node_labels) "white" else "lightblue",
            mark.groups = split(seq_len(vcount(g)), membership_i),
            mark.col = rgb(0.5, 0.5, 0.5, alpha = 0.1),
            mark.border = "red"
        )
        
        # Add title if requested
        if (add_title) {
            mtext(paste("Solution", i), side = 3, line = 2, cex = 0.8)
        }
        
        # Add probability label if requested
        if (add_prob_labels) {
            mtext(paste("p =", posterior_probs[i]), side = 3, line = 1, cex = 0.6)
        }
    }
    
    # Close PNG device if opened
    if (device == "png") {
        dev.off()
    }
    
    # Reset plotting parameters
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
    
    invisible(NULL)
}


#' Plot Solution Space Evolution Over Trials
#'
#' Visualizes how posterior probabilities of different solutions evolve across
#' trials. Shows convergence behavior and relative stability of solutions.
#'
#' @param ssp List output from \code{solutions_space_DM()}, must include
#'   \code{ssp$log$prob_long} with trial-by-trial probability estimates.
#' @param show_ci Logical. Display credible intervals as ribbons? Default: TRUE.
#' @param smooth Logical. Add LOESS smoothing curves? Default: FALSE.
#'
#' @return A ggplot object showing probability evolution.
#'
#' @details Requires \code{solutions_space_DM()} to be run with full logging
#'   enabled. Credible intervals show uncertainty in probability estimates.
#'   Smoothing can help identify trends in noisy convergence patterns.
#'
#' @examples
#' \dontrun{
#' ssp <- solutions_space_DM(g, n_trials = 200)
#' p <- plot_sol_space_evolution(ssp, show_ci = TRUE, smooth = TRUE)
#' print(p)
#' }
#'
#' @export
plot_sol_space_evolution <- function(ssp, show_ci = TRUE, smooth = FALSE) {
    # Validate input
    if (is.null(ssp$log$prob_long)) {
        stop("ssp$log$prob_long not found. Re-run solutions_space_DM() with full logging.")
    }
    
    suppressPackageStartupMessages({
        library(ggplot2)
        library(dplyr)
    })
    
    # Prepare data
    df <- ssp$log$prob_long %>%
        mutate(solution = factor(solution, levels = unique(solution)))
    
    # Build plot
    p <- ggplot(df, aes(x = trial, y = phat, color = solution)) +
        {
            if (show_ci) {
                geom_ribbon(
                    aes(ymin = pmin, ymax = pmax, fill = solution),
                    alpha = 0.15, color = NA
                )
            }
        } +
        geom_line(linewidth = 0.8, na.rm = TRUE) +
        {
            if (smooth) {
                geom_smooth(
                    se = FALSE, linewidth = 0.6, linetype = 3,
                    method = "loess", span = 0.4
                )
            }
        } +
        scale_y_continuous(
            limits = c(0, 1),
            expand = expansion(mult = c(0, 0.02))
        ) +
        labs(
            x = "Trial",
            y = "Posterior probability",
            color = "Solution",
            fill = "Solution",
            title = "Evolution of solution probabilities by trial"
        ) +
        theme_minimal(base_size = 13) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = "right"
        )
    
    return(p)
}