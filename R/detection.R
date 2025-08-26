
# bayesian_update_model <- function(prior, selected){
#     # a, b: Parameters of the prior Beta(a, b) distribution
#     # i: index of  solution that will be updates
#     posterior <- prior
#     for (i in 1:nrow(prior)){
#         if (i == selected){
#             # increment success count for the solution found
#             posterior$a[i] <- prior$a[i] + 1 
#         } else {
#             # increment failure count for all other solutions
#             posterior$b[i] <- prior$b[i]  + 1
#         }
#     }
#     return(posterior)
# }
# 
# bayesian_add_new_solution <- function(prior){
#     # updates the prior distribution adding a new solution
#     # the number of trials is embedded in the prior
#     t = prior$a[1] + prior$b[1] - 1
#     # add a new beta for the new solution
#     prior <- rbind(prior, data.frame(a = 1, b = t))
#     return(prior)
# }



#' Explore the Solution Space of Community Detection Algorithms
#'#' Explore the Solution Space of Community Detection Algorithms
#'
#' Repeatedly runs a community detection algorithm on a graph under randomized conditions
#' (e.g., vertex permutation) to explore the diversity of partitions generated. Tracks 
#' unique solutions, models posterior probabilities via a Dirichlet–Multinomial framework, 
#' and applies a precision stopping rule based on Beta credible intervals.
#'
#' @param g An `igraph` object. The graph whose community structure is to be analyzed.
#' @param n_trials Integer. Maximum number of trials to run (default: 200).
#' @param met Character. Community detection method. One of: `"IM"` (Infomap), `"WT"` (Walktrap), 
#'   `"LV"` (Louvain), `"LD"` (Leiden), `"LP"` (Label Propagation), `"EV"` (Leading Eigenvector), 
#'   `"EB"` (Edge Betweenness). Default is `"IM"`.
#' @param delta Numeric. Precision threshold: maximum allowed width of Beta credible intervals 
#'   for solution probabilities. Default is `0.05`.
#' @param confidence Numeric. Confidence level for Beta credible intervals (default: 0.95).
#' @param gamma0 Numeric. Symmetric Dirichlet prior parameter (default: 1).
#' @param seed Optional integer. If set, uses this value to seed the random number generator.
#' @param verbose Logical. If `TRUE`, prints progress at each trial. Default is `TRUE`.
#'
#' @return A list with three components:
#' \describe{
#'   \item{`partitions`}{A matrix (`nv` × `ns`) of unique partitions found, one per column. Rows are nodes.}
#'   \item{`probabilities`}{A data frame summarizing posterior probabilities for each solution:
#'     \code{id}, \code{count}, \code{phat}, \code{plower}, \code{pupper}.}
#'   \item{`log`}{A list of convergence diagnostics and trial-by-trial counts, including:
#'     \code{stop_trial}, \code{stop_reason}, \code{counts_matrix}, \code{counts_long},
#'     and the width of credible intervals at each step.}
#' }
#'
#' @details
#' - Uses name-based remapping to align partition membership vectors to the original graph order.
#' - Posterior inference is based on a symmetric Dirichlet prior and marginal Beta distributions.
#' - Stops early if the width of all marginal CIs falls below `delta`.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- make_ring(10)
#' result <- solutions_space_DM(g, n_trials = 100, met = "IM", seed = 42)
#' print(result$probabilities)
#' }
#'
#' @export
#' 

# ============================================================
# solutions_space(): Explore solution space  
# ============================================================
solutions_space_DM <- function(
        g,
        n_trials   = 200,
        met        = "IM",      # IM, WT, LV, LD, LP, EV, EB
        delta      = 0.05,      # precision rule: max marginal CI width
        confidence = 0.95,      # Beta marginal credible interval
        gamma0     = 1,         # symmetric Dirichlet prior
        seed       = NULL,
        verbose    = TRUE
) {
    if (!requireNamespace("igraph", quietly = TRUE))
        stop("Package 'igraph' is required.")
    if (!requireNamespace("aricode", quietly = TRUE))
        stop("Package 'aricode' is required.")
    if (!inherits(g, "igraph")) stop("`g` must be an igraph object.")
    nv <- igraph::vcount(g)
    if (nv < 2L) stop("Graph must have at least 2 vertices.")
    if (is.null(igraph::V(g)$name)) igraph::V(g)$name <- as.character(seq_len(nv))
    if (!is.null(seed)) set.seed(seed)
    
    # Internal settings (fixed)
    shuffle      <- TRUE
    im_trials    <- 10
    wt_steps     <- 3
    resolution   <- 1.0
    
    # Internal helpers
    eq_1 <- function(x, y) isTRUE(all.equal(aricode::ARI(x, y), 1))
    beta_ci <- function(counts) {
        ns <- length(counts)
        if (ns == 1L) {
            a <- counts + gamma0; b <- gamma0
            list(
                phat  = a / (a + b),
                lower = stats::qbeta((1 - confidence) / 2, a, b),
                upper = stats::qbeta((1 + confidence) / 2, a, b)
            )
        } else {
            a <- counts + gamma0
            A <- sum(a)
            list(
                phat  = a / A,
                lower = stats::qbeta((1 - confidence) / 2, a, A - a),
                upper = stats::qbeta((1 + confidence) / 2, a, A - a)
            )
        }
    }
    
    # State
    M <- matrix(numeric(0), nrow = nv, ncol = 0)
    counts <- numeric(0)
    names(counts) <- character(0)
    rownames(M) <- igraph::V(g)$name
    counts_by_trial <- vector("list", n_trials)
    max_ci_width <- numeric(n_trials)
    stop_trial <- NA_integer_
    stop_reason <- NA_character_
    
    for (t in seq_len(n_trials)) {
        perm <- if (shuffle) sample(nv) else seq_len(nv)
        gs <- igraph::permute(g, perm)
        
        memb_perm <- switch(
            met,
            "IM" = igraph::infomap.community(gs, nb.trials = im_trials)$membership,
            "WT" = igraph::walktrap.community(gs, steps = wt_steps)$membership,
            "LD" = igraph::cluster_leiden(gs, resolution_parameter = resolution)$membership,
            "LV" = igraph::cluster_louvain(gs)$membership,
            "LP" = igraph::label.propagation.community(gs)$membership,
            "EV" = igraph::cluster_leading_eigen(gs)$membership,
            "EB" = igraph::cluster_edge_betweenness(gs)$membership,
            stop("Unsupported method.")
        )
        
        # Remap using node names (always correct if V(g)$name is defined)
        membership <- memb_perm[match(igraph::V(g)$name, igraph::V(gs)$name)]
        
        matched <- FALSE
        if (ncol(M) > 0) {
            for (i in seq_len(ncol(M))) {
                if (eq_1(membership, M[, i])) {
                    counts[i] <- counts[i] + 1
                    matched <- TRUE
                    break
                }
            }
        }
        if (!matched) {
            M <- cbind(M, membership)
            colnames(M)[ncol(M)] <- sprintf("s%02d", ncol(M))
            counts <- c(counts, 1)
        }
        
        names(counts) <- colnames(M)
        counts_by_trial[[t]] <- counts
        
        # CI width
        ci <- beta_ci(counts)
        max_ci_width[t] <- if (length(counts) == 1L) ci$upper - ci$lower else max(ci$upper - ci$lower)
        if (verbose) message(sprintf("[t=%d] ns=%d, max_CI_width=%.4f", t, length(counts), max_ci_width[t]))
        
        if (max_ci_width[t] <= delta) {
            stop_trial <- t
            stop_reason <- "precision"
            break
        }
    }
    
    if (is.na(stop_trial)) {
        stop_trial <- n_trials
        stop_reason <- "budget_exhausted"
    }
    
     # ---- Final posterior summary ----
    ci <- beta_ci(counts)
    probs_df <- data.frame(
        id     = names(counts),
        count  = as.integer(counts),
        phat   = ci$phat,
        pmin   = ci$lower,
        pmax   = ci$upper,
        row.names = NULL
    )
    ord <- order(probs_df$phat, decreasing = TRUE)
    probs_df <- probs_df[ord, , drop = FALSE]
    M <- M[, match(probs_df$id, colnames(M)), drop = FALSE]
    
    # ---- Build trial-wise probability logs: phat / pmin / pmax ----
    all_ids <- colnames(M)
    n_sol   <- length(all_ids)
    
    phat_mat <- matrix(NA_real_, nrow = stop_trial, ncol = n_sol,
                       dimnames = list(trial = seq_len(stop_trial), solution = all_ids))
    pmin_mat <- phat_mat
    pmax_mat <- phat_mat
    
    for (t in seq_len(stop_trial)) {
        ct <- counts_by_trial[[t]]
        if (length(ct)) {
            # compute CI on the solutions known at trial t
            ci_t <- beta_ci(ct)
            # name vectors for alignment
            ph <- as.numeric(ci_t$phat); names(ph) <- names(ct)
            lo <- as.numeric(ci_t$lower); names(lo) <- names(ct)
            up <- as.numeric(ci_t$upper); names(up) <- names(ct)
            
            # write into matrices (only for solutions existing at trial t)
            cols <- intersect(names(ct), all_ids)
            phat_mat[t, cols] <- ph[cols]
            pmin_mat[t, cols] <- lo[cols]
            pmax_mat[t, cols] <- up[cols]
        }
    }
    
    # Also provide long-format log
    prob_long <- do.call(
        rbind,
        lapply(seq_len(stop_trial), function(t) {
            data.frame(
                trial    = t,
                solution = all_ids,
                phat     = phat_mat[t, all_ids],
                pmin     = pmin_mat[t, all_ids],
                pmax     = pmax_mat[t, all_ids],
                row.names = NULL
            )
        })
    )
    
    # ---- Return ----
    list(
        partitions    = M,
        probabilities = probs_df,
        log = list(
            stop_trial    = stop_trial,
            stop_reason   = stop_reason,
            delta         = delta,
            confidence    = confidence,
            gamma0        = gamma0,
            max_ci_width  = head(max_ci_width, stop_trial),
            # counts (kept for reference)
            counts_matrix = {
                counts_mat <- matrix(0L, nrow = stop_trial, ncol = length(all_ids),
                                     dimnames = list(trial = seq_len(stop_trial), solution = all_ids))
                for (t in seq_len(stop_trial)) {
                    ct <- counts_by_trial[[t]]
                    if (length(ct)) counts_mat[t, names(ct)] <- ct
                }
                counts_mat
            },
            # NEW: probabilities per trial and solution
            phat_matrix = phat_mat,
            pmin_matrix = pmin_mat,
            pmax_matrix = pmax_mat,
            prob_long   = prob_long
        )
    )
    
}


#' Classify Solution Space Taxonomy
#'
#' Given the result of \code{solutions_space_DM()}, this function assigns a
#' high-level label to the solution space, based on the number and posterior
#' credibility of distinct partitions found.
#'
#' @param ssp A list output from \code{solutions_space_DM()}. It must contain 
#'   at least the elements \code{partitions}, \code{probabilities}, and \code{log}.
#'
#' @return A single character string indicating the type of solution space. One of:
#' \describe{
#'   \item{\code{"Empty"}}{No valid partitions found (only degenerate ones).}
#'   \item{\code{"Single"}}{Only one solution is credible and clearly identified.}
#'   \item{\code{"Dominant"}}{One solution has a strictly dominant posterior.}
#'   \item{\code{"Multiple"}}{Multiple solutions with no clear dominance.}
#'   \item{\code{"Sparse"}}{Too many solutions relative to trials (unstable space).}
#'   \item{\code{"... uncertain..."}}{Atypical configuration not falling into other types.}
#' }
#'
#' @details
#' The function checks the number of credible solutions, posterior intervals, and convergence 
#' log data to categorize the structure of the solution space. Degenerate solutions 
#' (e.g., all nodes in one group or all in separate groups) are excluded from classification.
#'
#' @seealso \code{\link{solutions_space_DM}} for generating the input object.
#'
#' @examples
#' \dontrun{
#' g <- igraph::make_ring(10)
#' ssp <- solutions_space_DM(g, n_trials = 100, met = "IM", seed = 123)
#' solution_space_type(ssp)
#' }
#'
#' @export

# ============================================================
# solution_space_type(): classify SSP taxonomy
# ============================================================
solution_space_type <- function(ssp) {
    if (!is.list(ssp) || is.null(ssp$probabilities) || is.null(ssp$partitions) || is.null(ssp$log))
        stop("`ssp` must be the result of solutions_space().")
    
    probs <- ssp$probabilities
    M     <- ssp$partitions
    lg    <- ssp$log
    
    ns <- ncol(M)
    if (ns == 1L) return("Single")
    if (ns == 0L) return("Empty")
    
    # Degenerate partitions: all same label or all unique
    valid <- rep(TRUE, ns)
    for (i in seq_len(ns)) {
        k <- length(unique(M[, i]))
        if (k == 1L || k == nrow(M)) valid[i] <- FALSE
    }
    if (!any(valid)) return("Empty")
    
    pl <- probs$plower
    pu <- probs$pupper
    ph <- probs$phat
    top <- which.max(ph)
    max_width <- tail(lg$max_ci_width, 1)
    
    if (ns == 1L && is.finite(max_width) && max_width <= lg$delta) return("Single")
    if (length(ph) >= 2L && pl[top] > max(pu[-top])) return("Dominant")
    if (lg$stop_reason == "budget_exhausted" && ns / lg$stop_trial >= 0.20) return("Sparse")
    if (max(pl) < 0.5) return("Multiple")
    return("... uncertain...")
}



#' Build the pairwise agreement matrix Γ via explicit triple loops
#'
#' @description
#' Direct, publication-style implementation of
#'   γ_uv = Σ_i p_hat[i] * 1{u and v co-occur in the same community in partition i}.
#' This version is intentionally **SLOW but CLEAR** (O(ns * nv^2)), ideal for teaching,
#' reproducibility, and appendices. For large nv use a block/sparse approach.
#'
#' @param partitions Integer/character/factor **matrix** of size nv × ns
#'   (rows = nodes, columns = partitions), OR a **vector** of length nv for the
#'   single-partition (degenerate) case. `partitions[v, i]` is the label of node v
#'   in partition i. Labels may be integers, characters, or factors.
#' @param probs Optional numeric vector of length ns with nonnegative weights p̂[i]
#'   that sum (approximately) to 1. If `partitions` encodes only one partition
#'   (vector or one-column matrix) and `probs` is missing, it defaults to 1.
#'   If provided and not summing to 1 within tolerance, it is softly normalized
#'   with a warning.
#' @param tol Numeric tolerance for the sum-to-1 check on `probs`. Default 1e-9.
#'
#' @returns
#' An nv × nv numeric matrix Γ with entries in [0, 1] and unit diagonal.
#'
#' @details
#' • Complexity: O(ns * nv^2). **Slow but crystal-clear** (explicit indicator in
#'   the innermost loop).  
#' • `partitions` must refer to a fixed node ordering across columns.  
#' • Degenerate (single-partition) case: Γ[u, v] = 1 if u and v share the same
#'   community in that partition, else 0; diagonal forced to 1.
#'
#' @examples
#' # --- Multi-partition (nv = 5, ns = 2) as an nv × ns matrix ---
#' P_mat <- cbind(
#'   c(1,1,2,2,2),  # partition 1
#'   c(1,2,2,2,2)   # partition 2
#' )
#' Gamma1 <- build_gamma_loops_matrix(P_mat, probs = c(0.6, 0.4))
#'
#' # --- Single-partition as a vector (degenerate case) ---
#' P_vec <- c(1,1,2,2,2)  # length nv
#' Gamma2 <- build_gamma_loops_matrix(P_vec)  # probs defaults to 1
#'
#' # --- Single-partition as a one-column matrix (nv × 1) ---
#' P_one <- matrix(c(1,1,2,2,2), ncol = 1)
#' Gamma3 <- build_gamma_loops_matrix(P_one, probs = 1)

build_gamma_matrix <- function(partitions, probs = NULL, tol = 1e-9) {
    # ---- Accept matrix (nv × ns) OR vector (length nv) --------------------------
    if (is.null(partitions) || length(partitions) == 0L) {
        stop("`partitions` must be a non-empty matrix (nv × ns) or a vector (length nv).")
    }
    
    # Convert vector -> nv × 1 matrix for unified handling
    if (!is.matrix(partitions)) {
        # Expect a vector encoding a single partition
        if (!is.atomic(partitions) && !is.factor(partitions)) {
            stop("`partitions` must be a matrix or an atomic/factor vector.")
        }
        if (is.null(dim(partitions))) {
            partitions <- matrix(partitions, ncol = 1L)
        } else {
            stop("Unexpected `partitions` structure; provide a vector or an nv × ns matrix.")
        }
    }
    
    # Now guaranteed matrix
    nv <- nrow(partitions)
    ns <- ncol(partitions)
    if (nv < 1L || ns < 1L) {
        stop("`partitions` must have at least one row (node) and one column (partition).")
    }
    if (anyNA(partitions)) {
        stop("`partitions` contains NA values; please impute or remove them.")
    }
    
    # ---- Handle probs (including degenerate single-partition defaults) ----------
    if (is.null(probs)) {
        if (ns == 1L) {
            probs <- 1.0
        } else {
            stop("`probs` must be provided for multiple partitions (length ns).")
        }
    }
    if (!is.numeric(probs) || length(probs) != ns) {
        stop("`probs` must be a numeric vector of length equal to ncol(partitions).")
    }
    if (any(!is.finite(probs)) || any(probs < -1e-15)) {
        stop("`probs` must be finite and nonnegative (within numerical tolerance).")
    }
    s <- sum(probs)
    if (!is.finite(s) || s <= 0) {
        stop("Sum of `probs` must be positive and finite.")
    }
    if (abs(s - 1) > max(tol, 1e-12)) {
        warning("`probs` do not sum to 1; normalizing to sum exactly 1.")
        probs <- probs / s
    }
    
    # ---- Initialize Γ -----------------------------------------------------------
    Gamma <- matrix(0.0, nrow = nv, ncol = nv)
    
    # ---- Plain triple loop (slow but clear) -------------------------------------
    # For each partition i, for each pair (u, v), add weight if co-assigned.
    for (i in seq_len(ns)) {
        w <- probs[i]
        memb <- partitions[, i]  # labels for partition i (length nv)
        for (u in seq_len(nv)) {
            for (v in seq_len(nv)) {
                if (memb[u] == memb[v]) {
                    Gamma[u, v] <- Gamma[u, v] + w
                }
            }
        }
    }
    
    # ---- Numerical safety and invariants ----------------------------------------
    # Clip tiny excursions and enforce unit diagonal (each node co-occurs with itself).
    Gamma[Gamma < 0] <- 0
    Gamma[Gamma > 1] <- 1
    diag(Gamma) <- 1
    
    return(Gamma)
}

#' Calculate Normalized Co-occurrence Matrix
#'
#' The `co_occurrence` function computes a normalized co-occurrence matrix based on the 
#' solution space provided by the input object `ssp`. 
#' The output matrix indicates the extent to which pairs of nodes appear together in the same community 
#' across different trials, with weights normalized by the median values of the solutions.
#' 
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
    # Weighted co-occurrence (probability two nodes co-cluster across solutions)
    # - Uses ssp$partitions: rows = nodes, cols = solutions (partitions)
    # - Uses ssp$probabilities$phat as weights; defaults to equal weights
    # - Returns an n_nodes x n_nodes symmetric matrix in [0, 1] with diag = 1
    
    # --- inputs ---
    M <- ssp$partitions
    if (is.data.frame(M)) M <- as.matrix(M)
    if (!is.numeric(M)) stop("ssp$partitions must be a numeric matrix of community labels.")
    n_nodes <- nrow(M); n_solutions <- ncol(M)
    if (is.null(n_nodes) || is.null(n_solutions) || n_nodes < 1 || n_solutions < 1)
        stop("ssp$partitions must have at least one row (node) and one column (solution).")
    
    # --- weights ---
    w <- ssp$probabilities$phat
    if (is.null(w) || length(w) != n_solutions || any(!is.finite(w))) {
        w <- rep(1 / n_solutions, n_solutions)    # equal weights fallback
    } else {
        s <- sum(w)
        if (!isTRUE(all.equal(s, 1))) w <- w / s  # normalize to sum to 1
    }
    
    # --- initialize output ---
    if (is.null(rownames(M))) rownames(M) <- as.character(seq_len(n_nodes))
    D <- matrix(0, nrow = n_nodes, ncol = n_nodes,
                dimnames = list(rownames(M), rownames(M)))
    rownames(D) <- rownames(M) 
    colnames(D) <- rownames(M)
    
    # --- accumulate co-occurrences ---
    for (t in seq_len(n_solutions)) {
        labels_t <- M[, t]
        if (anyNA(labels_t)) next  # skip NA labels in this solution
        
        # Use unique labels (includes 0 if present) — avoids the 1:max(...) pitfall
        for (k in sort(unique(labels_t))) {
            comm_members <- which(labels_t == k)
            nc <- length(comm_members)
            if (nc >= 2) {
                # pairwise updates within the community
                for (i in 1:(nc - 1)) {
                    for (j in (i + 1):nc) {
                        ii <- comm_members[i]; jj <- comm_members[j]
                        D[ii, jj] <- D[ii, jj] + w[t]
                        D[jj, ii] <- D[ii, jj]           # keep symmetric
                    }
                }
            }
            # (nc == 1): no pairs to update
        }
    }
    
    # Every node co-occurs with itself with probability 1
    diag(D) <- 1
    
    return(D)
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


