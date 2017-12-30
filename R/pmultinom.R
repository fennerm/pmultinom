#' Megadaph mitochondrial variant calling R library

#' Expected Maximum Order Statistic Test
#'
#' Determine the probability of drawing a value as large as the observed maximum
#' of a numeric vector by sampling from a uniform distribution. I.e test whether
#' the maximum value is unexpectedly large.
#' @param method One of "sim", "corrado", "levin" TODO: explain
#' @export
multinomial_range_test <- function() {
    print()
}

#' Simulation based expected maximum order statistic test
#'
#' Contingency tables are repeatedly simulated with the same row and column
#' marginals as the input table. The proportion of these tables with a maximum
#' at least as extreme as the input table are returned.
#' @param ctab Matrix; the input contingency table (k x 2)
#' @param reps int; The number of repetitions. Will be rounded to the nearest
#'        10000
#' @return A p-value
sim_multinomial_range_test <- function(ctab, reps=10000) {
    # Calculate row and column marginals
    row_marginal <- rowSums(ctab)
    col_marginal <- colSums(ctab)

    if (col_marginal[1] == 0) {
        1
    } else {
        # Number of samples
        nsample <- nrow(ctab)

        # Observed maximum as a frequency
        obs_max <- sort(ctab[, 1] / row_marginal)[nsample]

        null_dist <- sim_null_dist(reps, row_marginal, col_marginal)

        # Proportion of simulated data at least as extreme as obs_max
        sim_more_extreme <- length(which(null_dist >= obs_max))

        # Calculate p-value. One is added so that p-values are never 0.
        p <- (sim_more_extreme + 1) / (reps + 1)

        p
    }
}

#' Permute contingency tables with given row and column marginals and return
#' a vector of observed maximum frequencies
sim_null_dist <- function(n, row_marginal, col_marginal) {
    if (n < 10000) {
        nsets <- 1
    } else {
        nsets <- round(n / 10000)
    }
    # Repeatedly simulate pools of 10000 tables to reduce memory overhead
    unlist(replicate(nsets, permute_max(10000, row_marginal, col_marginal),
               simplify = FALSE))
}

#' Permute contingency matrices with fixed row and column marginals and
#' calculate the maximal allele frequency.
#' @importFrom matrixStats colMaxs
#' @importFrom stats rmultinom
permute_max <- function(reps, row_marginal, col_marginal) {
    prob <- row_marginal / sum(row_marginal)
    # Permute column 1 counts with probability proportional to row
    # marginals.
    sim_mut_counts <- rmultinom(reps, size = col_marginal[1], prob = prob)
    sim_wt_counts <- rmultinom(reps, size = col_marginal[2], prob = prob)
    # Convert to proportions
    sim_prop <- sim_mut_counts / (sim_mut_counts + sim_wt_counts)
    rm(sim_mut_counts, sim_wt_counts)
    max_prop <- colMaxs(sim_prop)
    rm(sim_prop)
    max_prop
}


#' Calculate multinomial range probability using the sum of poisson random
#' variables
#' @export
levin_test <- function(ctab) {
    row_marginals <- rowSums(ctab)
    col_marginals <- colSums(ctab)
    nsamples <- nrow(ctab)

    # Observed mutant allele frequency
    obs <- sort(ctab[, 1] / row_marginals)
    prob <- row_marginals / sum(row_marginals)
    obs_max <- max(obs)
    cutoff <- ceiling(obs_max * row_marginals)
    nmuts <- col_marginals[1]
    if (nmuts == 0) {
        1
    } else {
        plevin(cutoff, nmuts, prob, lower.tail=FALSE)
    }
}

#' @importFrom truncdist dtrunc
dtpois <- function(nballs, prob, cutoff, lower, upper) {
    d <- dtrunc(lower:upper, lambda = nballs*prob, spec = "pois",
                a = lower, b = upper)
    d[(cutoff + 1):(upper+1)] <- 0
    d
}

#' Calculate multinomial range probability using the sum of poisson random
#' variables
#'
#' Can be used to calculate the probabilities of the various outcomes of
#' throwing n balls in k bins.
#'
#' Algorithm described in:
#' Levin,  B.:  A  representation  for  multinomial  cumulative  distribution
#' functions. Ann. Stat. 9, 1123–1126 (1981)
#' @param q Vector of quantiles; If q has length 1 then the same quantile is
#'          used for each bin. If q is length k, each bin is assigned its own
#'          quantile
#' @param n Number of balls to be placed in k baskets
#' @param prob A vector of length k, describing the probability of a ball being
#'             placed into bin b_k
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x],
#'                   otherwise, P[X > x].
#' @return The probability of
plevin <- function(q, n, prob, lower.tail = TRUE) {
    nbins <- length(prob)
    if (length(q) == 1) {
        q <- rep(q, nbins)
    }

    conv <- dtpois(n, prob[1], q[1], 0, q[1])
    dens <- conv
    if (nbins >= 2) {
        for (i in 2:nbins) {
            dens <- c(dtpois(n, prob[i], q[i], 0, q[i]),
                      rep(0, length(conv) - length(dens)))
            conv <- convolve(conv, rev(dens), type = "open")
        }
    }
    conv <- conv / sum(conv)
    conv_p <- conv[n + 1]
    if (conv_p < 0) {
        conv_p <- 0
    }

    1 - (sqrt(2 * pi * n)) *
        exp(sum(ppois(q - 1, n * prob, log.p = TRUE)) +
        log(conv_p))
}

#' Multinomial range test using stochastic matrices
#'
#' Algorithm described in:
#' Corrado, Charles J. “The Exact Distribution of the Maximum, Minimum and the
#' Range of Multinomial/Dirichlet and Multivariate Hypergeometric Frequencies.”
#' Statistics and Computing 21, no. 3 (July 1, 2011): 349–59.
#' https://doi.org/10.1007/s11222-010-9174-3.
#' @param ctab Matrix; A contingency table
#' @return A p-value
corrado_test <- function(ctab) {
    row_marginals <- rowSums(ctab)
    col_marginals <- colSums(ctab)
    nsamples <- nrow(ctab)
    obs <- sort(ctab[, 1] / row_marginals)
    prob <- row_marginals / sum(row_marginals)
    obs_max <- max(obs)
    cutoff <- ceiling(obs_max * row_marginals)

    # The p-value for the multinomial range is guaranteed to be smaller than or
    # equal to the p-value from the sum of independent binomial probabilities.
    # Therefore if the sum of independent probabilities is sufficiently small
    # (<2.2e-16) we don't need to bother computing the multinomial probability
    p_ind <- sum(pbinom(cutoff - 1, col_marginals[1], prob, lower.tail = F))

    # If there are no balls, p is always 1
    if (col_marginals[1] == 0) {
        p <- 1
    } else if (sum(sort(cutoff)[1:2]) > col_marginals[1]) {
        if (p_ind <= 0) {
            p <- .Machine$double.xmin
        } else {
            p <- p_ind
        }
    } else if (p_ind < 2.2e-16) {
        p <- 2.2e-16
    } else {
        choose_mat <- build_choose_matrix(col_marginals[1])
        cum_product <- build_stochastic_matrix(1, prob, choose_mat)
        cum_product <- cull(cum_product, cutoff[1])
        for (i in 2:nsamples) {
            M <- build_stochastic_matrix(i, prob, choose_mat)
            culled <- cull(M, cutoff[i])
            rm(M)
            if (i == nsamples) {
                culled <- matrix(culled[, ncol(culled)])
            }
            cum_product <- logmatrix_mult(cum_product, culled)
            rm(culled)
        }
        rm(choose_mat)
        p <- -expm1(cum_product)
    }

    if (p <= 0) {
        p <- 2.2e-16
    } else if (p > 1) {
        p <- 1
    }

    p
}

#' Compute the log of the sum of exponentials of input elements.
#' @param lx numeric vector
#' @return numeric
log_sum_exp <- function(lx) {
    ## We define our own function rather than using matrixStats::logSumExp
    ## because it was causing underflow errors
    max_idx <- which.max(lx)
    logsum <- log1p(sum(exp(lx[-max_idx] - lx[max_idx]))) + lx[max_idx]
    logsum[is.nan(logsum)] <- -Inf
    logsum
}

#' Normalize a log probability vector or matrix
normalize <- function(log_prob) {
    if (is.matrix(log_prob)) {
        sweep(log_prob, 1, apply(log_prob, 1, log_sum_exp))
    } else {
        log_prob - log_sum_exp(log_prob)
    }
}

corrado_prob <- function(j, i, nballs, log_pik, log_1_pik) {
    index_diff <- j-i

    if (index_diff >= 0) {
        if ((nballs-j == 0) && (log_pik == 0)) {
            p <- 0
        } else {
            p <- lchoose(nballs - i, index_diff) + index_diff * log_pik +
                ((nballs-j) * log_1_pik)
        }
    } else {
        p <- -Inf
    }

    p
}

#' @importFrom matrixStats colLogSumExps
build_choose_matrix <- function(nballs) {
    dim <- nballs+1
    choose_mat <- matrix(, ncol=dim, nrow=dim)
    choose_mat[dim, ] <- c(rep(-Inf, nballs), 0)
    for (i in nballs:1) {
        shifted <- c(choose_mat[i+1, 2:dim], 0)
        row_mat <- matrix(c(shifted[1:nballs],choose_mat[i+1, 1:nballs]),
                          ncol=nballs, byrow=TRUE)
        new_row <- colLogSumExps(row_mat)
        choose_mat[i, ] <- c(new_row, 0)
    }
    choose_mat
}

probability_proportion <- function(prob, k) {
    n <- length(prob)
    prob[k] / (sum(prob[k:n]))
}

build_stochastic_matrix <- function(k, prob, choose_mat) {
    nbins <- length(prob)
    dim <- ncol(choose_mat)
    nballs <- dim - 1
    pik <- probability_proportion(prob, k)
    log_pik <- log(pik)
    log_pik_inv <- log(1 - pik)
    if (k == 1) {
        stoch_mat <- choose_mat[1,] + (0:nballs) * log_pik +
            (nballs:0) * log_pik_inv
        stoch_mat <- matrix(stoch_mat, ncol=dim)
        stoch_mat <- normalize(stoch_mat)
    } else if (k == nbins) {
        stoch_mat <- rep(c(rep(-Inf, nballs), 0), dim)
        stoch_mat <- matrix(stoch_mat, ncol = dim, byrow = TRUE)
    } else {
        pik_diff <- (0:nballs) * log_pik
        pik_diff_inv <- (nballs:0) * (log_pik_inv)
        pik_mat <- matrix(, ncol=dim, nrow=dim)
        for (i in 1:dim) {
            pik_mat[i, ] <- pik_diff_inv + pik_diff
            pik_diff <- c(-Inf, pik_diff[1:(dim - 1)])
        }
        stoch_mat <- choose_mat + pik_mat
        stoch_mat <- normalize(stoch_mat)
    }
    stoch_mat
}

build_stochastic_matrices <- function(prob, nballs) {
    choose_mat <- build_choose_matrix(nballs)
    nbins <- length(prob)
    stoch_mats <- lapply(1:nbins, build_stochastic_matrix,
                         prob=prob, choose_mat=choose_mat)
    stoch_mats
}

#' @importFrom matrixStats colLogSumExps
logmatrix_mult <- function(mat1, mat2) {
    m <- nrow(mat1)
    n <- ncol(mat2)

    matrix(colLogSumExps(as.vector(mat1) + mat2))
}

cull <- function(mat, cutoff) {
    mat[col(mat) - row(mat) >= cutoff] <- -Inf
    mat
}
