#' Estimation of the network structure (edges)
#'
#' The network is estimated with the adaptive Graphical-Lasso, where weigths are inversely
#' proportional to the probability of connection estimated by the SBM
#'
#' @param S a d x d symmetric empirical covariance matrix
#' @param lambda postive real tuning the network sparsity
#' @param W a d x d matrix of weights
#' @importFrom glassoFast glassoFast
#' @export
estimate_network <- function(S, lambda, W = matrix(1, nrow(S), ncol(S))) {

  ## sanity checks
  stopifnot(lambda >= 0)
  stopifnot(all.equal(ncol(S), nrow(S), ncol(W), nrow(W)))
  stopifnot(all(isSymmetric(S), isSymmetric(W)))

  glasso_out <- glassoFast::glassoFast(S, lambda * W)
  Omega <- glasso_out$wi
  support <- (Omega != 0) * 1; diag(support) <- 0

  list(
    Omega = glasso_out$wi,
    Sigma = glasso_out$w,
    support = support
  )
}

