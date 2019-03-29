#' Adaptive graphical-Lasso
#'
#' Estimation of Stochastic Bloc Model under a Laplace emission law
#'
#' @param S a d x d symmetric empirical covariance matrix
#' @param Lambda a K x K matrix of connectivity
#' @param tau a d x K of posterior probabilities of memberships
#' @importFrom glassoFast glassoFast
#' @export
#' @examples
#' n <- 200; d <- 100; K <- 3
#' Lambda <- matrix(c(1, 2, 3, 2, 1, 4, 3, 4, 1), 3, 3)
#' alpha <- c(1/2, 1/4, 1/4)
#' K <- length(alpha)
#' sbm   <- rSBM_Laplace(d, Lambda, alpha)
#' cl <- igraph::V(sbm)$membership
#' Z <- matrix(0, d, K); Z[cbind(1:d, cl)] <- 1
#' X <- rmvnorm_SBM(n, sbm)
#' out <- estimate_precison(cov(X), 5 * n * Lambda, Z)
#' network <- out$wi; diag(network) <- 0
#' \dontrun{
#' image(Matrix(network[order(cl),order(cl)]))
#' }
estimate_precison <- function(S, Lambda, tau) {

  ## sanity checks
  stopifnot(all.equal(ncol(S), nrow(tau)))
  stopifnot(all.equal(ncol(tau), ncol(Lambda), nrow(Lambda)))
  stopifnot(all(isSymmetric(S), isSymmetric(Lambda)))

  glasso_out <- glassoFast::glassoFast(S, tau %*% (1/Lambda) %*% t(tau))
  glasso_out

}

