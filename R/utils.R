#' rSBM_Laplace
#'
#' Generate a symmetric matrix under the stochastic block model
#' with a Laplace distribution as the emission law of the edges.
#'
#' @param p number of nodes
#' @param Lambda a QxQ matrix of connectivity parameters
#' @param alpha a size-Q vector of mixture
#' @return an igraph object: undirected weigthed graph with a attribute "membership" for the vertices
#' @import Matrix igraph
#' @importFrom rmutil rlaplace
#' @export
rSBM_Laplace <- function(p, Lambda, alpha = c(1)) {
  Q <- length(alpha)
  Z <- t(rmultinom(p, 1, prob = alpha))
  A <- matrix(0, p, p)
  subset <- upper.tri(A)
  lambda <- (Z %*% Lambda %*% t(Z))[subset]
  x <- rmutil::rlaplace(n=p*(p-1)/2, m = 0, s = lambda)
  A[subset] <- x
  A <- A + t(A)
  mySBM <- graph_from_adjacency_matrix(A, weighted = TRUE, mode = "undirected")
  vertex_attr(mySBM, "membership") <- Z %*% 1:Q
  mySBM
}

#' SBM-structured multivariate Gaussian observation
#'
#' Draw observations from a multivariate Gaussian distribution with SBM-structured precision matrix
#'
#' @param n integer, the sample size
#' @param sbm an SBM, as outputed by the \code{\link{rSBM_Laplace}} function
#'
#' @importFrom igraph as_adj
#' @importFrom mvtnorm rmvnorm
#' @examples
#' n <- 200; d <- 100; K <- 3
#' Lambda <- matrix(c(1, 2, 3, 2, 1, 4, 3, 4, 1), 3, 3)
#' alpha <- c(1/2, 1/4, 1/4)
#' K <- length(alpha)
#' sbm   <- rSBM_Laplace(d, Lambda, alpha)
#' X <- rmvnorm_SBM(n, sbm)
#' @export
rmvnorm_SBM <- function(n, sbm) {
  Omega <- as_adj(sbm, attr = "weight")
  diag(Omega) <- colSums(abs(Omega))
  Sigma <- as.matrix(chol2inv(chol(Omega)))
  res <- mvtnorm::rmvnorm(n, sigma = Sigma)
  res
}

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
