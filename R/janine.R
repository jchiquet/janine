#' Janine: Just Another Network Infernce mEthod
#'
#' Iterate Adpative Graphical-Lasso with binary SBM estimation to recover the adaptive weights
#'
#' @param data a n x d matrix of multivariate Gaussian observation
#' @param penalty postive real tuning the network sparsity
#' @param n_blocks integer for the targetnumber of groups
#' @param epsilon convergence threshold (absolute relative change in objective function), default to 1e-3.
#' @param max_iter integer for the maximal number of iteration. Default to 20
#' @export
janine <- function(data, penalty, n_blocks, epsilon = 1e-4, max_iter = 20, penalize_diagonal = TRUE) {

  n <- nrow(data)
  d <- ncol(X)
  S <- cov(as.matrix(data))
  weights <- matrix(1, d, d)
  if (!penalize_diagonal) diag(weights) <- 0

  cond <- FALSE
  iter <- 1
  objective <- numeric(max_iter)
  objective[iter] <- Inf
  while(!cond) {
    ## M step (network structure estimation)
    net <- estimate_network(S, penalty, weights)
    sparsity <- 1 - sum(net$support) / (d**2)

    ## E step (latent block estimation)
    sbm <- estimate_block(net$support, n_blocks)
    weights <- (1 - sbm$connectProb)/sparsity
    if (!penalize_diagonal) diag(weights) <- 0

    ## Convergence assesment
    logdet <- determinant(net$Omega, logarithm = TRUE)$modulus
    objective[iter + 1] <- - (n/2) * (logdet - sum( diag( S %*% net$Omega )) ) + penalty * sum(abs(weights * net$Omega))
    cond <- abs(objective[iter + 1] - objective[iter])/abs(objective[iter + 1]) < epsilon | iter + 1 > max_iter
    iter <- iter + 1
  }

  list(network = net, membership = sbm, objective = objective[1:(iter - 1)])

}
