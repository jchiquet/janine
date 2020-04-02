#' Estimation of the latent block organisation of the network
#'
#' The underlying network is assumed to be drawn from a Stochastic Bloc Model.
#' This function uses a variational EM algorithm implemented in the package blockmodels to estimate such a model.
#'
#' @param adj_matrix a symmetric weighted adjacency matrix
#' @param n_blocks  integer, the number of blocks
#' @import blockmodels
#' @export
estimate_block <- function(adj_matrix, n_blocks){

  SBM_fits <- BM_bernoulli("SBM", adj_matrix, verbosity = 0, plotting="", explore_min = n_blocks)
  SBM_fits$estimate()

  list(
    blockProb   = SBM_fits$memberships[[n_blocks]]$Z,
    connectProb = SBM_fits$prediction(n_blocks),
    membership  = apply(SBM_fits$memberships[[n_blocks]]$Z, 1, which.max),
    parameters  = SBM_fits$model_parameters[[n_blocks]]
  )
}


# entropy <- function(distr) {
#   ## handle x log(x) = 0 when x = 0
#   -sum(distr * log(distr + 1*(distr == 0)))
# }
#
# estimate_structure <- function(adj_matrix, n_blocks, init_clustering = kmeans(adj_matrix, centers = n_blocks, nstart = 10)$cl, threshold = 1e-5, max_iterates = 50){
#
#   ## Initialization
#   p <- ncol(adj_matrix)
#   J <- vector("numeric", max_iterates)
#   zero  <- .Machine$double.threshold
#
#   ### variational lower bound
#   get_J <- function(theta, tau){
#     J <- sum( tau %*% log(theta$alpha) ) + entropy(tau) -
#        sum( tau %*%  log(2*theta$Lambda) %*% t(tau) + abs(adj_matrix) * tau %*% (1/theta$Lambda) %*% t(tau) )
#     J
#   }
#
#   ### M step: update theta (pi and alpha)
#   M_step <- function(tau){
#     Lambda <- (t(tau) %*% abs(adj_matrix) %*% tau) / (t(tau) %*% (1 - diag(1, p, p)) %*% tau)
#     alpha <- colMeans(tau)
#     alpha[alpha < zero] <- zero
#     list(Lambda = Lambda, alpha = alpha)
#   }
#
#   ### E step: update the clustering parameters (tau)
#   E_step <- function(theta, tau){
#     alpha  <- theta$alpha
#     Lambda <- theta$Lambda
#     tau <- matrix(log(alpha), p, n_blocks, byrow = TRUE) -
#               (1 - diag(1, p, p)) %*% tau %*% log(2 * theta$Lambda) - abs(adj_matrix) %*% tau %*% (1/theta$Lambda)
#     tau <- exp(tau - apply(tau, 1, max))
#     tau <- tau/rowSums(tau)
#     tau
#   }
#
#   VEM <- function(tau) {
#     cond <- FALSE; iter <- 0
#     while (!cond){
#       iter <- iter +1
#       theta <- M_step(tau)
#       # E step
#       tau <- E_step(theta, tau) # M step
#
#       J[iter] <- get_J(theta, tau) # assess convergence
#       if (iter > 1)
#         cond <- (iter > max_iterates) | (abs((J[iter] - J[iter-1])) < threshold)
#     }
#     return(list(theta = theta, tau = tau, J = J[1:iter]))
#   }
#
#   tau <- matrix(0,p,n_blocks); tau[cbind(1:p, init_clustering)] <- 1
#   best <- VEM(tau)
#
#   vBIC <- best$J[length(best$J)] - .5*(n_blocks*(n_blocks+1)/2)*log(p*(p - 1)/2) + (n_blocks - 1)*log(p)
#   vICL <- vBIC - entropy(best$tau)
#
#   list(theta = best$theta,
#        tau = best$tau,
#        membership = apply(best$tau, 1, which.max),
#        J = best$J, vICL = vICL, vBIC= vBIC)
# }
