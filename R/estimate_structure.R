entropy <- function(distr) {
 ## handle x log(x) = 0 when x = 0
 -sum(distr * log(distr + 1*(distr == 0)))
}

#' Estimation of the latent structure of the network
#'
#' The underlying network is assumed to be drawn from a Stochastic Bloc Model
#' with Laplace emission law for the edges. This function use a variational EM algorithm
#' to estimate such a model.
#'
#' @param adj_matrix a symmetric weighted adjacency matrix
#' @param n_blocks  integer, the number of blocks
#' @param init_clustering vector, an initial vector of memberships
#' @param threshold double, a convergence threshold
#' @param max_iterates integer,  the maximal number of itarates in the VEM algorithm
#' @importFrom stats kmeans rmultinom
#' @examples
#' Lambda <- matrix(c(1, 2, 3, 2, 1, 4, 3, 4, 1), 3, 3)
#' d <- 200
#' alpha <- c(1/2, 1/4, 1/4)
#' mySBM <- rSBM_Laplace(d, Lambda, alpha)
#' cl <- igraph::V(mySBM)$membership
#' adjMat <- mySBM %>% igraph::as_adj(attr = "weight")
#' vBlocks <- 1:5
#' out <- lapply(vBlocks, function(k) estimate_structure(adjMat, k))
#' @export
estimate_structure <- function(adj_matrix, n_blocks, init_clustering = kmeans(adj_matrix, centers = n_blocks, nstart = 10)$cl, threshold = 1e-5, max_iterates = 50){

  ## Initialization
  p <- ncol(adj_matrix)
  J <- vector("numeric", max_iterates)
  zero  <- .Machine$double.threshold

  ### variational lower bound
  get_J <- function(theta, tau){
    J <- sum( tau %*% log(theta$alpha) ) + entropy(tau) -
       sum( tau %*%  log(2*theta$Lambda) %*% t(tau) + abs(adj_matrix) * tau %*% (1/theta$Lambda) %*% t(tau) )
    J
  }

  ### M step: update theta (pi and alpha)
  M_step <- function(tau){
    Lambda <- (t(tau) %*% abs(adj_matrix) %*% tau) / (t(tau) %*% (1 - diag(1, p, p)) %*% tau)
    alpha <- colMeans(tau)
    alpha[alpha < zero] <- zero
    list(Lambda = Lambda, alpha = alpha)
  }

  ### E step: update the clustering parameters (tau)
  E_step <- function(theta, tau){
    alpha  <- theta$alpha
    Lambda <- theta$Lambda
    tau <- matrix(log(alpha), p, n_blocks, byrow = TRUE) -
              (1 - diag(1, p, p)) %*% tau %*% log(2 * theta$Lambda) - abs(adj_matrix) %*% tau %*% (1/theta$Lambda)
    tau <- exp(tau - apply(tau, 1, max))
    tau <- tau/rowSums(tau)
    tau
  }

  VEM <- function(tau) {
    cond <- FALSE; iter <- 0
    while (!cond){
      iter <- iter +1
      theta <- M_step(tau)
      # E step
      tau <- E_step(theta, tau) # M step

      J[iter] <- get_J(theta, tau) # assess convergence
      if (iter > 1)
        cond <- (iter > max_iterates) | (abs((J[iter] - J[iter-1])) < threshold)
    }
    return(list(theta = theta, tau = tau, J = J[1:iter]))
  }

  tau <- matrix(0,p,n_blocks); tau[cbind(1:p, init_clustering)] <- 1
  best <- VEM(tau)

  vBIC <- best$J[length(best$J)] - .5*(n_blocks*(n_blocks+1)/2)*log(p*(p - 1)/2) + (n_blocks - 1)*log(p)
  vICL <- vBIC - entropy(best$tau)

  return(
    list(theta = best$theta,
         tau = best$tau,
         membership = apply(best$tau, 1, which.max),
         J = best$J, vICL = vICL, vBIC= vBIC)
    )
}
