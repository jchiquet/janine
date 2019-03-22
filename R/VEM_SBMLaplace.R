entropy <- function(distr) {
 ## handle x log(x) = 0 when x = 0
 -sum(distr * log(distr + 1*(distr == 0)))
}

#' Estimation of a SBM Laplace
#'
#' This is the description of my function.
#' Estimation of Stochastic Bloc Model under a Laplace emission law
#'
#' @param Omega a matrix: the data
#' @param Q an integer: the number of clusters
#' @param cl.init a initial clustering
#' @param eps double convergence threshold
#' @param maxIter integer the maximal number of itarates in the VEM algorithm
#' @importFrom stats kmeans rmultinom
#' @examples
#' Lambda <- matrix(c(1, 2, 3, 2, 1, 4, 3, 4, 1), 3, 3)
#' p <- 200
#' alpha <- c(1/2, 1/4, 1/4)
#' mySBM <- rSBMLaplace(p, Lambda, alpha)
#' cl <- igraph::V(mySBM)$membership
#' Omega <- mySBM %>% igraph::as_adj(attr = "weight")
#' vBlocks <- 1:5
#' out <- lapply(vBlocks, function(q) VEM_SBM_laplace(Omega, q))
#' @export
VEM_SBM_laplace <- function(Omega, Q, cl.init = kmeans(Omega, centers = Q, nstart = 10)$cl, eps = 1e-5, maxIter = 50){

  ## Initialization
  p <- ncol(Omega)
  J <- vector("numeric", maxIter)
  zero  <- .Machine$double.eps

  ### variational lower bound
  get_J <- function(theta, tau){
    J <- sum( tau %*% log(theta$alpha) ) + entropy(tau) -
       sum( tau %*%  log(2*theta$Lambda) %*% t(tau) + abs(Omega) * tau %*% (1/theta$Lambda) %*% t(tau) )
    J
  }

  ### M step: update theta (pi and alpha)
  M_step <- function(tau){
    Lambda <- (t(tau) %*% abs(Omega) %*% tau) / (t(tau) %*% (1 - diag(1, p, p)) %*% tau)
    alpha <- colMeans(tau)
    alpha[alpha < zero] <- zero
    list(Lambda = Lambda, alpha = alpha)
  }

  ### E step: update the clustering parameters (tau)
  E_step <- function(theta, tau){
    alpha  <- theta$alpha
    Lambda <- theta$Lambda
    tau <- matrix(log(alpha), p, Q, byrow = TRUE) -
              (1 - diag(1, p, p)) %*% tau %*% log(2 * theta$Lambda) - abs(Omega) %*% tau %*% (1/theta$Lambda)
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
        cond <- (iter > maxIter) | (abs((J[iter] - J[iter-1])) < eps)
    }
    return(list(theta = theta, tau = tau, J = J[1:iter]))
  }

  tau <- matrix(0,p,Q); tau[cbind(1:p, cl.init)] <- 1
  best <- VEM(tau)

  vBIC <- best$J[length(best$J)] - .5*(Q*(Q+1)/2)*log(p*(p - 1)/2) + (Q - 1)*log(p)
  vICL <- vBIC - entropy(best$tau)
  return(list(theta = best$theta,
              tau = best$tau, membership = apply(best$tau, 1, which.max),
              J = best$J, vICL = vICL, vBIC= vBIC))
}
