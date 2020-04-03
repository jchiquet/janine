#' Janine: Just Another Network Infernce mEthod
#'
#' Iterate Adaptive Graphical-Lasso with binary SBM estimation to recover the adaptive weights
#'
#' @param data a n x d matrix of multivariate Gaussian observation
#' @param n_blocks integer for the targetnumber of groups
#' @param penalties a vector of postive real number in decreasing order tuning the network sparsity. The default (NULL) generates
#' an hopefully appropriate collection of penalties.
#' @param control_penalties a list controling how \code{penalties} is generated, with three entries:
#' a double \code{min_ratio} (default 0.1), a integer \code{length} (default 25) and a logical \code{diagonal} (default FALSE)
#' indicating weither the diaognal should be penalized or not.
#' @param control_optim a list controling how the alternate optimization between adaptive graphical-Lasso and SBM is conducted
#' for each penalty level. Contains three entries: a convergence threshold \code{epsilon} (default to 1e-3),
#' a maximal number of iteration \code{max_iter} (default 20) and verbosity level \code{trace} (default 1).
#' @export
janine <- function(data, n_blocks, penalties = NULL,
                   control_optim = list(epsilon = 1e-4, max_iter = 20, trace = 1),
                   control_penalties = list(min_ration = 0.1, length = 25, diagonal = FALSE)
                   ) {

  n <- nrow(data)
  d <- ncol(X)
  S <- cov(as.matrix(data))

  ## this function willl be the method of a janine_fit object
  optim_janine <- function(penalty) {

    if (control_optim$trace == 1) {
      cat("\tsparsifying penalty =", penalty, "\r")
      flush.console()
    }

    weights <- matrix(1, d, d)
    if (!control_penalties$diagonal) diag(weights) <- 0
    cond <- FALSE
    iter <- 1
    objective <- numeric(control_optim$max_iter)
    objective[iter] <- Inf
    while(!cond) {
      ## M step (network structure estimation)
      net <- estimate_network(S, penalty, weights)
      sparsity <- 1 - sum(net$support) / (d**2)

      ## E step (latent block estimation)
      sbm <- estimate_block(net$support, n_blocks)
      weights <- (1 - sbm$connectProb)/sparsity
      if (!control_penalties$diagonal) diag(weights) <- 0

      ## Convergence assesment
      logdet <- determinant(net$Omega, logarithm = TRUE)$modulus
      objective[iter + 1] <- - (n/2) * (logdet - sum( diag( S %*% net$Omega )) ) + penalty * sum(abs(weights * net$Omega))
      cond <- abs(objective[iter + 1] - objective[iter])/abs(objective[iter + 1]) < control_optim$epsilon | iter + 1 > control_optim$max_iter
      iter <- iter + 1
    }

    list(network = net, membership = sbm, weights = weights, objective = objective[1:(iter - 1)])
  }


  ## default vector of penalties
  if (is.null(penalties)) {
    max_pen   <- max(abs(S))
    penalties <- 10^seq(from = log10( max(abs(S)) ),
                        to   = log10( max(abs(S))*control_penalties$min_ratio ),
                        len  = control_penalties$length)
  }

  ## call to optim
  res <- lapply(penalties, optim_janine)

  ## send back the results as a list if several penalties, a single fit other
  if (length(res) == 1) res <- res[[1]]
  res

}
