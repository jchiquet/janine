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
#' @importFrom stats cov
#' @importFrom utils flush.console tail
#' @export
janine <- function(data, n_blocks, penalties = NULL,
                   control_optim = list(epsilon = 1e-4, max_iter = 20, trace = 1, n_cores = 4),
                   control_penalties = list(min_ratio = 0.1, length = 20, diagonal = TRUE)
                   ) {

  n <- nrow(data)
  d <- ncol(data)
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
      if (sparsity < 1) {
        sbm <- estimate_block(net$support, n_blocks, control_optim$n_cores)
        weights <- (1 - sbm$connectProb)/sparsity
        if (!control_penalties$diagonal) diag(weights) <- 0
      } else {
        sbm <- NA
      }

      ## Convergence assesment
      loglik <- (n/2) * (determinant(net$Omega, logarithm = TRUE)$modulus - sum( diag( S %*% net$Omega )) )
      objective[iter + 1] <- - loglik + penalty * sum(abs(weights * net$Omega))
      cond <- abs(objective[iter + 1] - objective[iter])/abs(objective[iter + 1]) < control_optim$epsilon | iter + 1 > control_optim$max_iter
      iter <- iter + 1
    }

    n_edges <- sum(net$support)/2
    BIC   <- loglik - .5 * log(n) * n_edges
    EBIC  <- BIC - .5 * ifelse(n_edges > 0, n_edges * log(.5 * d*(d - 1)/n_edges), 0)

    structure(
      list(network = net, membership = sbm, weights = weights,
         objective = objective[1:(iter - 1)], loglik = loglik, BIC = BIC, EBIC = EBIC), class = "janine_fit")
  }

  ## default vector of penalties
  if (is.null(penalties)) {
    max_pen   <- max(abs(S[upper.tri(S, diag = control_penalties$diagonal)]))
    penalties <- 10^seq(from = log10( max_pen ),
                        to   = log10( max_pen*control_penalties$min_ratio ),
                        len  = control_penalties$length)
  }

  if (control_optim$trace > 0) cat("\n Adjusting", length(penalties), "SBM-structured GGM with sparse adaptive regularisation \n")

  ## call to optim
  models <- lapply(penalties, optim_janine)

  # formatting output
  criteria <- data.frame(
    penalty   = penalties,
    sparsity  = sapply(models, function(model) 1 - sum(model$net$support) / (d**2)),
    n_edges   = sapply(models, function(model) sum(model$net$support)/2),
    loglik    = sapply(models, function(model) sum(model$loglik)),
    objective = sapply(models, function(model) tail(model$objective, 1)),
    BIC       = sapply(models, function(model) sum(model$BIC)),
    EBIC      = sapply(models, function(model) sum(model$EBIC))
  )

  structure(list(models = models, criteria = criteria), class = "janine_collection")
}

