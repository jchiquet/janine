### TODO
## warm start with pprevious estimator in glasso
## check model selection criteria with l2 term

#' Janine: Just Another Network Infernce mEthod
#'
#' Iterate Adaptive Graphical-Lasso with binary SBM estimation to recover the adaptive weights
#'
#' @param data a n x d matrix of multivariate Gaussian observation
#' @param n_blocks integer for the target number of groups. If NULL (the default), will be chosen automatically internally by ICL in the SBM fit.
#' @param penalties a vector of postive real number in decreasing order tuning the network sparsity. The default (NULL) generates
#' an hopefully appropriate collection of penalties.
#' @param alpha a positive scalar tuning the mixture between the weighted-sparse penlaty and the trace-Laplacian regularisation.
#' @param control_penalties a list controling how \code{penalties} is generated, with three entries:
#' a double \code{min_ratio} (default 0.1), a integer \code{length} (default 20) and a logical \code{diagonal} (default TRUE)
#' indicating weither the diaognal should be penalized or not.
#' @param control_optim a list controling how the alternate optimization between adaptive graphical-Lasso and SBM is conducted
#' for each penalty level. Contains three entries: a convergence threshold \code{epsilon} (default to 1e-3),
#' a maximal number of iteration \code{max_iter} (default 20) and verbosity level \code{trace} (default 1).
#' @importFrom stats cov
#' @importFrom utils flush.console tail
#' @importFrom igraph laplacian_matrix graph_from_adjacency_matrix
#' @examples
#' ## Network settting
#' nNodes  <- 60
#' blockProp <- c(1/3, 1/3, 1/3)   # group proportions
#' nbBlock   <- length(blockProp) # number of blocks
#' connectParam <- diag(.4, nbBlock) + 0.01 # connectivity matrix: affiliation network
#' mySBM <- rggm::rSBM(nNodes, connectParam, blockProp)
#' Omega <- rggm::graph2prec(mySBM, cond_var = rep(1, nNodes), neg_prop = 0.5)
#' ## Multivariate Gaussian Vector generation
#' n <- 300
#' X <- rggm::rmgaussian(n, means = rep(0, nNodes), solve(Omega))
#' ## Network inference
#' fits <- janine(X, penalties = 0.1, control_optim = list(n_cores = 1))
#' plot(fits$models[[1]])
#' @export
janine <- function(data, partition = NULL, n_blocks = NULL, penalties = NULL, alpha = 0,
                   control_optim = list(), control_penalties = list()
                   ) {

  ## Initialization
  n <- nrow(data)
  d <- ncol(data)
  S <- cov(as.matrix(data))

  ctrl_optim <- list(epsilon = 1e-4, max_iter = 20, trace = 1, n_cores = 4)
  ctrl_optim[names(control_optim)] <- control_optim
  ctrl_penalties <- list(min_ratio = 0.1, length = 20, diagonal = TRUE, weighted = TRUE)
  ctrl_penalties[names(control_penalties)] <- control_penalties

  ## this function will be the optimization method of a janine_fit object
  optim_janine <- function(lambda) {

    if (ctrl_optim$trace == 1) {
      cat("\tAmount of regularisation =", lambda, "\r")
      flush.console()
    }

    weights <- matrix(1, d, d)
    if (!ctrl_penalties$diagonal) diag(weights) <- 0
    Laplacian <- matrix(0, d, d)
    cond <- FALSE
    iter <- 1
    objective <- numeric(ctrl_optim$max_iter)
    objective[iter] <- Inf
    while(!cond) {
      ## M step (network structure estimation)
      net <- estimate_network(S + lambda * alpha * Laplacian, lambda * (1-alpha), weights)
      sparsity <- 1 - sum(net$support) / (d**2)

      ## E step (latent block estimation)
      if (sparsity < 1) {
        sbm <- estimate_block(net$support, n_blocks, ctrl_optim$n_cores)
        predicted <- sbm$blockProb %*% sbm$connectParam %*% t(sbm$blockProb)
        if (ctrl_penalties$weighted) {
          weights <- (1 - predicted)/sparsity
          if (!ctrl_penalties$diagonal) diag(weights) <- 0
        }
        Laplacian <- predicted %>%
          graph_from_adjacency_matrix(weighted = TRUE, diag = FALSE, mode = "undirected") %>%
          laplacian_matrix(sparse = FALSE)
      } else {
        sbm <- NA
      }

      ## Convergence assesment
      loglik <- (n/2) * (determinant(net$Omega, logarithm = TRUE)$modulus - trace(S, net$Omega))
      sparse_reg <- sum(abs(weights * net$Omega))
      graph_reg  <- trace(Laplacian,  net$Omega)
      objective[iter + 1] <- - (2/n) * loglik + lambda * ( (1-alpha) * sparse_reg + alpha * graph_reg )
      cond <- abs(objective[iter + 1] - objective[iter])/abs(objective[iter + 1]) < ctrl_optim$epsilon | iter + 1 > ctrl_optim$max_iter
      iter <- iter + 1
    }

    n_edges <- sum(net$support)/2
    BIC   <- -2 * loglik + log(n) * n_edges
    EBIC  <- BIC + ifelse(n_edges > 0, n_edges * log(.5 * d*(d - 1)/n_edges), 0)

    structure(
      list(network = net, block = sbm, weights = weights, penalty = lambda,
         objective = objective[1:(iter - 1)], loglik = loglik, BIC = BIC, EBIC = EBIC), class = "janine_fit")
  }

  ## default vector of penalties
  if (is.null(penalties)) {
    max_pen   <- max(abs(S[upper.tri(S, diag = ctrl_penalties$diagonal)]))
    penalties <- 10^seq(from = log10( max_pen ),
                        to   = log10( max_pen*ctrl_penalties$min_ratio ),
                        len  = ctrl_penalties$length)
  }

  if (ctrl_optim$trace > 0) cat("\n Adjusting", length(penalties), "SBM-structured GGM with sparse adaptive regularisation \n")

  ## call to optim
  models <- lapply(penalties, optim_janine)

  # formatting output
  criteria <- data.frame(
    penalty   = penalties,
    sparsity  = sapply(models, function(model) 1 - sum(model$net$support) / (d**2)),
    n_edges   = sapply(models, function(model) sum(model$net$support)/2),
    loglik    = sapply(models, function(model) sum(model$loglik)),
    objective = sapply(models, function(model) utils::tail(model$objective, 1)),
    BIC       = sapply(models, function(model) sum(model$BIC)),
    EBIC      = sapply(models, function(model) sum(model$EBIC))
  )

  structure(list(models = models, criteria = criteria), class = "janine_collection")
}

