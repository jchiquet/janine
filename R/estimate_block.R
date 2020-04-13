#' Estimation of the latent block organisation of the network
#'
#' The underlying network is assumed to be drawn from a Stochastic Bloc Model.
#' This function uses a variational EM algorithm implemented in the package blockmodels to estimate such a model.
#'
#' @param adj_matrix a symmetric weighted adjacency matrix
#' @param partition a factor indicating a known partition of the variables to be respected during the clustering.
#' If NULL (the default), to predfined partition is consider.
#' @param n_blocks integer, the number of blocks. Default is NULL, in which case the best model in terms of ICL is returned.
#' @param n_cores integer, the number of cores used for initializing the SBM exploration
#' @import blockmodels
#' @import GREMLIN
#' @importFrom Matrix symmpart
#' @export
estimate_block <- function(adj_matrix, partition = NULL, n_blocks = NULL, n_cores = 1){

  stopifnot(sum(adj_matrix) != 0)

  ## standard stochastic blockmodel
  if (is.null(partition)) {
    SBM_fits <- BM_bernoulli(
      membership_type = "SBM",
      adj             = adj_matrix,
      verbosity       = 0,
      plotting        = "",
      explore_min     = ifelse(is.null(n_blocks), 4, n_blocks),
      ncores          = n_cores
    )

    ## call to block model avoiding any output
    {
      sink("/dev/null"); SBM_fits$estimate() ; sink()
    }

    n_blocks <- ifelse(is.null(n_blocks), which.max(SBM_fits$ICL), n_blocks)

    res <- list(
      blockProb    = SBM_fits$memberships[[n_blocks]]$Z,
      membership   = factor(apply(SBM_fits$memberships[[n_blocks]]$Z, 1, which.max), levels = 1:n_blocks),
      connectParam = Matrix::symmpart(SBM_fits$model_parameters[[n_blocks]]$pi)
    )
  } else {
  ## Multipartite blockmodel is sued to handle the fixed partition

    networks <- adj_matrix2GREMLIN_format(adj_matrix, partition)
    n_grp <- nlevels(partition)

    fit <- GREMLIN::multipartiteBM(
      networks,
      v_distrib = rep('bernoulli', length(networks)),
      namesFG = levels(partition),
      nbCores = n_cores,
      verbose = FALSE
    )$fittedModel[[1]]$paramEstim

    cols <- vector("list", n_grp)
    for (i in 1:n_grp) cols[[i]] <- do.call(rbind, fit$list_theta[(1:n_grp) + n_grp*(i-1)])
    connectParam <- do.call(cbind, cols)

    membership <- fit$Z
    for (i in 2:n_grp) membership[[i]] <-  membership[[i]] + cumsum(fit$v_K)[i-1]
    membership <- unlist(membership)

    Tau <- as.matrix(Matrix::bdiag(fit$tau))

    res <- list(
      blockProb    = Tau,
      membership   = membership,
      connectParam = connectParam
    )
  }
  res
}

