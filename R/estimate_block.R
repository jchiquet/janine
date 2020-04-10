#' Estimation of the latent block organisation of the network
#'
#' The underlying network is assumed to be drawn from a Stochastic Bloc Model.
#' This function uses a variational EM algorithm implemented in the package blockmodels to estimate such a model.
#'
#' @param adj_matrix a symmetric weighted adjacency matrix
#' @param n_blocks  integer, the number of blocks. Default is NULL, in which case the best model in terms of ICL is returned.
#' @param n_cores  integer, the number of cores used for initializing the SBM exploration
#' @import blockmodels
#' @export
estimate_block <- function(adj_matrix, n_blocks = NULL, n_cores = 1){

  stopifnot(sum(adj_matrix) != 0)

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
    blockProb   = SBM_fits$memberships[[n_blocks]]$Z,
    connectProb = (SBM_fits$prediction(n_blocks) + t(SBM_fits$prediction(n_blocks)))/2,
    membership  = factor(apply(SBM_fits$memberships[[n_blocks]]$Z, 1, which.max), levels = 1:n_blocks),
    parameters  = SBM_fits$model_parameters[[n_blocks]]
  )
}

#' #' @import GREMLIN
#' estimate_block.list <- function(list_adj_matrix, n_blocks = NULL, n_cores = 1){
#'
#'   stopifnot(sum(adj_matrix) != 0)
#'
#'   S <- length(list_adj_matrix)
#'
#'   group_names <- names(list_adj_matrix)
#'
#'   if (is.null(group_names)) group_names <- paste0("group_", 1:S)
#'
#'   list_network <- mapply(FUN = function(mat, name) { GREMLIN::defineNetwork(mat, "adj", name) })
#'
#'   ## call to GREMLIN avoiding any screen output
#'   mp_SBM <- ifelse(is.null(n_blocks), multipartiteBM, multipartiteBMFixedModel)
#'   {
#'     sink("/dev/null")
#'     res <- mp_SBM(
#'             list_Net  = list_network,
#'             v_distrib = rep("bernoulli",S),
#'             verbose   = FALSE,
#'             nbCores   = n_cores)
#'     sink()
#'   }
#'
#'   n_blocks <- ifelse(is.null(n_blocks), which.max(SBM_fits$ICL), n_blocks)
#'
#'   tau <- res$fittedModel[[1]]$paramEstim$tau
#'
#'   Matrix::bdiag(tau)
#'
#'   res <- list(
#'     blockProb   = SBM_fits$memberships[[n_blocks]]$Z,
#'     connectProb = (SBM_fits$prediction(n_blocks) + t(SBM_fits$prediction(n_blocks)))/2,
#'     membership  = factor(apply(SBM_fits$memberships[[n_blocks]]$Z, 1, which.max), levels = 1:n_blocks)
#'   )
#' }
