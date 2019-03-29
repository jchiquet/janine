context("test-sbm-laplace")

library(igraph)
library(aricode)

error <- function(theta_hat, theta, sort = FALSE) {
  if (sort) {
    err <- sum( (sort(theta_hat) -  sort(theta) )^2)
  } else {
    err <- sum( (theta_hat -  theta )^2)
  }
  err
}

set.seed(97857)
Lambda <- matrix(c(1, 2, 3, 2, 1, 4, 3, 4, 1), 3, 3)
p <- 200
alpha <- c(1/2, 1/4, 1/4)
K <- length(alpha)
sbm   <- rSBM_Laplace(p, Lambda, alpha)

clustering  <- V(sbm)$membership
data_matrix <- sbm %>% as_adj(attr = "weight")


test_that("SBM Laplace is working", {

  tol_ARI <- .9
  tol_ref <- 1e-2

  cl0 <- kmeans(data_matrix, centers = K, nstart = 10)$cl
  out <- estimate_structure(data_matrix, K, cl0)

  ## expect SBM clustering is working
  expect_gt(ARI(out$membership, clustering), tol_ARI)

  ## expect SBM clustering is working better thant initialization
  expect_gt(ARI(out$membership, clustering), ARI(cl0, clustering))

  ## expect estimation works for mixture parameters
  expect_lt(error(out$theta$alpha, alpha, sort = TRUE), tol_ref)

  ## expect estimation works for connectivity parameters
  expect_lt(error(out$theta$Lambda, Lambda, sort = TRUE), tol_ref)

})
