library(rggm)
library(GREMLIN)
library(blockmodels)
library(glassoFast)
library(Matrix)
library(corrplot)

nNodes  <- 90
blockProp <- c(1/3, 1/3, 1/3)   # group proportions
nbBlock   <- length(blockProp) # number of blocks
connectParam <- diag(.4, nbBlock) + 0.01 # connectivity matrix: affiliation network

## Graph Sampling
mySBM <- rSBM(nNodes, connectParam, blockProp)

## Sampling Gaussian weights
mu_within   <- 4 ; sigma_within  <- .5
mu_between  <- 2 ; sigma_between <- .5
theta <- list()
theta$mu    <- matrix(mu_between   , nbBlock, nbBlock); diag(theta$mu)    <- mu_within    # means
theta$sigma <- matrix(sigma_between, nbBlock, nbBlock); diag(theta$sigma) <- sigma_within # sd
mySBM_Gaussian <- rWeightSBM(mySBM, "gaussian", theta)

## Building precision and covariance matrices
Omega <- graph2prec(mySBM_Gaussian, cond_var = rep(1, nNodes))
Sigma <- solve(Omega)

## Sampling Gaussian data
n <- 1000
means <- rep(0, ncol(Sigma))
X <- rmgaussian(n, means, Sigma)
Sigma_hat <- cov(as.matrix(X))


out <- janine(X, 0.075, nbBlock)


####
## JANINE: Smart SBM via GREMLIN + GLASSO
##
glasso_out <- estimate_network(S = Sigma_hat, lambda = 0.075)
Omega_GL <- as.matrix(glasso_out$wi)
Sigma_GL <- as.matrix(glasso_out$w)

par(mfrow = c(3,2))
corrplot(Omega_GL, is.corr = FALSE, tl.pos = "n", method = 'color', type = "upper", diag = FALSE)
corrplot(Sigma_GL, is.corr = FALSE, tl.pos = "n", method = 'color', type = "upper", diag = FALSE)
hist(Omega_GL[upper.tri(Omega_GL) * Omega_GL != 0 ], main = "histogram of non-null entries in Omega hat")
hist(Sigma_GL[upper.tri(Sigma_GL) * Sigma_GL != 0 ], main = "histogram of non-null entries in Sigma hat")
hist(Omega[upper.tri(Omega) * Omega != 0 ], main = "histogram of non-null entries in Omega star")
hist(Sigma[upper.tri(Sigma) * Sigma != 0 ], main = "histogram of non-null entries in Sigma star")

# net_support <- defineNetwork((Omega_GL != 0)*1, "adj", "block", "block")
# net_weights <- defineNetwork(Omega_GL, "adj", "block", "block")
# SBM_hat <- multipartiteBMFixedModel(list(net_support), c("bernoulli"), c("block"), v_K = 3)

net_support <- (Omega_GL != 0)*1; diag(net_support) <- 0
structure <- estimate_block(net_support, nbBlock)
corrplot(structure$connectProb, method = "color", is.corr = FALSE, tl.pos = "n", cl.pos = "n")
density <- sum(net_support)/ nNodes**2
hist((1 - structure$connectProb)/(1-density))

glasso_out <- estimate_network(S = Sigma_hat, lambda = 0.075, W = 1-structure$connectProb)
Omega_GL <- as.matrix(glasso_out$wi)
Sigma_GL <- as.matrix(glasso_out$w)
par(mfrow = c(3,2))
corrplot(Omega_GL, is.corr = FALSE, tl.pos = "n", method = 'color', type = "upper", diag = FALSE)
corrplot(Sigma_GL, is.corr = FALSE, tl.pos = "n", method = 'color', type = "upper", diag = FALSE)
hist(Omega_GL[upper.tri(Omega_GL) * Omega_GL != 0 ], main = "histogram of non-null entries in Omega hat")
hist(Sigma_GL[upper.tri(Sigma_GL) * Sigma_GL != 0 ], main = "histogram of non-null entries in Sigma hat")
hist(Omega[upper.tri(Omega) * Omega != 0 ], main = "histogram of non-null entries in Omega star")
hist(Sigma[upper.tri(Sigma) * Sigma != 0 ], main = "histogram of non-null entries in Sigma star")

