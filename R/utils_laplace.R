#' rSBMLaplace
#'
#' Generate a symmetric matrix under the stochastic block model
#' with a Laplace distribution as the emission law of the edges.
#'
#' @param p number of nodes
#' @param Lambda a QxQ matrix of connectivity parameters
#' @param alpha a size-Q vector of mixture
#' @return an igraph object: undirected weigthed graph with a attribute "membership" for the vertices
#' @import Matrix igraph
#' @importFrom rmutil rlaplace
#' @export
rSBMLaplace <- function(p, Lambda, alpha = c(1)) {
  Q <- length(alpha)
  Z <- t(rmultinom(p, 1, prob = alpha))
  A <- matrix(0, p, p)
  subset <- upper.tri(A)
  lambda <- (Z %*% Lambda %*% t(Z))[subset]
  x <- rmutil::rlaplace(n=p*(p-1)/2, m = 0, s = lambda)
  A[subset] <- x
  A <- A + t(A)
  mySBM <- graph_from_adjacency_matrix(A, weighted = TRUE, mode = "undirected")
  vertex_attr(mySBM, "membership") <- Z %*% 1:Q
  mySBM
}


#
# ############# SIMPLE EXAMPLE ON A COMMUNITY NETWORK
# p_in <- 0.5
# p_out <- 0.01
# pi.comm <- matrix(p_out, 3, 3); diag(pi.comm) <- p_in
# n <- 100
# alpha <- c(.5, .25, .25)
#
# G <- rNetwork(n, pi, alpha)
# plot(G, vertex.color=V(G)$membership, vertex.label=NA)
#
# ## run VEM for SBM for various number of clusters
# seq.Q <- 2:6
# library(pbmcapply)
# VEM_out <- pbmclapply(seq.Q, function(Q) {
#   VEM_SBM(G, Q)
# })
# vBIC <- sapply(VEM_out, function(model) model$vBIC)
# vICL <- sapply(VEM_out, function(model) model$vICL)
#
# plot(seq.Q,vBIC, type="l", col="red")
# lines(seq.Q,vICL, col="blue")
# legend("bottomright", legend=c("vBIC","vICL"), lty=c(1,1), col=c("blue", "red"))
# best.model <- VEM_out[[which.max(vICL)]]
# aricode::ARI(best.model$membership, V(G)$membership) # yes !
#
# ############# ANALYSIS of the FRENCH POLITICAL BLOGSPHERE
# library(sand)
# fblog <- upgrade_graph(fblog)
# seq.Q <- 2:15
# library(pbmcapply)
# VEM_out <- pbmclapply(seq.Q, function(Q) {
#   VEM_SBM(fblog, Q)
# })
# vBIC <- sapply(VEM_out, function(model) model$vBIC)
# vICL <- sapply(VEM_out, function(model) model$vICL)
# plot(seq.Q,vBIC, type="l", col="red")
# lines(seq.Q,vICL, col="blue")
# legend("bottomright", legend=c("vBIC","vICL"), lty=c(1,1), col=c("blue", "red"))
# best.model.fblog <- VEM_out[[which.max(vICL)]]
# aricode::ARI(best.model.fblog$membership, V(fblog)$PolParty) # yes !
#
# plot(fblog, vertex.color = best.model.fblog$membership, vertex.label=NA)



