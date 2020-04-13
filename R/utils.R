#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @importFrom GREMLIN defineNetwork
NULL

# trace operator
trace <- function(A, B) { sum ( diag( A %*% B ) ) }


# turn an adjacency matrix with a partition (functional groups) to a list of networks
adj_matrix2GREMLIN_format <- function(adj_matrix, partition)  {

  groups <- levels(partition)
  pairs  <- expand.grid(groups, groups)
  networks <- vector("list", nrow(pairs))

  for (i in 1:nrow(pairs)) {
    networks[[i]] <- GREMLIN::defineNetwork(
      adj_matrix[partition  == pairs[i, 1], partition == pairs[i, 2]],
      ifelse (pairs[i, 1] == pairs[i, 2], "adj", "inc"),
      as.character(pairs[i, 1]),
      as.character(pairs[i, 2])
    )
  }
  networks
}
