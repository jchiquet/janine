#' #' Extract and plot the network (pArtial correlation, support or inverse covariance) from a \code{janine} object
#' #'
#' #' @name plot.janine_fit
#' #'
#' #' @param x an S3 object with class janine_fit
#' #' @param type character. Value of the weigth of the edges in the network, either "partial_cor" (partial correlation) or "support" (binary). Default is \code{"partial_cor"}.
#' #' @param output the type of output used: either 'igraph' or 'corrplot'. Default is 'igraph'.
#' #' @param edge.color 2D numeric. Color for positive/negative edges. Default is c("#F8766D", "#00BFC4"). Only relevant for igraph output.
#' #' @param node.labels vector of character. The labels of the nodes. The Default will use the column names ot the data matrix.
#' #' @param remove.isolated if \code{TRUE}, isolated node are remove before plotting. Only relevant for igraph output.
#' #' @param layout an optional igraph layout. Only relevant for igraph output.
#' #' @param plot logical. Should the final network be displayed or only sent back to the user. Default is \code{TRUE}.
#' #' @param ... Not used (S3 compatibility).
#' #'
#' #' @return Send back an invisible object (igraph or Matrix, depending on the output chosen) and optionaly displays a graph (via igraph or corrplot for large ones)
#' #' @export
#' plot.janine_fit <-
#'   function(x,
#'            type            = c("partial_cor", "support"),
#'            output          = c("igraph", "corrplot"),
#'            edge.color      = c("#F8766D", "#00BFC4"),
#'            remove.isolated = FALSE,
#'            node.labels     = NULL,
#'            layout          = layout_in_circle,
#'            plot            = TRUE, ...) {
#'     invisible(x$plot_network(type, output, edge.color, remove.isolated, node.labels, layout, plot))
#'   }
