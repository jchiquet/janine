
is_janine_fit <- function(x) (class(x) == "janine_fit")
is_janine_collection <- function(x) (class(x) == "janine_collection")

#' Extract and plot the parameters (partial correlation, support or inverse covariance) from a \code{janine_fit} object
#'
#' @name plot.janine_fit
#'
#' @param x an S3 object with class janine_fit
#' @param type character. Value of the weigth of the edges in the network, either "partial_cor" (partial correlation), "precision", "covariance" or "support" (binary). Default is \code{"partial_cor"}.
#' @param title an optional string for the title of the plot.
#' @param ... Not used (S3 compatibility).
#' @import corrplot
#' @return Send back an invisible object (igraph or Matrix, depending on the output chosen) and optionaly displays a graph (via igraph or corrplot for large ones)
#' @export
plot.janine_fit <- function(x, type= c("partial_cor", "precision", "covariance", "support"), title = paste(match.arg(type), "matrix"), ...) {

  stopifnot(is_janine_fit(x))
  type   <- match.arg(type)
  net <- switch(
      match.arg(type),
      "support"     = x$network$support,
      "covariance"  = x$network$Sigma,
      "precision"   = x$network$Omega,
      "partial_cor" = {
        tmp <- -x$network$Omega / tcrossprod(sqrt(diag(x$network$Omega))); diag(tmp) <- 1
        tmp
      }
    ) %>% as.matrix()

  if (ncol(net) > 50) colnames(net) <- rownames(net) <- rep(" ", ncol(net))
  corrplot(as.matrix(net), method = "color", is.corr = FALSE, tl.pos = "td", cl.pos = "n", tl.cex = 0.5, type = "upper", diag = FALSE, title = title)

}

#' Model selection from a collection of fits
#'
#' @param x an object with class janine_collection
#' @param crit a character for the criterion used to performed the selection. Either
#' "BIC", "EBIC" or "loglik". Default is \code{BIC}.
#' @return  Send back an object with class \code{janine_fit}
#' @export
select_model <- function(x, crit = c("BIC", "EBIC", "loglik")) {UseMethod("select_model", x)}
#' @export
select_model.janine_collection <- function(x, crit = c("BIC", "EBIC", "loglik")){
  stopifnot(is_janine_collection(x))
  crit <- match.arg(crit)
  id <- 1
  if (length(x$criteria[[crit]]) > 1) {
    id <- which.min(x$criteria[[crit]])
    if (crit == "loglik")  id <- which.max(x$criteria[[crit]])
  }
  model <- x$models[[id]]
  model
}

#' Plot model selection criteria for a collection of fits
#'
#' @param x an S3 object with class janine_collection
#' @param criteria vector of characters. The criteria to plot in c("loglik", "BIC", "EBIC").
#' @param log.x logical: should the x-axis be repsented in log-scale? Default is \code{TRUE}.
#' @param ... use for S3 compatibility
#' @importFrom dplyr select group_by mutate
#' @importFrom tidyr gather
#' @importFrom tidyselect all_of
#' @import ggplot2
#' @export
plot.janine_collection <- function(x, criteria = c("loglik", "BIC", "EBIC"), log.x = TRUE, ...) {
  stopifnot(is_janine_collection(x))
  dplot <- x$criteria %>%
    mutate(loglik = -2 * .data$loglik) %>%
    dplyr::select(c("penalty", all_of(criteria))) %>%
    tidyr::gather(key = "criterion", value = "value", -.data$penalty) %>%
    dplyr::group_by(.data$criterion)
  p <- ggplot(dplot, aes_string(x = "penalty", y = "value", group = "criterion", colour = "criterion")) +
    geom_line() + geom_point() + ggtitle("Model selection criteria") + theme_bw() + xlab("penalty")
  if (log.x) p <- p + ggplot2::coord_trans(x = "log10")
  p
}


#' Extract the regularization path of a collection of fits
#'
#' @name coefficient_path
#' @param x an object with class \code{janine_collection}
#' @param type a character, should we extract the path of covariance, precision or partial correlation coefficients ? Default is \code{partial_cor}.
#' @importFrom dplyr filter bind_rows
#' @importFrom stats setNames
#' @return  Send back a tibble/data.frame.
#' @export
coefficient_path <- function(x, type = c("partial_cor", "precision", "covariance")) {UseMethod("coefficient_path", x)}
#' @importFrom rlang .data
#' @export
coefficient_path.janine_collection <- function(x, type = c("partial_cor", "precision", "covariance")) {
  type <- match.arg(type)
  lapply(x$models, function(model) {
    G <- switch(
      type,
      "covariance"  = model$network$Sigma,
      "precision"   = model$network$Omega,
      "partial_cor" = {
        tmp <- -model$network$Omega / tcrossprod(sqrt(diag(model$network$Omega))); diag(tmp) <- 1
        tmp
      }
    ) %>% as.matrix()
    if(is.null(colnames(G))) colnames(G) <- 1:ncol(G)
    if(is.null(rownames(G))) rownames(G) <- 1:ncol(G)

    setNames(
      cbind(
        expand.grid(colnames(G), rownames(G)),
        as.vector(G)), c("Node1", "Node2", "Coeff")
    ) %>%
      mutate(Penalty = model$penalty,
             Node1   = as.character(.data$Node1),
             Node2   = as.character(.data$Node2),
             Edge    = paste0(.data$Node1, "|", .data$Node2)) %>%
      dplyr::filter(.data$Node1 < .data$Node2)
  }) %>% bind_rows()
}

