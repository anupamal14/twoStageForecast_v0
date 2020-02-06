#' @title  Metric for model selection
#'
#' @description Finds psi = n x ln(MSE) + 2 x p, an alternative to AIC.
#'     psi, or penalized lack-of-fit, is used for model selection.
#'
#' @param MSE Mean Square Error of the fitted model.
#'
#' @param n length of the data.
#'
#' @param p number of parameters estimated by the model.
#'
#'
#' @examples
#'


twoStage.penalizedLikelihood <- function(n, MSE, p) {
  psi <- n*log(MSE) + 2*p
}
