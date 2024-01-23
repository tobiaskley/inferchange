
#' Find suitable grid of lambda values
#'
#' @param lambda_max TODO Add Description
#' @param n TODO Add Description
#' @param p TODO Add Description
#' @param K TODO Add Description
#'
#' @return TODO Add Description
#' @export
#'
#' @examples
#' # TODO Add Example or remove this
lambdapath <- function(lambda_max, n, p, K = 100) {

  epsilon <- ifelse(n < p, 0.01, 1e-04)
  lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon),
                              length.out = K)), digits = 10)
  return(lambdapath)
}
