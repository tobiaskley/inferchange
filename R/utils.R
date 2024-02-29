
#' Find suitable grid of lambda values
#'
#' Generates K lambda values from `lambda_max` to 0.01 (if `n < p`) or 0.0001
#' (if `n >= p`) in decreasing order that are equaly spaced on the log scale.
#' This is the grid of lambda values employed as in `glmnet`.
#' In applications lambda_max is chosen as
#'
#' @param lambda_max maximum value of lambda values on grid
#' @param n value `n`
#' @param p value `p`
#' @param K number of lambda values
#'
#' @return vector with lambda values
#' @export
#'
#' @examples
#' lambdapath(10, 100, 200, K = 5)
lambdapath <- function(lambda_max, n, p, K = 100) {

  epsilon <- ifelse(n < p, 0.01, 1e-04)
  lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon),
                              length.out = K)), digits = 10)
  return(lambdapath)
}

#' @export
print.inferchange.cp <- function(x, ...) {
  cat("Detected", length(x$cp), "change points at", x$cp, "\n")
}

#' @export
print.inferchange.ci <- function(x, ...) {
  cat(paste(100 * (1 - x$alpha), "% simultaneous confidence intervals not containing 0", sep = ''), "\n")
  positives <- which(x$ci[, 1] * x$ci[, 2] > 0)
  print(x$ci[positives, ])
}

