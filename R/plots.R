#' @title Plotting the confidence intervals for the differential parameters
#' @method plot inferchange.ci
#' @description Plotting method for S3 objects of class \code{inferchange.ci}.
#' Displays the de-sparsified estimators of the differential parameters and the simultaneous confidence intervals for the respective parameters
#' @param x \code{inferchange.ci} object
#' @param ... additional arguments
#' @return A plot produced as per the input arguments
#' @seealso \link[inferchange]{ci.delta}
#' @examples
#' \donttest{
#' }
#' @export
plot.inferchange.ci <- function(x, ...){

  ci <- x$ci
  p <- nrow(ci)
  plot(1:p, rep(NA, p), ylim = range(ci), xlab = "variables", ylab = "", main = "confidence intervals")
  arrows(1:p, ci[, 1], 1:p, ci[, 2], length = 0.05, angle = 90, code = 3, col = 8)
  points(delta.check, col = 1)
  # points(delta.hat, col = 4, pch = 3)
  positives <- which(ci[, 1] * ci[, 2] > 0)
  if(length(positives) > 0) arrows(positives, ci[positives, 1], positives, ci[positives, 2], length = 0.05, angle = 90, code = 3, col = 2)
  abline(h = 0, col = 8)
}
