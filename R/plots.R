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


#' @param x \code{inferchange.cp} object
#' @param ... additional arguments
#' @importFrom graphics abline matplot par
#' @export
plot.inferchange.cp <- function(x, ...) {
  Xy = attr(x,"X") * matrix(attr(x,"y"), nrow = nrow(attr(x,"X")),
                            ncol = ncol(attr(x,"X")))
  matplot(Xy, type = "l", col = "gray",
          ylab = "Covariance of X and y",
          xlab = "Sample index", lty = 1,
          main = "Estimated change points")
  abline(v = x$cp, col = "blue", lty = "dotted")
  if (!is.null(attr(x, "solution_path"))) {
    sopa = attr(x, "solution_path")
    plt <- ggplot(data.frame(value = unlist(sopa[,2]),
                             ncp = sapply(sopa[,3], length)),
                  aes(x=ncp, y=value)) +
      geom_line() + geom_count(alpha = 0.5) +
      geom_point(data = data.frame(ncp = length(unlist(sopa[id,"cps"])),
                                   value = unlist(sopa[id,"val"])),
                 shape = 4, size = 3, color = "red")+
      ylab("Evidence of undetected change points") +
      xlab("Number of detected change points") +
      guides(size = guide_legend(title = "Number of solutions")) +
      scale_size_continuous(breaks = round) +
      theme_minimal() +
      theme(axis.line = element_line())
    print(plt)
  }
}
