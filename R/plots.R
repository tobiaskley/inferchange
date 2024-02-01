#' @title Plotting the confidence intervals for the differential parameters
#' @method plot inferchange.ci
#' @description Plotting method for S3 objects of class \code{inferchange.ci}.
#' Displays the de-sparsified estimator of differential parameter for each coefficient (x-axis) and the simultaneous confidence interval constructed around it
#' @param x \code{inferchange.ci} object
#' @param ... additional arguments
#' @return A plot produced as per the input arguments
#' @seealso \link[inferchange]{ci.delta}
#' @importFrom graphics plot arrows points abline
#' @examples
#' \donttest{
#' }
#' @export
plot.inferchange.ci <- function(x, ...){

  ci <- x$ci
  p <- nrow(ci)
  plot(1:p, rep(NA, p), ylim = range(ci), xlab = "Variables", ylab = "", main = paste( 100 * x$alpha, "% simultaneous confidence intervals", sep = ''))
  arrows(1:p, ci[, 1], 1:p, ci[, 2], length = 0.05, angle = 90, code = 3, col = 8)
  points(x$delta.check, col = 1)
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
  df = expand.grid(sample_index=1:n, coordinate=1:p)
  df$covariance = as.vector(Xy)
  plt = ggplot(df, aes(sample_index, coordinate, fill = covariance)) +
    geom_tile() + scale_fill_viridis() +
    theme_minimal() + geom_vline(xintercept = x$cp)
  print(plt)
  if (!is.null(attr(x, "solution_path"))) {
    sopa = attr(x, "solution_path")
    id   = attr(x, "selected_solution")
    plt  = ggplot(data.frame(value = unlist(sopa[,2]),
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
