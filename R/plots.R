#' @title Plotting the confidence intervals for the differential parameters
#' @method plot inferchange.ci
#' @description Plotting method for S3 objects of class \code{inferchange.ci}.
#' Displays the original and de-sparsified estimators of differential parameter for each coefficient (x-axis) and the simultaneous confidence interval constructed around it
#' @param x \code{inferchange.ci} object
#' @param ... additional arguments
#' @return A plot produced as per the input arguments
#' @seealso \link{ci_delta}
#' @importFrom graphics plot arrows points abline legend
#' @export
plot.inferchange.ci <- function(x, ...){

  ci <- x$ci
  p <- nrow(ci)
  rn <- range(ci)
  plot(1:p, rep(NA, p), ylim = c(rn[1] - 2, rn[2]), xlab = "Variables", ylab = "", main = paste( 100 * (1 - x$alpha), "% simultaneous confidence intervals", sep = ''))
  arrows(1:p, ci[, 1], 1:p, ci[, 2], length = 0.05, angle = 90, code = 3, col = 8)
  points(x$delta.check, col = 1)
  points(x$delta.hat, col = 2, pch = "x")
  positives <- which(ci[, 1] * ci[, 2] > 0)
  if(length(positives) > 0) arrows(positives, ci[positives, 1], positives, ci[positives, 2], length = 0.05, angle = 90, code = 3, col = 6)
  abline(h = 0, col = 8)
  legend("bottomright", pch = c("o", "x"), col = 1:2, legend = c("de-sparsified", "original"), bty = "n")

}

#' @title Plotting the segmentation result
#'
#' @description
#' Plotting the segmentation result, which is of class \code{inferchange.cp},
#' e.g. returned by function \code{McScan}
#'
#'
#' @param x \code{inferchange.cp} object
#' @param ... additional arguments
#'
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @export
plot.inferchange.cp <- function(x, ...) {
  X = attr(x, "X")
  y = attr(x, "y")
  n = nrow(X)
  p = ncol(X)
  Xy =  X * matrix(y, nrow = n, ncol = p)
  Xy = apply(Xy, 2, cumsum)
  k  = 1:n
  csXy = matrix(sqrt((n - k)/n/k) + sqrt(k/n/(n - k)), nrow = length(k),
                    ncol = ncol(Xy)) * Xy[k, ]
  csXy = csXy - matrix(sqrt(k/n/(n - k)), ncol = 1) %*% Xy[n, ]
  df   = utils::stack(as.data.frame(csXy, row.names = NULL, col.names = NULL))
  df$x = rep(seq_len(nrow(csXy)), ncol(csXy))
  values <- ind <- ncp <- value <- NA  # dummy codes to avoid notes
  plt = ggplot(df, aes(x = x, y = values, group = ind)) +
    geom_line(color = "gray", na.rm = TRUE)  +
    geom_vline(xintercept = x$cp, linetype = "dashed", color = "red") +
    xlab("Index of sample") + ylab("CUSUM of Cov(X, y)") +
    theme_minimal() +
    ggtitle("Estimated change points (red dashed vertical lines)")
  if (!is.null(attr(x, "solution_path"))) {
    sopa = attr(x, "solution_path")
    sopa = sopa[is.finite(unlist(sopa[,2])), ]
    id   = attr(x, "selected_solution")
    plt2 = ggplot(data.frame("value" = unlist(sopa[,2]),
                             "ncp" = sapply(sopa[,3], length)),
                  aes(x = ncp, y = value)) +
      geom_line() + geom_count(alpha = 0.5) +
      geom_point(data = data.frame("ncp" = length(unlist(sopa[id,"cps"])),
                                   "value" = unlist(sopa[id,"val"])),
                 shape = 4, size = 3, color = "red")+
      ylab("Evidence of undetected change points") +
      xlab("Number of detected change points") +
      guides(size = guide_legend(title = "Number of solutions")) +
      scale_size_continuous(breaks = round) +
      theme_minimal() +
      theme(axis.line = element_line())
    gridExtra::grid.arrange(plt, plt2, nrow = 2)
  } else { print(plt) }
}
