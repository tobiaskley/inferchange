

#' Function for change point detection and inference about changes
#'
#' This is the `workhorse` function that implements all the methodologies
#' developed in Cho, Kley & Li (2024) that detect multiple change points
#' and infer about the changes in a high-dimensional linear regression setting.
#'
#' The following steps are applied:
#' \itemize{
#'  \item McScan algorithm, cf. Section 2.1 in Cho, Kley & Li (2024),
#'  \item LOPE estimator \eqn{\hat\delta_j}, cf. Section 3.1 in Cho, Kley & Li (2024),
#'  \item confidence intervals \eqn{C_{ij}(0.1)}, cf. eqn (14),
#'    Section 3.2.2 in Cho, Kley & Li (2024).
#' }
#'
#' @param X Design matrix n x p
#' @param y Response vector n
#' @param ... additional parameters that will be forwarded to \code{[McScan]}
#'
#' @return a named list with elements
#'    \itemize{
#'      \item \code{cp} detected change points, length q
#'      \item \code{delta} n times q matrix, j-th column is \eqn{\delta_j}
#'      \item \code{ci} n times 2 times q array, \code{ci[,1,j]} and
#'        \code{ci[,1,j]} are lower and upper bounds of the confidence intervals,
#'        respectively
#'    }
#' @name inferchange
#' @export
#'
#' @examples
#' old_seed <- .Random.seed
#' set.seed(12345)
#' data <- dgp_gauss_sparse(n = 200, p = 20, z = 100, s = 3, rho = 1, sigma = 1)
#' X <- data$X
#' y <- data$y
#' res <- inferchange(X, y)
#' .Random.seed <- old_seed
inferchange <- function(X, y, ...) {

  n <- nrow(X)
  p <- ncol(X)
  # Step 1: obtain McScan estimates for change locations
  Theta <- McScan(X, y)$cp # , ...

  a <- c()
  b <- c()
  Delta <- c()
  hat_delta <- matrix(nrow = p, ncol = length(Theta))
  colnames(hat_delta) <- 1:length(Theta)
  ci_delta <- array(dim = c(p, 2, length(Theta)))
  for (j in 1:length(Theta)) {
    # compute a_j, b_j and Delta_j from eq (9)
    thjm1 <- c(0, Theta, n)[j]
    thj   <- c(0, Theta, n)[j+1]
    thjp1 <- c(0, Theta, n)[j+2]
    Delta[j] <- min(thj - floor(2 * thjm1 / 3 + thj / 3),
                    ceiling(thj / 3 + 2 * thjp1 / 3) - thj)
    a[j] <- thj - Delta[j]
    b[j] <- thj + Delta[j]

    # estimate differential parameter
    X0 <- X[(a[j]+1):b[j],]
    y0 <- y[(a[j]+1):b[j]]
    k0 <- as.integer(Theta[j] - a[j])
    hat_delta[, j] <- lope(X0, y0, k0)
    ci_delta[, , j] <- ci_delta(X = X0, y = y0, k = k0)$ci
  }
  ret <- list(cp    = Theta,
              delta = hat_delta,
              ci = ci_delta)
  class(ret) <- "inferchange"
  return(ret)

}

#' @importFrom utils head
#' @export
print.inferchange <- function(x, ...) {
  cat("Detected", length(x$cp), "change points at", x$cp, "\n")
  for(j in 1:length(x$cp)) {
    cat("First entries for LOPE estimate and CI for change point ", j, "\n")
    M <- cbind(x$delta[,j], x$ci[,,j])
    colnames(M) <- c("LOPE", "CI-lower", "CI-upper")
    print(head(M))
  }
}
