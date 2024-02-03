#' Generate a simulated data set
#'
#' @param n number of observations
#' @param p dimension
#' @param z length of the first segment
#' @param s sparsity
#' @param rho L2 norm of delta
#' @param sigma standard deviation of errors
#'
#' @return a list with elements X and y
#' @export
#'
#' @importFrom stats rnorm runif
#'
#' @examples
#' data <- dgp_gauss_sparse(n = 20, p = 20, z = 10, s = 3, rho = 1, sigma = 1)
dgp_gauss_sparse <- function(n, p, z, s, rho = 1, sigma = 1) {

  X  <- matrix(rnorm(n * p), nrow = n, ncol = p)
  mu <- matrix(0, nrow = p)

  Delta <- rep(0, p)

  # generate random indices to change;
  # entries 1 to s of a random permutation of (1 .. p)
  idx <- order(runif(p))[1:s]

  Delta[idx] <- rnorm(s)
  Delta[idx] <- rho * Delta[idx] / sqrt(sum(Delta[idx]^2))

  beta1 <- mu - Delta/2
  beta2 <- mu + Delta/2

  eps <- matrix(rnorm(n, sd = sigma), nrow = n)

  y <- matrix(nrow = n, ncol = 1)
  y[1:z] <- X[1:z, ] %*% beta1 + eps[1:z]
  y[(z+1):n] <- X[(z+1):n, ] %*% beta2 + eps[(z+1):n]

  return(list(X = X, y = as.vector(y)))
}
