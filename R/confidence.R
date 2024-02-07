#' Simultaneous confidence intervals for the differential parameter
#'
#' Generate a de-sparsified estimator of the differential parameter
#' representing the change in the regression coefficients before and
#' after a change point and, based on a Gaussian approximation result,
#' produces a simultaneous confidence interval at a given level.
#'
#' See Cho, Kley & Li (2024) for further details.
#'
#' @param X design matrix with the rows containing the observations
#' @param y vector of the responses
#' @param k index of a single change point; must be an integer between \code{1} and \code{nrow(X) - 1}
#' @param standardize boolean; if \code{standardize = TRUE}, each column of \code{X} is divided by its L2-norm
#' @param delta.hat an estimator of the differential parameter; if \code{delta.hat = NULL}, it is generated via \link[inferchange]{lope} with a tuning parameter generated via cross validation
#' @param nfolds number of folds for the cross validation when producing an estimator of the precision matrix of \code{X}
#' @param nlambdas size of the grid of the tuning parameter for the cross validation when producing an estimator of the precision matrix of \code{X}
#' @param alpha 1 - confidence level between 0 and 1
#' @param M number of random samples to be generated from the target Gaussian distribution
#' @param do.split boolean; whether to split the sample in confidence interval generation or not. When the sample size is moderate, it is recommended not to split the sample.
#' @return an S3 object of class \code{inferchange.ci}, which contains the following fields:
#' \item{delta.check}{de-sparsified estimator of the differential parameter}
#' \item{delta.hat}{original, biased estimator of the differential parameter}
#' \item{ci}{matrix containing the lower and the upper bound of the simultaneous confidence intervals}
#' \item{omega}{estimator of the precision matrix of \code{X}}
#' \item{alpha}{input argument}
#' @examples
#' \donttest{
#' old_seed <- .Random.seed
#' set.seed(1234)
#' data <- dgp_gauss_sparse(n = 300, p = 50, z = 75, s = 5, rho = 2, sigma = 1)
#' X <- data$X
#' y <- data$y
#' ci <- ci_delta(X, y, 75)
#' print(ci)
#' plot(ci)
#' }
#' @importFrom stats cov quantile
#' @export
ci_delta <- function(X, y, k, standardize = FALSE,
                     delta.hat = NULL, nfolds = 3, nlambdas = 50,
                     alpha = .1,  M = 999, do.split = FALSE){

  stopifnot(is.matrix(X))
  n <- dim(X)[1]; p <- dim(X)[2]

  if(length(y) != n) { stop("Input X should be a matrix of dimensions n x p, and y a vector of length n!") }
  stopifnot(is.numeric(k) && (k == round(k)))
  if(!do.split && (k <= 0 || k >= n)) { stop('The change point location k should be between 1 and n - 1') }
  if(do.split && (k <= 1 || k >= n - 1)) { stop('The change point location k should be between 2 and n - 2 if do.split = TRUE') }
  if(alpha < 0 || alpha > 1) { stop('The confidence level alpha should be between 0 and 1') }

  if(standardize) X <- sweep(X, 2, sqrt(colSums(X^2)), "/")
  # if(intercept)

  if(do.split){
    ind0 <- which(1:n %% 2 == 0)
    ind1 <- setdiff(1:n, ind0)
  } else ind0 <- ind1 <- 1:n

  omega <- omega_est(X[ind0,, drop = FALSE], nfolds = nfolds, nlambdas = nlambdas)

  if(is.null(delta.hat)) delta.hat <- lope(X[ind0,, drop = FALSE], y[ind0], which.min(abs(ind0 - k)), lambdapath = NULL, nfolds = nfolds)

  Sigma <- t(X[ind1,, drop = FALSE]) %*% X[ind1,, drop = FALSE] / length(ind1)
  gamma_l <- c(y[ind1[ind1 <= k]] %*% X[ind1[ind1 <= k],, drop = FALSE] / sum(ind1 <= k))
  gamma_r <- c(y[ind1[ind1 > k]] %*% X[ind1[ind1 > k],, drop = FALSE] / sum(ind1 > k))

  delta.check <- delta.hat - omega %*% (Sigma %*% delta.hat - gamma_r + gamma_l)

  u <- X[ind1, ] * 0
  u[ind1 <= k, ] <- X[ind1[ind1 <= k], ] * c(y[ind1[ind1 <= k]] + k/n * X[ind1[ind1 <= k], ] %*% delta.hat)
  u[ind1 > k, ] <- X[ind1[ind1 > k], ] * c(y[ind1[ind1 > k]] - (n - k)/n * X[ind1[ind1 > k], ] %*% delta.hat)
  gamma <- k/n * cov(u[ind1 <= k, ]) + (n - k)/n * cov(u[ind1 > k, ])
  sv <- svd(gamma, nv = 0)
  v.cov <- omega %*% sv$u %*% diag(sqrt(pmax(sv$d, 0)))

  tmp <- v.cov %*% matrix(stats::rnorm(M * p), nrow = p)
  qq <- c(apply(abs(tmp), 2, max), max(abs(delta.check)))
  rr <- quantile(qq, 1 - alpha/2) * sqrt(n/k/(n - k))
  ci <- cbind(delta.check - rr, delta.check + rr)
  colnames(ci) <- c('lower', 'upper')

  out <- list(delta.check = delta.check, delta.hat = delta.hat, ci = ci, omega = omega, alpha = alpha)
  attr(out, "class") <- "inferchange.ci"
  return(out)

}

#' @keywords internal
omega_est <- function(x, nfolds = 3, nlambdas = 50){

  n <- dim(x)[1]
  p <- dim(x)[2]

  xx <- t(x)%*%x/n
  lambda_seq <- exp(seq(log(max(abs(xx))), log(max(abs(xx)) * .01), length.out = nlambdas))
  ind <- sample(1:n, n)

  cv_err <- matrix(0, nrow = p, ncol = nlambdas)
  for(kk in 1:nfolds){
    ind1 <- ind[1:n %% nfolds == kk - 1]
    ind0 <- setdiff(1:n, ind1)
    S <- t(x[ind1, ])%*%x[ind1, ]/length(ind1)
    res <- psm_for_clime_cv(x[ind0, ], S = S, nlambdas = nlambdas, min.lambda = min(lambda_seq))
    for(ii in 1:p){
      for(ll in 1:nlambdas){
        if(lambda_seq[ll] < min(res[ii, , 1])) jj <- nlambdas else jj <- min(which(res[ii,, 1] <= lambda_seq[ll]))
        cv_err[ii, ll] <- cv_err[ii, ll] + res[ii, jj, 2]
      }
    }
  }
  ind <- apply(cv_err, 1, which.min)
  omega <- psm_for_clime(x, nlambdas = nlambdas, min.lambda = min(lambda_seq), lambda = lambda_seq[ind])

  return(omega)

}

#' @importFrom PRIMAL PSM_solver
#' @keywords internal
psm_for_clime <- function(x, nlambdas = 50, min.lambda, lambda){

  n <- dim(x)[1]; p <- dim(x)[2]
  xx <- t(x) %*% x / n

  zeros <- rep(0, 2 * p)
  A <- cbind(cbind(rbind(xx, -xx), - rbind(xx, -xx)), diag(rep(1, 2 * p)))
  c <- c(rep(-1, 2 * p), rep(0, 2 * p))
  c_bar <- rep(0, 4 * p)
  b_bar <- rep(1, 2 * p)
  B_init <- seq(2 * p, 4 * p - 1)

  omega <- matrix(0, p, p)
  for(ii in 1:p){
    b <- zeros
    b[c(ii, ii + p)] <- c(1, -1)

    out <- PRIMAL::PSM_solver(A, b, b_bar, c, c_bar, B_init, max_it = nlambdas, lambda_threshold = min.lambda)

    if(lambda[ii] < min(out$lambda)) jj <- length(out$lambda) else jj <- min(which(out$lambda <= lambda[ii]))
    omega[ii, ] <- out$beta[1:p, jj] - out$beta[p + 1:p, jj]
  }
  omega

}

#' @importFrom PRIMAL PSM_solver
#' @keywords internal
psm_for_clime_cv <- function(x, S, nlambdas = 50, min.lambda){

  n <- dim(x)[1]; p <- dim(x)[2]
  xx <- t(x) %*% x / n

  zeros <- rep(0, 2 * p)
  A <- cbind(cbind(rbind(xx, -xx), - rbind(xx, -xx)), diag(rep(1, 2 * p)))
  c <- c(rep(-1, 2 * p), rep(0, 2 * p))
  c_bar <- rep(0, 4 * p)
  b_bar <- rep(1, 2 * p)
  B_init <- seq(2 * p, 4 * p - 1)

  res <- array(0, dim = c(p, nlambdas, 2))
  for(ii in 1:p){
    b <- zeros
    b[c(ii, ii + p)] <- c(1, -1)

    out <- PRIMAL::PSM_solver(A, b, b_bar, c, c_bar, B_init, max_it = nlambdas, lambda_threshold = min.lambda)
    res[ii, 1:length(out$lambda), 1] <- out$lambda
    res[ii, 1:length(out$lambda), 2] <- apply(t(as.matrix(out$beta[1:p, ] - out$beta[p + 1:p, ])) %*% S, 1, function(zz){ max(abs(zz - b[1:p])) })
    res[ii, out$beta[ii, ] - out$beta[p + ii, ] == 0, 2] <- Inf

  }
  res

}
