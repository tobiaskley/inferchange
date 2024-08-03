
#' Compute CLOM estimator
#'
#' Computes the CLOM estimator \eqn{\hat\delta_{0,n}(k)}
#' defined in eqn. (8) of Section 3.1 in Cho, Kley & Li (2024).
#'
#' By default cross-validation is performed and the cross-validation estimate
#' is returned. for cross-validation the estimates are obtained for all
#' values of lambda on the lambdapath obtained with [lambdapath()].
#'
#' If \code{nfolds = 0} is set then the argument \code{lambdapath} has to be
#' a vector of positive numerics and the estimates for these values of lambda
#' will be returned. In preprocessing, \code{lambdapath} is sorted to be in
#' decreasing order.
#'
#' @param X covariates; n x p matrix
#' @param y response; n vector
#' @param k index that splits the n observations into {1, ..., k} and
#'    {k+1, ..., n}; k takes values in {1, ..., n-1}
#' @param standardize boolean; if \code{standardize = TRUE}, each column of
#'    \code{X} is divided by its L2-norm
#' @param lambdapath A user supplied lambda sequence,
#'    only used when no cross validation is done / ignored when nfolds > 1;
#'    vector of non-negative numerics
#' @param nfolds number of folds to use in cross validation;
#'    if 0 then no cross validation is used; non-negative integer
#' @param nlambda number of values on lambdapath for cross validation;
#'    default 100
#'
#' @return a matrix with p rows and with one (in case of cross-validation) or
#'    length(lambdapath) colummns;
#'    each column corresponds to the estimate obtained for one value of lambda
#'    the values of lambda are available as colnames.
#' @export
#'
#' @importFrom PRIMAL Dantzig_solver
#'
#' @examples
#' # TODO Add Example or remove this
clom <- function( X, y, k, standardize = FALSE,
                  lambdapath = NULL, nfolds = 5, nlambda = 100 ) {

  # Check inputs
  stopifnot(is.matrix(X))
  n = nrow(X)
  p = ncol(X)
  if (length(y) != n) { stop("Input X should be of dim n x p, and y of n!") }
  stopifnot(is.numeric(k) && (k == round(k)) && k >= 1 && k <= n)

  # Pre-computation
  if (standardize) { X = sweep(X, 2, sqrt(colSums(X^2)), "/") }

  # helper function to pic betas from Dantzig solver solution path
  select_beta <- function(lambda0, D) {
    if (lambda0 > max(D$lambda)) {
      res <- rep(0, p)
    } else {
      res <- D$beta[ , max(which(D$lambda >= lambda0))]
    }
  }
  select_beta <- Vectorize(select_beta, vectorize.args = "lambda0")

  Ytilde <- c(rep(-1/k, k), rep(1/(n-k), n-k)) * y
  if (nfolds == 0) {
    # compute solutions from all of the data
    D <- PRIMAL::Dantzig_solver(X, Ytilde, max_it = 1000,
                        lambda_threshold = min(lambdapath) / sqrt(k*(n-k)/n))

    delta1.all <- select_beta(lambdapath / sqrt(k*(n-k)/n), D)
    res <- n * delta1.all
    colnames(res) <- lambdapath

  } else {

    # determine lambda_max
    lambdapath0 <- lambdapath(max(abs(t(X) %*% Ytilde)), n, p, nlambda)

    y0 <- y
    X0 <- X
    y0[1:k] = y[1:k]*(n/k)
    X0[1:k, ] = -X[1:k, ]
    y0[(k+1):n] = y[(k+1):n]*(n/(n-k))
    X0[(k+1):n, ] = X[(k+1):n, ]

    cv_err <- rep(0, nlambda)
    for(kk in 1:nfolds){
      ind1 <- (1:n)[0:(n-1) %% nfolds == kk - 1]
      ind0 <- setdiff(1:n, ind1)
      D <- PRIMAL::Dantzig_solver(X[ind0, ], Ytilde[ind0],
                          max_it = 1000, lambda_threshold = min(lambdapath0))
      B <- length(ind0) * select_beta(lambdapath0, D)
      cv_err <- cv_err +
        colSums( ( X0[ind1, ] %*% B - y0[ind1] )^2 )
    }
    lambda_CV <- lambdapath0[which.min(cv_err)]

    D <- PRIMAL::Dantzig_solver(X, Ytilde,
                        max_it = 1000, lambda_threshold = min(lambdapath0))

    delta1.all <- select_beta(lambda_CV, D)
    res <- matrix(n * delta1.all, ncol = 1)
    colnames(res) <- lambda_CV * sqrt(k*(n-k)/n)

  }

  return(res)

}

