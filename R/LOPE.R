
#' Compute LOPE estimator
#'
#' Computes the LOPE estimator \eqn{\hat\delta_{0,n}(k)}
#' defined in eqn. (7) of Section 3.1 in Cho et al. (2024).
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
#' @param nfolds number of folds to use in cross-validation;
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
#' @importFrom glmnet glmnet cv.glmnet
#'
#' @examples
#' old_seed <- .Random.seed
#' set.seed(12345)
#' data <- dgp_gauss_sparse(n = 20, p = 20, z = 10, s = 3, rho = 1, sigma = 1)
#' X <- data$X
#' y <- data$y
#' res <- lope(X, y, 10)
#' .Random.seed <- old_seed
lope <- function(X, y, k, standardize = FALSE,
                 lambdapath = NULL, nfolds = 5, nlambda = 100) {

  # Check inputs
  stopifnot(is.matrix(X))
  n = nrow(X)
  p = ncol(X)
  if (length(y) != n) { stop("Input X should be of dim n x p, and y of n!") }
  stopifnot(is.numeric(k) && (k == round(k)) && k >= 1 && k <= n)

  # Pre-computation
  if (standardize) { X = sweep(X, 2, sqrt(colSums(X^2)), "/") }

  # proper scaling
  y[1:k] = y[1:k]*(n/k)
  X[1:k, ] = -X[1:k, ]
  y[(k+1):n] = y[(k+1):n]*(n/(n-k))
  X[(k+1):n, ] = X[(k+1):n, ]

  if (nfolds == 0) {
    # compute solutions from all of the data
    lasso1.all = glmnet::glmnet(X, y, alpha = 1, standardize = FALSE, intercept = FALSE,
                        lambda = lambdapath)
    delta1.all <- as.matrix(lasso1.all$beta)
    colnames(delta1.all) <- sort(lambdapath, decreasing = TRUE)
  } else {
    cv.mod = glmnet::cv.glmnet(X, y, alpha = 1, standardize = FALSE, intercept = FALSE,
                       foldid = rep(1:nfolds, length.out = n))
    lasso.mod = glmnet::glmnet(X, y, alpha = 1, standardize = TRUE, intercept = FALSE,
                       lambda = cv.mod$lambda.min)
    delta1.all <- matrix(lasso.mod$beta, ncol = 1)
    colnames(delta1.all) <- cv.mod$lambda.min

  }

  return(delta1.all)
}
