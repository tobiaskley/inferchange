
#' Compute LOPE estimator
#'
#' Computes the LOPE estimator \eqn{\\hat\\delta_{0,n}(k)}
#' defined in eqn. (7) of Section 3.1 in Cho et al. (2024).
#'
#' By default cross-validation is performed and the cross-validation estimate
#' is returned. for cross-validation the estimates are obtained for all
#' values of lambda on the lambdapath obtained with [lambdapath()].
#'
#' If \code{nfolds = 0} is set then the argument \code{lambdapath} has to be
#' a vector of positive numerics and the estimates for these values of lambda
#' will be returned.
#'
#' @param X covariates; n x p matrix
#' @param y response; n vector
#' @param k index that splits the n observations into {1, ..., k} and
#'    {k+1, ..., n}; k takes values in {1, ..., n-1}
#' @param lambdapath A user supplied lambda sequence,
#'    only used when no cross validation is done / ignored when nfolds > 1;
#'    vector of non-negative numerics
#' @param nfolds number of folds to use in cross-validation;
#'    if 0 then no cross validation is used; non-negative integer
#' @param nlambda number of values on lambdapath for cross validation;
#'    default 100
#'
#' @return a vector of length p with cross-validation estimate or a
#'    p x length(lambdapath) matrix with the j-th column corresponding to the
#'    j-th largest value in lambdapath
#' @export
#'
#' @importFrom glmnet glmnet cv.glmnet
#'
#' @examples
#' # TODO Add Example or remove this
lope <- function(X, y, k, lambdapath = NULL, nfolds = 5, nlambda = 100) {

  n <- nrow(X)
  p <- ncol(X)

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
  } else {
    cv.mod = glmnet::cv.glmnet(X, y, alpha = 1, standardize = FALSE, intercept = FALSE,
                       foldid = rep(1:nfolds, length.out = n))
    lasso.mod = glmnet::glmnet(X, y, alpha = 1, standardize = TRUE, intercept = FALSE,
                       lambda = cv.mod$lambda.min)
    delta1.all <- as.numeric(lasso.mod$beta)
  }

  return(delta1.all)
}
