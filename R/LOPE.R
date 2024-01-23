
#' Compute LOPE estimator
#'
#' TODO: Description of Details
#'
#' @param X covariates; n x p matrix
#' @param y response; n vector
#' @param k index that splits the n observations into {1, ..., k} and
#'    {k+1, ..., n}; k takes values in {1, ..., n-1}
#' @param lambdapath values of lambda to be used when no cross validation is
#'    is done; vector on non-negative numerics
#' @param nfolds number of folds to use in cross validation;
#'    if 0 then no cross validation is used; non-negative integer
#' @param nlambda number of values on lambdapath for cross validation;
#'    default 100
#'
#' @return TODO Add Description
#' @export
#'
#' @examples
#' # TODO Add Example or remove this
lope <- function(X, y, k, lambdapath, nfolds = 0, nlambda = 100) {

  n <- nrow(X)
  p <- ncol(X)

  # proper scaling
  y[1:k] = y[1:k]*(n/k)
  X[1:k, ] = -X[1:k, ]
  y[(k+1):n] = y[(k+1):n]*(n/(n-k))
  X[(k+1):n, ] = X[(k+1):n, ]

  if (nfolds == 0) {
    # compute solutions from all of the data
    lasso1.all = glmnet(X, y, alpha = 1, standardize = FALSE, intercept = FALSE,
                        lambda = lambdapath)
    delta1.all <- as.matrix(lasso1.all$beta)
  } else {
    cv.mod = cv.glmnet(X, y, alpha = 1, standardize = FALSE, intercept = FALSE,
                       foldid = rep(1:nfolds, length.out = n))
    lasso.mod = glmnet(X, y, alpha = 1, standardize = TRUE, intercept = FALSE,
                       lambda = cv.mod$lambda.min)
    delta1.all <- as.numeric(lasso.mod$beta)
  }

  return(delta1.all)
}
