
#' Compute CLOM estimator
#'
#' Computes the CLOM estimator \eqn{\\hat\\delta_{0,n}(k)}
#' defined in eqn. (8) of Section 3.1 in Cho et al. (2024).
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
#' @param nfolds number of folds to use in cross validation;
#'    if 0 then no cross validation is used; non-negative integer
#' @param nlambda number of values on lambdapath for cross validation;
#'    default 100
#'
#' @return a vector of length p with cross-validation estimate or a
#'    p x length(lambdapath) matrix with the j-th column corresponding to the
#'    j-th value in lambdapath
#' @export
#'
#' @importFrom PRIMAL Dantzig_solver
#'
#' @examples
#' # TODO Add Example or remove this
clom <- function( X, y, k, lambdapath = NULL, nfolds = 5, nlambda = 100 ) {

  n <- nrow(X)
  p <- ncol(X)

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
    D <- Dantzig_solver(X, Ytilde, max_it = 1000, lambda_threshold = min(lambdapath))

    delta1.all <- select_beta(lambdapath, D)
    res <- n * delta1.all

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
    res <- n * delta1.all

  }

  return(res)

}

