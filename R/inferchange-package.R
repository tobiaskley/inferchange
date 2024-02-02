#' @keywords internal
"_PACKAGE"

################################################################################
#' Multiscale Covariance Scanning and Related Algorithms
#'
#' Implements function for inference of changes in high-dimensional
#' linear regression; cf. TODO refer to arxiv preprint.
#'
#' @details
#'  \tabular{ll}{
#'    \cr Package: \tab inferchange
#'    \cr Type:    \tab Package
#'    \cr Version: \tab 0.0.0.9000
#'    \cr Date:    \tab 2024-02-??
#'    \cr License: \tab GPL (>= 3)
#'  }
#'
#' @section Contents:
#' The \pkg{inferchange} package contains ...
#'
#' More text.
#'
#' \enumerate{
#'   \item A
#'   \item \code{B}
#' }
#'
#' @section Another section:
#' To improve ...
#'
#' \itemize{
#'  \item What functions to look at (cf. \code{\link{McScan}})
#'  \item Workhorse function,
#'         (cf. \code{\link{inferchange}})
#' }
#'
#' @name inferchange-package
#' @author Haeran Cho, Tobias Kley, Housen Li
#'
#' @references
#' Cho, H., Kley, T., and Li, H. (2024).
#' Detection and inference of changes in high-dimensional linear regression with non-sparse structures
#' (cf. \url{http://arxiv.org/abs/????.?????})
#'
NULL

# Taken from quantreg-package and adapted.
".onAttach" <- function(lib, pkg) {
  if(interactive() || getOption("verbose"))
    packageStartupMessage("Package inferchange loaded.\n     To cite, see citation(\"inferchange\").\n     For demos, see demo(package = \"inferchange\").")
}
