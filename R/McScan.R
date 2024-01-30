# Logs:
# Housen, April 10, 2023, created
# Housen, Jan 30, 2024, simplified codes

#' Post-processing
#' @keywords internal
.post_process <- function(est, cid, CS, post) {
  if (cid == 0) { return(data.frame()) }
  est = est[1:cid, ]
  est = est[order(est$cp, decreasing = FALSE), ]
  if (post) {
    spl = c(0, est$cp, nrow(CS))
    for (i in 1:cid) {
      if (spl[i+2] - spl[i] > 2*est$bnd[i]) {
        est$st[i] = spl[i]
        est$ed[i] = spl[i+2]
        est$bnd[i] = max(est$bnd[i],
                         floor(min(spl[i+1]-spl[i], spl[i+2]-spl[i+1])/2-1))
      }
    }
    for (i in 1:cid) {
      est[i,] = .find_single_cp(CS, est$st[i], est$ed[i], est$bnd[i])
    }
  }
  return(est)
}


#' Find a single change point on (st, ed] with 'bnd' apart from boundaries.
#' @keywords internal
.find_single_cp <- function(CS, st, ed, bnd) {
  stopifnot(ed - st > 2*bnd)
  k  = (st + bnd + 1):(ed - bnd)
  kl = sqrt((ed - k)/(ed - st)/(k - st))
  kr = sqrt((k - st)/(ed - st)/(ed - k))
  cov_diff = matrix(kl + kr, nrow = length(k), ncol = ncol(CS)) * CS[k, ]
  cov_diff = cov_diff - matrix(kr, ncol = 1) %*% CS[ed, ]
  if (st > 0) {
    cov_diff = cov_diff - matrix(kl, ncol = 1) %*% CS[st, ]
  }
  score = apply(abs(cov_diff), 1, max) # default: "inf"
  score[1] <- score[length(k)] <- 0 # Down rate the score on the boundaries
  cp  = which.max(score)
  val = score[cp]
  return(data.frame("cp" = k[cp], "val" = val,
                    "st" = st, "ed" = ed, "bnd" = bnd))
}


#' @title Multiscale covariance scanning for data segmentation described in
#' Section 2.1 of Cho et al. (2024)
#'
#' @param X Design matrix n x p
#' @param y Response vector n
#' @param ncp Specifies the number of estimated change points, and
#'            overrides the input \code{thd}
#' @param thd Stopping threshold, default log(n*p)/2
#' @param method \code{"not"}  NOT (default) with fixed threshold \code{thd}
#'               \code{"auto"} NOT with automatic threshold (overrides \code{ncp} and \code{thd})
#'               \code{"wbs"}  wild binary segmentation with seeded intervals
#'               \code{"bs"}   binary segmentation
#' @param standardise Logical, standardising every coordinate or not
#'                    Default TRUE when \code{thd} is in action
#' @param bnd At least \code{bnd} away from boundaries of intervals
#' @param post Logical, post processing (default) or not
#'
#' @return A list of class \code{inferchange.cp} with \code{print} and \code{plot}
#' @export
#'
#' @importFrom stats mad
#' @examples
#' # TODO Add Example or remove this
McScan <- function(X, y, ncp, thd, method = c("not", "auto", "wbs", "bs"),
                   bnd = max(round(min(2*log(length(X)), length(y)/2)), 5),
                   standardise = FALSE, post = TRUE) {
  # Check inputs
  stopifnot(is.matrix(X))
  n = nrow(X)
  p = ncol(X)
  if (length(y) != n) { stop("Input X should be of dim n x p, and y of n!") }
  stopifnot(n > 2*bnd)
  method = match.arg(method)
  if (method == "auto") {
    if (!missing(ncp)) { warning("Input 'ncp' is ignored!") }
    if (!missing(thd)) { warning("Input 'thd' is ignored!") }
  } else {
    if (missing(ncp)) {
      if (missing(thd)) {
        thd = log(n*p)/2
        cat('Default threshold log(n*p)/2 is used!\n')
      }
      standardise = TRUE
    } else {
      stopifnot(ncp > 0 && ncp == round(ncp))
      if (!missing(thd)) { warning("Input 'thd' is ignored!") }
      thd = 0
      if (ncp == 1) { method = "bs" }
    }
  }

  # Pre-computation
  M  = X * matrix(as.array(y), nrow = n, ncol = p)
  if (standardise) {
    M = sweep(M, 2, apply(M, 2, function(x) {mad(diff(x))/sqrt(2)}), "/")
  }
  CS = apply(M, 2, cumsum)

  ret = data.frame()
  if (method == "bs") { # Binary segmentation
    cid  = 0
    intv = data.frame("st" = 0, "ed" = n)
    intv = intv[intv$ed - intv$st > 2*bnd, ]

    if (missing(ncp)) { # BS with a fixed 'thd'
      i = 1
      while (i <= nrow(intv)) {
        tmp = .find_single_cp(CS, intv[i, 1], intv[i, 2], bnd)
        i = i + 1
        if (tmp$val >= thd) {
          ret = rbind(ret, tmp)
          cid = cid + 1
          intv = rbind(intv, data.frame("st" = c(ret$st[cid], ret$cp[cid]),
                                        "ed" = c(ret$cp[cid], ret$ed[cid])))
          intv = intv[intv$ed - intv$st > 2*bnd, ]
        }
      }
    } else { # Find the threshold giving 'ncp' change points
      while (cid < ncp && nrow(intv) > 0) {
        tmp = data.frame()
        for (i in 1:nrow(intv)) {
          tmp = rbind(tmp, .find_single_cp(CS, intv[i, 1], intv[i, 2], bnd))
        }
        imax = which.max(tmp$val)
        ret = rbind(ret, tmp[imax, ])
        cid = cid + 1

        intv = intv[-imax, ]
        intv = rbind(intv, data.frame("st" = c(ret$st[cid], ret$cp[cid]),
                                      "ed" = c(ret$cp[cid], ret$ed[cid])))
        intv = intv[intv$ed - intv$st > 2*bnd, ]
      }
    }
  } else if (method == "not") {
    intv = seeded_interval(n, minl = 2*bnd+2)
    for (i in 1:nrow(intv)) {
      ret = rbind(ret, .find_single_cp(CS, intv[i, 1], intv[i, 2], bnd))
    }
    ret = ret[order(ret$ed - ret$st, decreasing = FALSE), ]
    if (missing(ncp)) { # NOT with a fixed 'thd'
      ret = ret[ret$val >= thd, ]
      cid = 0
      while (cid < nrow(ret)) {
        cid = cid + 1
        if (cid < nrow(ret)) { # Remove the s containing the found cp
          tmp = ret[(cid+1):nrow(ret), ]
          # tmp = tmp[!(tmp$st < ret$cp[cid] - bnd & tmp$ed >= ret$cp[cid] + bnd), ]
          tmp = tmp[!(tmp$st < ret$cp[cid] & tmp$ed >= ret$cp[cid]), ]
          ret = rbind(ret[1:cid, ], tmp)
        }
      }
    } else { # Find the threshold giving 'ncp' change points
      thds = sort(ret$val, decreasing = TRUE)
      for (i in ncp:length(thds)) {
        thd = thds[i]
        aux = ret[ret$val >= thd, ]
        cid = 0
        while (cid < ncp && cid < nrow(aux)) {
          cid = cid + 1
          if (cid < nrow(aux)) { # Remove the s containing the found cp
            tmp = aux[(cid+1):nrow(aux), ]
            # tmp = tmp[!(tmp$st < aux$cp[cid] - bnd & tmp$ed >= aux$cp[cid] + bnd), ]
            tmp = tmp[!(tmp$st < aux$cp[cid] & tmp$ed >= aux$cp[cid]), ]
            aux = rbind(aux[1:cid, ], tmp)
          }
        }
        if (cid == ncp) { break }
      }
      ret = aux
    }
  } else if (method == "wbs") {
    intv = seeded_interval(n, minl = 2*bnd+2)
    for (i in 1:nrow(intv)) {
      ret = rbind(ret, .find_single_cp(CS, intv[i, 1], intv[i, 2], bnd))
    }
    ret = ret[ret$val >= thd, ]
    ret = ret[order(ret$val, decreasing = TRUE), ]
    cid = 0
    while (cid < ncp && cid < nrow(ret)) {
      cid = cid +1
      if (cid < nrow(ret)) { # Remove the s containing the found cp
        tmp = ret[(cid+1):nrow(ret), ]
        # tmp = tmp[!(tmp$st < ret$cp[cid] - bnd & tmp$ed >= ret$cp[cid] + bnd), ]
        tmp = tmp[!(tmp$st < ret$cp[cid] & tmp$ed >= ret$cp[cid]), ]
        ret = rbind(ret[1:cid, ], tmp)
      }
    }
  } else if (method == "auto") {
    ## compute NOT solution path
    intv = seeded_interval(n, minl = 2*bnd+2)
    losc = data.frame() # record local scan
    for (i in 1:nrow(intv)) {
      losc = rbind(losc, .find_single_cp(CS, intv[i, 1], intv[i, 2], bnd))
    }
    losc = losc[order(losc$ed - losc$st, decreasing = FALSE), ]
    thds = unique(sort(losc$val, decreasing = TRUE))
    sopa = list() # solution path
    for (thd in thds) {
      vid = which(losc$val >= thd)
      rem = rep(TRUE, nrow(losc))
      cps = numeric()
      while (length(vid) > 0) {
        cp  = losc$cp[vid[1]]
        # rem = rem & !(losc$st < cp - bnd & losc$ed >= cp + bnd)
        rem = rem & !(losc$st < cp & losc$ed >= cp)
        vid = which((losc$val >= thd) & rem)
        cps = c(cps, cp)
      }
      sopa = rbind(sopa, list("thresh" = thd,
                              "val" = max(losc$val[rem], -Inf),
                              "cps" = sort(cps)))
    }

    ## Heuristic based model selection
    ncp_ = sapply(sopa[,"cps"], length)
    nc   = unique(ncp_)
    lb   = numeric(length(nc))
    ub   = numeric(length(nc))
    for (i in seq_along(nc)) {
      lb[i] = min(unlist(sopa[ncp_ == nc[i], "val"]))
      ub[i] = max(unlist(sopa[ncp_ == nc[i], "val"]))
    }
    srng = numeric(length(nc))
    srng[1] = ub[1]-min(lb[1],ub[2])
    if (length(nc) > 1) {
      for (i in 2:(length(nc)-1)) {
        srng[i] = max(max(lb[i-1], ub[i]) - min(ub[i+1], lb[i]),
                      -Inf, na.rm = TRUE)
      }
    }
    # plot(nc, srng, type="b")
    i_ = min(which(diff(srng) <= 0)) + 1
    id = which(ncp_ == nc[i_])
    id = id[which.min(unlist(sopa[id, "val"]))]
    ret = data.frame("cp" = unlist(sopa[id, "cps"]),
                     "val" = NA,
                     "st" = NA,
                     "ed" = NA,
                     "bnd" = bnd)
    cid = nrow(ret)
  }
  ret = .post_process(ret, cid, CS, post)
  if (exists("sopa")) { # when solution path is available
    attr(ret, "solution_path") = sopa
  }
  attr(ret, "X") = X
  attr(ret, "y") = y
  class(ret) = "inferchange.cp"
  return(ret)
}
