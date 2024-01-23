# Logs:
# Housen, April 10, 2023, created

# if (!require("ChangePoints", quietly = TRUE)) {
#   devtools::install_github('kovacssolt/ChangePoints')
# }

#' Bootstrap of covariance scanning stat for a single change point
#'
#' TODO: Description of Details
#'
#' @param X TODO Add Description
#' @param y TODO Add Description
#' @param nb TODO Add Description
#' @param st TODO Add Description
#' @param ed TODO Add Description
#' @param bnd TODO Add Description
#'
#' @return TODO Add Description
#' @export
#'
#' @examples
#' # TODO Add Example or remove this
bts_cov_stat <- function(X, y, nb = max(dim(X)), st = 0, ed = length(y), bnd = 5) {
  n = nrow(X)
  p = ncol(X)
  k = (st + bnd):(ed - bnd)
  kl = sqrt((ed - k)/(ed - st)/(k - st))
  kr = sqrt((k - st)/(ed - st)/(ed - k))
  stat = matrix(nrow = nb, ncol = p)
  for (r in 1:nb) {
    # if (r %% 1e3 == 0) {
    #   print(sprintf("%d / %d", r, nb))
    # }
    id = sample.int(n, replace = TRUE)
    Xb = X[id,]
    yb = y[id]
    M  = Xb * matrix(as.array(yb), nrow = n, ncol = p)
    CS = apply(M, 2, cumsum)
    cov_diff = matrix(kl + kr, nrow = length(k), ncol = ncol(CS)) * CS[k, ]
    cov_diff = cov_diff - matrix(kr, ncol = 1) %*% CS[ed, ]
    if (st > 0) {
      cov_diff = cov_diff - matrix(kl, ncol = 1) %*% CS[st, ]
    }
    stat[r,] = apply(abs(cov_diff), 2, max)
  }
  return(stat)
}

#' Covariance Scanning
#'
#' TODO: Description of Details
#'
#' @param n TODO Add Description
#' @param nIntv TODO Add Description
#' @param minLen TODO Add Description
#'
#' @return TODO Add Description
#' @export
#'
#' @examples
#' # TODO Add Example or remove this
rnd_interval <- function(n, nIntv, minLen = 1) { # generate random interval
  st = sample.int(n, size = nIntv, replace = TRUE)
  ed = sample.int(n, size = nIntv, replace = TRUE)

  tmp = pmax(st, ed)
  st  = pmin(st, ed)
  ed  = tmp
  ret = data.frame("st" = st, "ed" = ed)
  ret = ret[ret$ed - ret$st >= minLen, ]
  if (nrow(ret) == 0) { ret = data.frame("st" = 0, "ed" = n) }
  return(ret[!duplicated(ret), ])
}


#' Title
#'
#' @param est  TODO Add Description
#' @param cid  TODO Add Description
#' @param CS  TODO Add Description
#' @param agg  TODO Add Description
#' @param alpha  TODO Add Description
#' @param stat.bts  TODO Add Description
#' @param post  TODO Add Description
#' @param visual  TODO Add Description
#'
#' @return TODO Add Description
#' @export
#'
#' @examples
#' # TODO Add Example or remove this
post_process <- function(est, cid, CS, agg, alpha, stat.bts, post, visual) {
  if (cid == 0) { return(data.frame()) }
  est = est[1:cid, ]
  est = est[order(est$cp, decreasing = FALSE), ]
  if (post) {
    spl = est$cp
    spl = c(0, round((spl[1:(cid-1)] + spl[2:cid])/2), nrow(CS))
    for (i in 1:cid) {
      if (spl[i+1] - spl[i] > 2*est$bnd[i]) {
        est$st[i] = spl[i]
        est$ed[i] = spl[i+1]
      }
    }
  }
  if (post || visual) {
    for (i in 1:cid) {
      est[i,] = find_single_cp(CS, est$st[i], est$ed[i], est$bnd[i], agg,
                               est$s[i], alpha, stat.bts, visual)
    }
  }
  return(est)
}

#' Find a single change point on given interval
#'
#' Find a single change point on (st, ed] with 'bnd' apart from boundaries.
#'
#' @param CS TODO Add Description
#' @param st TODO Add Description
#' @param ed TODO Add Description
#' @param bnd TODO Add Description
#' @param agg TODO Add Description
#' @param s TODO Add Description
#' @param alpha TODO Add Description
#' @param stat.bts TODO Add Description
#' @param visual TODO Add Description
#'
#' @return TODO Add Description
#'
#' @examples
#' # TODO Add Example or remove this
find_single_cp <- function(CS, st, ed, bnd, agg, s, alpha, stat.bts,
                           visual = FALSE) {
  stopifnot(ed - st > 2*bnd)
  k  = (st + bnd + 1):(ed - bnd)
  kl = sqrt((ed - k)/(ed - st)/(k - st))
  kr = sqrt((k - st)/(ed - st)/(ed - k))
  cov_diff = matrix(kl + kr, nrow = length(k), ncol = ncol(CS)) * CS[k, ]
  cov_diff = cov_diff - matrix(kr, ncol = 1) %*% CS[ed, ]
  if (st > 0) {
    cov_diff = cov_diff - matrix(kl, ncol = 1) %*% CS[st, ]
  }
  if (agg == "inf") {
    score = apply(abs(cov_diff), 1, max)
  } else if (agg == "l2") {
    score = apply(cov_diff^2, 1, sum)
  } else if (agg == "2s") {
    score = t(apply(t(abs(cov_diff)), 2, sort))
    score = score[, ncol(CS):(ncol(CS)-s+1), drop = FALSE]^2
    score = apply(score, 1, sum)
  } else if (agg == "2sFWER") {
    p = ncol(CS)
    alphaQnt = function(x) { return(quantile(x, 1-alpha/p)) }
    thd.feature = apply(stat.bts, 2, alphaQnt)
    crd.max = apply(abs(cov_diff), 2, max)
    id = which(crd.max >= pmin(thd.feature, max(crd.max)))
    score = sqrt(apply(cov_diff[, id, drop = FALSE]^2, 1, mean))
    s = length(id)
  } else if (agg == "2sFDR") {
    p = ncol(CS)
    crd.max = apply(abs(cov_diff), 2, max)
    pVal = apply(stat.bts > matrix(rep(1,nrow(stat.bts)), ncol = 1) %*%
                   matrix(crd.max, nrow = 1), 2, mean)
    id = order(pVal)
    id = id[1:max(c(1,which(pVal[id] <= alpha * (1:p)/p)))]
    score = sqrt(apply(cov_diff[, id, drop = FALSE]^2, 1, mean))
    s = length(id)
  }
  score[1] <- score[length(k)] <- 0 # Down rate the score on the boundaries
  cp  = which.max(score)
  val = score[cp]
  if (visual) {
    lines(k, score, col = "gray")
    abline(v = k[cp], col = "blue")
  }
  return(data.frame("cp" = k[cp], "val" = val, "st" = st, "ed" = ed,
                    "bnd" = bnd, "s" = s))
}

#' Covariance Scanning with inf norm
#'
#' Multiscale covariance scanning for data segmentation described in
#' Section 2.1 of Cho et al. (2024)
#'
#' @param X TODO Add Description
#' @param y TODO Add Description
#' @param type TODO Add Description
#' @param agg TODO Add Description
#' @param alpha TODO Add Description
#' @param nIntv TODO Add Description
#' @param s TODO Add Description
#' @param ncp TODO Add Description
#' @param thd TODO Add Description
#' @param standardise TODO Add Description
#' @param stat.bts TODO Add Description
#' @param bnd TODO Add Description
#' @param post TODO Add Description
#' @param visual TODO Add Description
#'
#' @return TODO Add Description
#' @export
#'
#' @examples
#' # TODO Add Example or remove this
cp_scan <- function(X, # design matrix n x p
                    y, # response vector n
                    # multiple change point approaches
                    type = c("bs",        # binary segmentation
                             "wbs.rand",  # WBS with random intervals
                             "not.rand",  # NOT with random intervals
                             "wbs.seed",  # WBS with seeded intervals
                             "not.seed"), # NOT with seeded intervals
                    # aggregation approaches
                    agg = c("inf",    # infinity norm
                            "l2",     # 2-norm, i.e., sum of squares
                            "2s",     # (2,s)-norm, with user-specified "s"
                            "2sFWER", # selection of s, via FWER
                            "2sFDR"), # selection of s, via FDR
                    alpha = 0.1,     # error rate (FDR, or FWER)
                    nIntv = 1e4,     # number of random intervals
                    s = 1,           # only used if agg = "2s"
                    ncp = length(y), # stop if "ncp" changes are found
                    thd = sqrt(log(max(dim(X)))), # stopping threshold,
                                                  #   overridden by "ncp"
                    standardise = FALSE, # standardise the columns
                    stat.bts = NULL, # bootstrap statistics
                    bnd = max(round(min(2*log(length(X)), length(y)/2)), 5),
                    # at least "bnd" away from boundary
                    post,            # post processing or not (default: FALSE
                                     #   if estimated # cp = 1, else TRUE)
                    visual = FALSE) { # plot or not (default)
  # Check inputs
  stopifnot(is.matrix(X))
  n = nrow(X)
  p = ncol(X)
  stopifnot(n > 2*bnd)
  stopifnot(ncp > 0)
  if (length(y) != n) { stop("Input X should be of dim n x p, and y of n!") }
  if (ncp < n) { thd = 0 }
  agg  = match.arg(agg)
  type = match.arg(type)
  if (standardise) { X = sweep(X, 2, sqrt(colSums(X^2)), "/") }

  # Pre-computation
  M  = X * matrix(as.array(y), nrow = n, ncol = p)
  CS = apply(M, 2, cumsum)

  if (type == "bs") { # Binary segmentation
    ret = data.frame()
    cid = 0

    intv = data.frame("st" = 0, "ed" = n)
    intv = intv[intv$ed - intv$st > 2*bnd, ]

    if (ncp < n) { # Find the threshold giving 'ncp' change points
      while (cid < ncp && nrow(intv) > 0) {
        tmp = data.frame()
        for (i in 1:nrow(intv)) {
          tmp = rbind(tmp, find_single_cp(CS, intv[i, 1], intv[i, 2],
                                          bnd, agg, s, alpha, stat.bts))
        }
        imax = which.max(tmp$val)
        ret = rbind(ret, tmp[imax, ])
        cid = cid + 1

        intv = intv[-imax, ]
        intv = rbind(intv, data.frame("st" = c(ret$st[cid], ret$cp[cid]),
                                      "ed" = c(ret$cp[cid], ret$ed[cid])))
        intv = intv[intv$ed - intv$st > 2*bnd, ]
      }
    } else { # BS with a fixed 'thd'
      i = 1
      while (i <= nrow(intv)) {
        tmp = find_single_cp(CS, intv[i, 1], intv[i, 2],
                             bnd, agg, s, alpha, stat.bts)
        i = i + 1
        if (tmp$val >= thd) {
          ret = rbind(ret, tmp)
          cid = cid + 1
          intv = rbind(intv, data.frame("st" = c(ret$st[cid], ret$cp[cid]),
                                        "ed" = c(ret$cp[cid], ret$ed[cid])))
          intv = intv[intv$ed - intv$st > 2*bnd, ]
        }
      }
    }


  } else if (type %in% c("not.seed", "not.rand")) {
    if (type == "not.seed") {
      intv = interval(n, minl = 2*bnd+1)
    } else {
      intv = rnd_interval(n, nIntv, minLen = 2*bnd+1)
    }
    ret  = data.frame()
    for (i in 1:nrow(intv)) {
      ret = rbind(ret, find_single_cp(CS, intv[i, 1], intv[i, 2],
                                       bnd, agg, s, alpha, stat.bts))
    }
    ret = ret[order(ret$ed - ret$st, decreasing = FALSE), ]
    if  (ncp < n) { # Find the threshold giving 'ncp' change points
      thds = sort(ret$val, decreasing = TRUE)
      for (i in ncp:length(thds)) {
        thd = thds[i]
        aux = ret[ret$val >= thd, ]
        cid = 0
        while (cid < ncp && cid < nrow(aux)) {
          cid = cid + 1
          if (cid < nrow(aux)) { # Remove the intervals containing the found cp
            tmp = aux[(cid+1):nrow(aux), ]
            tmp = tmp[!(tmp$st < aux$cp[cid] - bnd & tmp$ed >= aux$cp[cid] + bnd), ]
            aux = rbind(aux[1:cid, ], tmp)
          }
        }
        if (cid == ncp) {
          break
        }
      }
      ret = aux
    } else { # NOT with a fixed 'thd'
      ret = ret[ret$val >= thd, ]
      cid = 0
      while (cid < ncp && cid < nrow(ret)) {
        cid = cid + 1
        if (cid < nrow(ret)) { # Remove the intervals containing the found cp
          tmp = ret[(cid+1):nrow(ret), ]
          tmp = tmp[!(tmp$st < ret$cp[cid] - bnd & tmp$ed >= ret$cp[cid] + bnd), ]
          ret = rbind(ret[1:cid, ], tmp)
        }
      }
    }
  } else if (type %in% c("wbs.seed", "wbs.rand")) {
    if (type == "wbs.seed") {
      intv = interval(n, minl = 2*bnd+1)
    } else {
      intv = rnd_interval(n, nIntv, minLen = 2*bnd+1)
    }
    ret  = data.frame()
    for (i in 1:nrow(intv)) {
      ret = rbind(ret, find_single_cp(CS, intv[i, 1], intv[i, 2],
                                      bnd, agg, s, alpha, stat.bts))
    }
    ret = ret[ret$val >= thd, ]
    ret = ret[order(ret$val, decreasing = TRUE), ]
    cid = 0
    while (cid < ncp && cid < nrow(ret)) {
      cid = cid +1
      if (cid < nrow(ret)) { # Remove the intervals containing the found cp
        tmp = ret[(cid+1):nrow(ret), ]
        tmp = tmp[!(tmp$st < ret$cp[cid] - bnd & tmp$ed >= ret$cp[cid] + bnd), ]
        ret = rbind(ret[1:cid, ], tmp)
      }
    }
  }
  if (visual == TRUE && cid > 0) {
    # put visualization here
    plot(NULL, xlim = c(0,n), ylim = range(c(0,ret$val)), xlab = "", ylab = "",
         main = paste("COV", agg, sep = "-"))
  }
  if (missing(post)) { post = ifelse(cid > 1, TRUE, FALSE) }
  ret = post_process(ret, cid, CS, agg, alpha, stat.bts, post, visual)
  return(ret)
}
