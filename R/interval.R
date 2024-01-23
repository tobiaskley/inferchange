#' Interval Function in Fortran
#'
#' @param n Interval length > 2
#' @param dec Decay rate: Rate of how fast the different layers will decrease in size. Default is set to sqrt(2).
#' @param minl Minimal theoretical interval-length: The minimal size of intervals to be considered. Default is set to 2.
#' @return A 2-column matrix with all the start- and endpoints of the generated intervals
#' @examples
#' \donttest{
#' interval(10)
#' }
#' @export
#' @useDynLib inferchange


interval <-function(n,dec=sqrt(2),minl=2){
  if (!is.integer(n)) {storage.mode(n) <- 'integer'}
  if(n < 3){stop("n should be at least 3")}
  if (!is.double(dec)) {storage.mode(dec) <- 'double'}
  dep <- floor(log(n)/log(dec))
  ilen <- rep(0,dep-1)
  nint <- rep(0L,dep-1)
  for(i in 1:(dep-1)){
    ilen[i] <- n*((1/dec)**i)
    if(ilen[i] < (minl-1)){ dep <- dep-1 }
    else{ nint[i] <- ceiling(n/ilen[i]-1e-8)*2-1 }
  }
  dep <- as.integer(dep)
  nsum <- as.integer(1+sum(nint))
  nint <- as.integer(nint[1:dep])
  ilen <- as.double(ilen[1:dep])
  return(.Call("intervalC",n,dec,dep,ilen,nint,nsum))
}
