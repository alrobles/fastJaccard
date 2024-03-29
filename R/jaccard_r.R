#' Compute a Jaccard/Tanimoto similarity coefficient
#'
#' @param x a binary vector (e.g., fingerprint)
#' @param y a binary vector (e.g., fingerprint)
#' @param center whether to center the Jaccard/Tanimoto coefficient by its expectation
#' @param px probability of successes in \code{x} (optional)
#' @param py probability of successes in \code{y} (optional)
#'
#' @return \code{jaccard.test.bootstrap} returns an expected value.
#'
#' @export jaccard_r
#'
#' @examples
#' set.seed(1234)
#' x = rbinom(100,1,.5)
#' y = rbinom(100,1,.5)
#' jaccard_r(x,y)
jaccard_r <- function(x, y, center=FALSE, px=NULL, py=NULL) {
  if(length(x) != length(y)) {
    stop("Two fingerprints (x and y) must be of the same length.")
  }
  
  if(is.null(px) | is.null(py)){
    px <- mean(x)
    py <- mean(y)
  }
  
  sumxy <- sum(x & y)
  unionxy <- sum(x)+sum(y)-sumxy
  if(unionxy == 0) {
    j <- (px*py)/(px+py-px*py)
  } else {
    j <- sumxy/unionxy
  }
  if(center == FALSE) {
    return(j)
  } else {
    return(j - (px*py)/(px+py-px*py))
  }
}