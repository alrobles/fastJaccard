# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rcpp_hello_world <- function() {
    .Call(`_fastJaccard_rcpp_hello_world`)
}

#' Multiplies two doubles
#'
#' @param mat Binary matrix to be evaluated
#' @return Jaccard simmilarity matrix
rcpp_parallel_jaccard_distance <- function(mat) {
    .Call(`_fastJaccard_rcpp_parallel_jaccard_distance`, mat)
}
