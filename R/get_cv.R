#' Perform a row-wise relative standard deviation calculation and return a numeric vector of the coefficient of variation
#'
#' @param mat input data matrix
#'
#' @return A numeric vector of the coefficient of variation
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(50),5,10)
#' get_cv(mat)
get_cv <- function(mat){

  # perform row-wise cv calculation
  cv <- apply(mat, 1, function (x) calc_cv(x))
  return(cv)
}
