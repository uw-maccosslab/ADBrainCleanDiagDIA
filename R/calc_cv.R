#' Computes the coefficient of variation of a numeric vector
#'
#' @param x input numeric vector
#'
#' @return A numeric value of the relative standard deviation
#' @export
#'
#' @examples
#' calc_cv(1:10)
calc_cv <- function(x){

  # compute relative standard deviation
  stats::sd(x, na.rm = T)/mean(x, na.rm = T)
}
