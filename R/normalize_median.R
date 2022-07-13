#' Median location Normalization
#'
#' @param df Numeric matrix or dataframe
#'
#' @return Global median location normalized numeric matrix
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(50),5,10)
#' normalize_median(mat)
normalize_median = function (df){
  df_med = apply(as.matrix(df), 2, function(x) median(x, na.rm = T))
  df_med_loc = df_med - median(df_med)
  df.median = t(t(df)-df_med_loc)
  return(df.median)
}
