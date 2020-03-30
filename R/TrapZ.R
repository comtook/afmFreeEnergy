#' @title Simple Integration
#'
#' @description
#' This function estimates integration
#'
#' @usage TrapZ(x, y)
#' @param x Z coordinate in nm
#' @param y Force values in pN
#'
#' @return work in pNnm
#'
#' @examples
#'
#' @export

TrapZ <- function(x, y) {
  wrk <- 0.0
  if (length(x) > 1) {
    idx <- 2:length(x)
    wrk <- as.double((x[idx] - x[idx - 1]) %*%
                     (y[idx] + y[idx - 1])) * 0.5
  }
  return(wrk)
}
