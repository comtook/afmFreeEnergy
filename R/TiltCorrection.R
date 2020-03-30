#' @title Perform a baseline correction to an AFM F-Z curve
#'
#' @description
#' This function performs the baseline correction to an AFM F-Z curve
#' setting force at points in bulk to be zero
#'
#' @usage TiltCorrection(Z, Force, Alpha)
#' @param (Z, Force) reversely indexed
#' @param Alpha percent of points in bulk
#'
#'
#' @return aligned force and zero point
#' @importFrom stats lm predict
#' @examples
#'
#' @export

TiltCorrection <- function(Z, Force, Alpha) {
  # set bulk force zero
  n <- length(Force)
  BR <- as.integer(Alpha * n)
  fz.lm <- lm(y ~ x, data = data.frame(x = Z[1:BR], y = Force[1:BR]))
  fz.pred <- predict(fz.lm, data.frame(x = Z))
  Force <- Force - fz.pred
  
  return(list(Force = Force, sigma = sd(Force[1:BR])))
}
