#' @title Perform a baseline correction to an AFM F-z curve
#'
#' @description
#' This function performs the baseline correction to an AFM F-z curve
#' using the minimum rupture point, higher index closer to surface
#'
#' @usage AlignForceZ(Z, Force, BR)
#' @param (Z, Force) reversely indexed
#' @param BR minimum detach point in the retract segment of the curve
#' that defines the bulk region
#'
#' @return aligned force and zero point
#' @importFrom stats lm predict
#' @examples
#'
#' @export

AlignForceZ <- function(Z, Force, BR) {
  # set bulk force zero
  fz.lm <- lm(y ~ x, data = data.frame(x = Z[1:BR], y = Force[1:BR]))
  fz.pred <- predict(fz.lm, data.frame(x = Z))
  Force <- Force - fz.pred
  
  # find zero force point
  nForce <- which(Force < 0)
  if ((length(nForce) == 0) || (length(nForce) == length(Force)))
  {
    jZP <- c()
  } else {
    i0 <- max(nForce)
    i1 <- i0 + 1
    ZO <- Z[i0] - Force[i0] * (Z[i1] - Z[i0]) / (Force[i1] - Force[i0])
    # slight adjustment to fit experimental data point
    # zero point
    jZP <- which.min(abs(Z - ZO))
  }
  # return adjusted force, zero force point
  return(list(Force = Force, jZP = jZP))
}
