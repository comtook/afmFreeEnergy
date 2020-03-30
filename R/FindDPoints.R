#' @title Locate discontinuous points on an AFM F-Z curve
#'
#' @description
#' This function locates discontinuous points on an AFM F-Z curve
#' discontinuous points correspond to positive peaks of the derivatives
#'
#' @usage FindDPoints(Z, Force, Alpha, LOD)
#' @param (Z, Force) reversely indexed
#'
#' @return jDP vector of discontinuous points
#' @examples
#'
#' @export

FindDPoints <- function(Z, Force, Alpha = 0.5, LOD = 4) {
  # discontinuous points correspond to maximum peaks of derivatives
  tiltC <- TiltCorrection(Z, Force, Alpha)
  Force <- tiltC$Force
  cutoff <- -(LOD * tiltC$sigma)
  n <- length(Force)
  app <- matrix(c(rep(1, n), Z, Force), nrow = n, ncol = 3)
  maxWidth <- as.integer(0.5 * Alpha * n)
  width.try <- seq.int(from = 5, to = maxWidth, 3)
  jBP.old <- n
  jBP <- c()
  for (width in width.try) {
    # calculate force derivative 
    bRoll <- CalcSlope(app, width)
    delta <- c(rep(bRoll[1], width), bRoll, rep(bRoll[length(bRoll)], width))
    # retain the sharpest 
    jDP <- which.max(delta)
    # find maxima of delta
    deltas <- diff(delta,1)
    turns <- which(deltas[-1] * deltas[-length(deltas)] < 0) + 1
    turnsMax <- subset(turns, ((delta[turns] - delta[turns - 1]) > 0) & 
                            ((delta[turns] - delta[turns + 1]) > 0))
    # keep about top 25% candidates
    uppBound <- quantile(bRoll)[4] + 0.1 * IQR(bRoll)
    turnsMax <- subset(turnsMax, delta[turnsMax] > uppBound)
    k <- length(turnsMax)
    if (k > 0) {  
      jBP <- FindMinima(Force, turnsMax, width)
      jBP <- subset(jBP, Force[jBP] < cutoff)
      jBP.new <- length(jBP)
    } else {
      break
    } 
    if (jBP.new == jBP.old) {
      break
    } else {
      jBP.old <- jBP.new
    }
  }
  return(list(jDP = jBP, width = width))
}
