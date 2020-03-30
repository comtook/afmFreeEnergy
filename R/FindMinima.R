#' @title Locate all minima on an AFM F-Z curve
#'
#' @description
#' This function locates all minima on an AFM F-Z curve
#'
#' @usage FindBPoints(Force, jDP, width)
#' @param Force higher index closer to surface
#' @param jDP discontinuous points from function FindDPoints
#' @param width radius of a window
#'
#' @return jBP vector of minimum points
#'
#' @export

FindMinima <- function(Force, jDP, width) {
  jDP <- sort(unique(jDP))
  lngth.jDP <- length(jDP)
  n <- length(Force)
  jDP <- c(jDP, n)
  jBP <- c()
  for (i in 1:lngth.jDP) {
    DP <- jDP[i]
    BR <- min(DP + 1 + width, jDP[i + 1] - 1)
    tRP <- which.min(Force[DP:BR]) + DP - 1
    jBP <- c(jBP, tRP)
  }
  return(jBP)
}

