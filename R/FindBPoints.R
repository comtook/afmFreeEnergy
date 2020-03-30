#' @title Locate significant minima on an AFM F-Z curve
#'
#' @description
#' This function locates significant minima on an AFM F-Z curve
#'
#' @usage FindBPoints(Force, jDP, width, noise)
#' @param Force higher index closer to surface
#' @param jDP discontinuous points from function FindDPoints
#' @param width radius of a window
#' @param noise signal threshold
#'
#' @return jBP vector of binding points
#'
#' @export

FindBPoints <- function(Force, jDP, width, noise) {
  lngth.jDP <- length(jDP)
  n <- length(Force)
  jDP <- c(jDP, n)
  jBP <- c()
  for (i in 1:lngth.jDP) {
    DP <- jDP[i]
    BR <- min(DP + 1 + width, jDP[i + 1] - 1)
    tRP <- which.min(Force[DP:BR]) + DP - 1
    # get rid of dangling points
    if (max(Force[tRP:(tRP + 2)]) < noise) {
      jBP <- c(jBP, tRP)
    }
  }
  return(sort(unique(jBP)))
}

