#' @title Find binding points, zero force points, quasi-equilibrium points, do alignment
#'
#' @description
#' This function locates binding points, zero force points, quasi-equilibrium points
#' and performs the baseline correction on an AFM F-Z curve
#'
#' @usage FindKeyPoints(myData, Alpha, LOD, mulBase)
#' @param Alpha percent of points in bulk
#' @param LOD noise level in sigma unit
#' @param mulBase baseline perturbation level
#'
#' @return jDP, jZP, jBP, jQP and aligned force
#'
#' @importFrom stats lm predict
#' @examples
#'
#' @export

FindKeyPoints <- function(myData, Alpha = 0.5, LOD = 4, mulBase = 2) {
  # locate rupture, binding and contact points
  # initialize
  myData$jDP <- c()           # discontinuous points
  myData$jBP <- c()           # all force peaks location
  myData$jQP <- c()           # all quasi-equilibrium points
  myData$Label <- 0           # 1 indicates successful pulling

  Z <- myData$Z
  Force <- myData$Force
  n <- length(Z)
  Direction <- Z[n] - Z[1]
  # want smaller index corresponding to flat curve
  if ( Direction > 0 ) {
    Z <- rev(Z)
    Force <- rev(Force)
  }

  # find discontinuous points
  DP <- FindDPoints(Z, Force, Alpha, LOD)
  jDP <- DP$jDP
  if ((length(jDP) == 0) || (min(jDP) < (Alpha * n))) {
    # no jDP or jDP in bulk
    return(myData)
  }

  # do alignment
  width <- DP$width
  BR <- min(jDP) - width
  alignFZ <- AlignForceZ(Z, Force, BR)
  Force <- alignFZ$Force
  jZP <- alignFZ$jZP
  sigma <- sd(Force[1:BR])

  # find binding points
  threshold <- -(LOD * sigma)
  jBP <- FindBPoints(Force, jDP, width, threshold)

  # find quasi-equilibrium points
  noise <- mulBase * sigma
  jQP <- FindQPoints(Force, jBP, noise)

  # change back
  if ( Direction > 0 ) {
    Z <- rev(Z)
    Force <- rev(Force)
    jDP <- n - jDP + 1
    jZP <- n - jZP + 1
    jBP <- n - jBP + 1
    jQP <- n - jQP + 1
  }

  myData$width <- width
  myData$jDP <- jDP             # discontinuous points
  myData$jBP <- jBP             # all force peaks location
  myData$jQP <- jQP             # all quasi-equilibrium points
  myData$sigma <- sigma         # standard deviation of force in bulk
  myData$jZP <- jZP             # zero point
  myData$Force <- Force         # aligned force

  return(myData)
}

