#' @title Locate quasi-equilibrium points on an AFM F-Z curve
#'
#' @description
#' This function locates quasi-equilibrium points before pulling
#' based on zero force criterion within error bar
#'
#' @usage FindQPoints(Z, Force, jBP, noise)
#' @param (Z, Force) retraction force high index close to surface
#' @param jBP index corresponding to significant force minima
#' @param noise lower bound of baseline
#'
#' @return jQP vector of quasipoints
#'
#' @export
#'

FindQPoints <- function(Force, jBP, noise) {
  # detect quasi-equilibrium points
  jQP <- c()
  if (length(jBP) < 2) {
    return(jQP)
  } else {
    length.jBP <- length(jBP) - 1
  }
  
  for (i in 1:length.jBP) {
    jL <- jBP[i]
    jR <- jBP[i + 1] - 1
    jRForce <- abs(Force[jL:jR])
    # quasi-points
    jMin <- min(which.min(jRForce)) + jL - 1

    if (abs(Force[jMin]) < noise) {
      jQP <- c(jQP, jMin)
    }
  }
  
  # add contact point
  if (length(jQP) > 0) {
    jL <- max(jQP)
    jR <- max(jBP)
    tmp <- which(abs(Force[jL:jR]) < noise) + jL - 1
    jQP <- c(jQP, tmp)
  }
  
  return(unique(jQP))
}

