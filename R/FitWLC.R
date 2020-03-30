#' @title Fit into Worm-Like Chain Model
#'
#' @description
#' This function fits force-separation curves into worm-like chain
#'
#' @usage FitWLC(myData, nlsEr = 4)
#' @param myData stored jBP, jQP, Z, Q and Force
#' @param nlsEr tolerence for fitting error
#'
#' @return iBP, iQP, iCP, WLC
#'
#' @importFrom stats nls predict
#'
#' @examples
#'
#' @export

FitWLC <- function(myData, nlsEr = 3) {
  Z <- myData$Z
  Force <- myData$Force
  sigma <- myData$sigma
  jBP <- myData$jBP
  jQP <- myData$jQP
  iZP <- myData$jZP
  width <- myData$width
  
  noise <- sigma * nlsEr          # tolerence for nls error
  
  # rupture points
  iRP <- max(jBP)
  jCP <- jQP[(jQP - iRP) < 0]

  if (length(jCP) == 0) {
    # no contact point
    return(myData)
  } else {
    iQP <- max(jCP)
  }

  stretch <- 11
  if ((iRP - iQP) < stretch) {
    # impose minimum stretch 11 steps
    return(myData)
  } else {
    iRP <- which.min(Force[iQP:iRP]) + iQP - 1
  }

  # fit into worm-like model
  selZ <- Z[iZP:iRP]
  tmpForce <- Force - Force[iQP]
  selForce <- c(rep(0, iQP - iZP + 1), tmpForce[(iQP + 1):iRP])
  dat <- data.frame(x = selZ - Z[iZP], y = selForce)
  # theoretical persistence length
  startLp <- 0.37   # Anna Rita 2011
  # allow 0.1 error
  minLp <- 0.9 * startLp
  maxLp <- 1.1 * startLp
  # extension
  minL <- Z[iRP] - Z[iZP]
  maxL <- 1.20 * minL
  startL <- Z[iRP + 1] - Z[iZP]
  # original WLC model
  fitWLC <- tryCatch(nls(y ~ I(-(4.143 / Lp) * ((0.25 / ((1.0 - (x / L))^2)) -
                                               0.25 + (x / L))),
                data = dat, start = list(Lp = startLp, L = startL),
                lower = c(minLp, minL), upper = c(maxLp, maxL),
                algorithm = 'port'),
                error = function(e) e)
  if (!any(class(fitWLC) == "error")) {
    if ((summary(fitWLC)$sigma < noise)) {
      Lp <- summary(fitWLC)$coefficients[1, 1]
      Lcnt <- summary(fitWLC)$coefficients[2, 1]
      predZ <- selZ - Z[iZP]
      predWLC <- predict(fitWLC, data.frame(x = predZ))
      WLC <- data.frame(x = predZ + Z[iZP], y = predWLC + Force[iQP])
      myData$Label <- 1
      myData[["Lp"]] <- Lp
      myData[["L"]] <- Lcnt
      myData[["iRP"]] <- iRP             # rupture point
      myData[["iQP"]] <- iQP             # quasi point before pulling
      myData[["iCP"]] <- iZP             # contact point
      myData[["WLC"]] <- WLC             # fitted curve
    }
  }
  return(myData)
}
