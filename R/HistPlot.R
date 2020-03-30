#' @title Plot histograms of all variables
#'
#' @description This function extracts all parameters and do histogram plot
#'
#' @param outDir Directory is where results.RData locates.
#' @param Temp temperature in Kelvin
#'
#' @export
#'

HistPlot <- function(outDir, Temp, springK) {
  # output selected curves
  Cls <- paste0(outDir, "/Selected")
  if ( dir.exists(Cls) ) {unlink(Cls, recursive <- TRUE)}
  dir.create(Cls)

  # load data
  load(file=file.path(outDir, 'results.RData'))

  # loop all trajectories
  numData <- length(myOutput)
  fName <- names(myOutput)
  rupW <- c()     # work
  rupF <- c()     # force
  rupQ <- c()     # separation
  shape <- c()    # parabolic length
  Sigma <- c()
  selectedName <- c()
  for ( j in 1:numData ) {
    myData <- myOutput[[j]]
    if (myData$Label == 1) {
      Z <- myData$Z
      Force <- myData$Force
      Q <- Z + Force / springK      # separation distance
      iQP <- myData$iQP
      iRP <- myData$iRP
      Qr <- Q[iRP]
      Qc <- Q[iQP]
      Wr <- TrapZ(Z[iQP:iRP], Force[iQP] - Force[iQP:iRP])
      rupW = c(rupW, Wr)
      rupF <- c(rupF, Force[iRP])
      rupQ <- c(rupQ, Qr)
      shape <- c(shape, Qr - Qc)
      Sigma <- c(Sigma, myData$sigma)
      selectedName = c(selectedName, fName[j])
      PlotQFWLC(myData, Cls, fName[j])
    }
  }

  if (length(rupW) > 10) {
  # convert pNnm to kT
  kT <- Temp * 0.00198709     # kcal/mol
  # nm * pN = 0.143953 kcal/mol
  nmpN <- 0.143953
  betA <- nmpN / kT
  rupW <- betA * rupW
  x <- data.frame(names=selectedName, pN=rupF, nm=rupQ, kT=rupW, sigma=Sigma)
  save(x, file = file.path(outDir, 'work.RData'))
  # plot histogram
  jpeg(file.path(outDir, "HistWork.jpeg"),
       units = "in", width = 3.5, height = 3.5, quality = 100, res = 300)
  mar.default <- c(5, 4, 4, 2) + 0.1
  par(oma=c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  hist(rupW, main = "", xlab = expression(paste("Work (k"[B], "T)")), ylab = "Frequency")
  dev.off()

  jpeg(file.path(outDir, "HistForce.jpeg"),
       units = "in", width = 3.5, height = 3.5, quality = 100, res = 300)
  mar.default <- c(5, 4, 4, 2) + 0.1
  par(oma=c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  hist(rupF, main = "", xlab = "Force (pN)", ylab = "Frequency")
  dev.off()

  jpeg(file.path(outDir, "HistQ.jpeg"),
       units = "in", width = 3.5, height = 3.5, quality = 100, res = 300)
  mar.default <- c(5, 4, 4, 2) + 0.1
  par(oma=c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  hist(rupQ, main = "", xlab = "Separation (nm)", ylab = "Frequency")
  dev.off()

  jpeg(file.path(outDir, "HistShape.jpeg"),
       units = "in", width = 3.5, height = 3.5, quality = 100, res = 300)
  mar.default <- c(5, 4, 4, 2) + 0.1
  par(oma=c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  hist(shape, main = "", xlab = "Displacement (nm)", ylab = "Frequency")
  dev.off()

  }
  return(rupW)
}
