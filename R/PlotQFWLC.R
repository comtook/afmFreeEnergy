#' @title Plot F-Q curves
#'
#' @description This function plots force vs distance
#'
#' @param outDir result directory
#'
#' @export
#'
#'

PlotQFWLC <- function(myData, outDir, rootFileName) {
  figName <- paste0(rootFileName, '.jpeg')
  jpeg(file.path(outDir, figName),
       units = "in", width = 3.5, height = 3.5, quality = 100, res = 300)
  mar.default <- c(5, 4, 4, 2) + 0.1
  par(oma=c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  Z <- myData$Z
  Force <- myData$Force
  iRP <- myData$iRP
  iQP <- myData$iQP
  plot(Z, Force, type = "l", 
           xlab="Piezo displacement (nm)",
           ylab="Force (pN)", main="")
  abline(v=Z[c(iRP, iQP)], col='blue')
  lines(myData$WLC, col='red')
  dev.off()
}



