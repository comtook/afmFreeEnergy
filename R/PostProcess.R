#' @title Get rid of outliers
#'
#' @description This function deletes all outliers
#'
#' @param outDir Directory is where work.RData locates.
#'
#' @export
#'

PostProcess <- function(outDir, QL=0, QR=200, FL=-1000, FR=0, WL=0, WR=1000) {
  # load data
  load(file=file.path(outDir, 'work.RData'))
  x <- x[(x$nm < QR),]
  x <- x[(x$nm > QL),]
  x <- x[(x$pN < FR),]
  x <- x[(x$pN > FL),]
  x <- x[(x$kT > WL),]
  x <- x[(x$kT < WR),]
  
  save(x, file = file.path(outDir, 'selectedData.RData'))
  
  # plot histogram
  jpeg(file.path(outDir, "tunedHistWork.jpeg"), 
       units = "in", width = 3.5, height = 3.5, quality = 100, res = 300)
  mar.default <- c(5, 4, 4, 2) + 0.1
  par(oma=c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  hist(x$kT, main = "", xlab = expression(paste("Work (k"[B], "T)")), ylab = "Frequency")
  dev.off()
  
  jpeg(file.path(outDir, "tunedHistForce.jpeg"), 
       units = "in", width = 3.5, height = 3.5, quality = 100, res = 300)
  mar.default <- c(5, 4, 4, 2) + 0.1
  par(oma=c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  hist(x$pN, main = "", xlab = "Force (pN)", ylab = "Frequency")
  dev.off()
  
  jpeg(file.path(outDir, "tunedHistQ.jpeg"), 
       units = "in", width = 3.5, height = 3.5, quality = 100, res = 300)
  mar.default <- c(5, 4, 4, 2) + 0.1
  par(oma=c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  hist(x$nm, main = "", xlab = "Separation (nm)", ylab = "Frequency")
  dev.off()
  
  return(x$kT)
}
