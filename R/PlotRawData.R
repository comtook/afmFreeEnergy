#' @title Plot raw data in images directory
#' @export
#'

  PlotRawData <- function(datDir, outDir, fileType) {
  # main loop
  datFiles <- list.files(path = file.path(datDir), pattern = fileType)
  numData <- length(datFiles)

  for (i in 1:numData) {
    fileName <- datFiles[i]
    myData <- ReadNanoscopeData(fileName, datDir)
    
    jpeg(file.path(outDir, paste0(fileName, ".jpeg")),
         units = "in", width = 3.5, height = 3.5, res = 300)
    mar.default <- c(5, 4, 4, 2) + 0.1
    par(oma=c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
    
    plot(myData$Z, myData$Force,
           xlab="Piezo displacement (nm)",
           ylab="Force (pN)", main="")
    lines(myData$Zapp, myData$Fapp, col='red')
    
    dev.off()
  }
}

