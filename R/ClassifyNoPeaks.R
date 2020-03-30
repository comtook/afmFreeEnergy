#' @title The main function classifies data into five classes
#'        according to number of binding points
#'
#' @description This function allows you to classify F-Z curves into five classes.
#' @param datDir data directory
#' @param desDir result directory
#' @param Alpha percent of data points in bulk region
#' @param LOD user given level of noise in standard deviation unit
#' @param nlsEr RMSD tolerence for nls fit in standard deviation unit
#' @param mulBase user given level of baseline perturbation in standard deviation unit
#'
#' @import tools
#' @export
#'

ClassifyNoPeaks <- function(datDir, desDir, Alpha=0.5, LOD=3,
                            nlsEr=4,  mulBase = 2) {
  # make image output directories
  MakeDirNoPeaks(desDir, 4)
  # output container
  myOutput <- NULL
  # main loop
  datFiles <- list.files(path = file.path(datDir))
  numData <- length(datFiles)

  for (i in 1:numData) {
    # read data
    fileName <- datFiles[i]
    rootFileName <- tools::file_path_sans_ext(fileName)
    myData <- NULL
    dat <- read.table(file = file.path(datDir, fileName), header = TRUE)
    # keep experimental record
    myData$Z <- dat$Z
    myData$Force <- dat$Force
    myData$rawForce <- myData$Force
    # print(fileName)
    # find binding events and do alignment
    myData <- FindKeyPoints(myData, Alpha, LOD, mulBase)
    # plot curves according to number of peaks
    k <- length(myData$jBP)
    if ( k < 4 ) {
      Cls <- paste0(desDir, "/peakNo", as.character(k))
      PlotExRt(myData, rootFileName, Cls)
    } else {
      Cls <- paste0(desDir, "/peakNo4")
      PlotExRt(myData, rootFileName, Cls)
    }
    # identify single binding events
    j <- length(myData$jQP)
    if ((k > 1) && (j > 0)) {
      # continue to do worm-like chain fitting
      myData <- FitWLC(myData, nlsEr)
    }
    # record results
    myOutput[[rootFileName]] <- myData
  }
  # output results
  save(myOutput, file = file.path(desDir, 'results.RData'))
}
