#' @title Make subdirectories in output directory according to number of peaks
#' @param desDir output directory
#' @export

MakeDirNoPeaks <- function(desDir, n = 4) {
  for (k in 0:n) {
    Cls <- paste0(desDir, "/peakNo", as.character(k))
    if ( dir.exists(Cls) ) {unlink(Cls, recursive <- TRUE)}
    dir.create(Cls)
  }  
}
