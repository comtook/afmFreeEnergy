#' @title Estimate Free Energy Error Using Bootstrap for Norm
#'
#' @description This function predicts mean and sd of free energy.
#' @param data external work in kT.
#' @keywords norm
#' @import boot
#' @export
#'

bootNorm <- function(data, replicates = 1000) {

  boot.fn <- function(data, index) {
    x <- data[index]
    x.mean <- mean(x)
    x.sd <- sd(x)
    FE <- (x.mean - 0.5 * x.sd * x.sd)
    return(FE)
  }

  output <- boot(data, boot.fn, replicates)
  return(list(mean = mean(output$t), sd = sd(output$t)))
}
