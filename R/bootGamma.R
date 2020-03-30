#' @title Estimate Free Energy Error Using Bootstrap for Gamma
#'
#' @description This function predicts mean and sd of free energy.
#' @param data external work in kT.
#' @keywords gamma
#' @import boot
#' @export
#'

bootGamma <- function(data, replicates = 1000) {

  boot.fn <- function(data, index) {
    x <- data[index]
    x.gamma.mean <- mean(x)
    x.gamma.var <- var(x)
    rate.est <-  x.gamma.mean / x.gamma.var
    shape.est <- x.gamma.mean^2 / x.gamma.var
    FE <- shape.est * log((rate.est + 1) / rate.est)
    return(FE)
  }

  output <- boot(data, boot.fn, replicates)
  return(list(mean = mean(output$t), sd = sd(output$t)))
}
