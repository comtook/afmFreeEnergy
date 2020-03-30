#' @title Map F-D curves to Derivative Space
#'
#' @description
#' The derivative are approximated by the slope of a linear fit to the data in a window of
#' a given radius width.
#'
#' @usage CalcSlope(X, width)
#' @param X a matrix on the form [1 z Force],
#' @param width radius of the window for the local regression (in vector position units)
#' @return A vector containing the derivatives.
#'
#' @examples
#' n <- 100
#' x <- seq(0, 2 * pi, length.out = n)
#' y = sin(x) + 0.1 * rnorm(n)
#' X <- matrix(c(rep(1, n), x, y), nrow = n, ncol = 3)
#' width <- 5
#' b <- CalcSlope(X, width)
#' plot(x[(width + 1) : (n - width)], b, xlab = "x", ylab = "y", type = "l")
#' lines(x, y, col = "red")
#' legend("bottomleft", c("Slopes","Signal"), col = c(1, 2), lty = 1)
#' @importFrom stats lm.fit
#' @export
#'

CalcSlope <- function(X, width) {
  lenX <- dim(X)[1]
  SEQ1 <- seq(width + 1, lenX - width)
  SEQ2 <- lapply(SEQ1, function(x) (x - width):(x + width))
  OUT <- lapply(SEQ2, function(a) coef(lm.fit(X[a, 1:2], X[a, 3]))[2])
  OUT <- simplify2array(OUT, higher = TRUE)
  return(OUT)
}


