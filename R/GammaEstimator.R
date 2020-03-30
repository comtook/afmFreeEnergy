#' @title Estimate Free Energy Assuming Underlying Gamma Distribution
#'
#' @description This function tests underlying gamma distribution of work and estimates free energy.
#' @param wrkDir Directory for image output.
#' @keywords gamma
#' @import fitdistrplus
#' @export
#'

GammaEstimator <- function(wrkDir, data) {

  # ml estimate
  fit.gamma <- fitdist(data, distr = "gamma", method = "mle")
  alpha <- fit.gamma$estimate[1]   # shape
  rate <- fit.gamma$estimate[2]   # rate
  alpha.sd <- fit.gamma$sd[1]
  rate.sd <- fit.gamma$sd[2]
  FE <- alpha * log((rate + 1) / rate)   # unit kT

  # estimate error
  FE.boot <- bootGamma(data, 10000)

  # Kolmogorov-Smirnov test
  kS <- ks.test(data, "pgamma", shape = alpha, rate = rate)

  output <- list(FE.boost = FE.boot$mean,
                 FE.boost.sd = FE.boot$sd,
                 FE.mle = FE,
                 shape = alpha,
                 rate = rate,
                 shape.sd = alpha.sd,
                 rate.sd = rate.sd,
                 KS_D=kS$statistic,
                 KS_P=kS$p.value)

  # standard plot from fitdistrplus
  minX <- min(data)
  maxX <- max(data)
  # cumulative density function
  jpeg(file.path(wrkDir, 'cdfGamma.jpeg'), units = "in",
     width = 3.25, height = 3.25, quality = 100, res = 300)
  # default margin
  mar.default <- c(5,4,4,2) + 0.1
  # set custom margin
  par(oma = c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  cdfcomp(fit.gamma, lwd = 2, cex = 0.8,
        horizontals = TRUE, verticals = TRUE,
        xlim = c(minX, maxX), ylim = c(0, 1), main = "",
        xlab = expression(paste("Work (k"[B], "T)")))
  dev.off()

  # Q-Q plot
  jpeg(file.path(wrkDir, 'qqGamma.jpeg'), units = "in",
     width = 3.25, height = 3.25, quality = 100, res = 300)
  # default margin
  mar.default <- c(5,4,4,2) + 0.1
  # set custom margin
  par(oma = c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  qqcomp(fit.gamma)
  dev.off()

  # P-P plot
  jpeg(file.path(wrkDir, 'ppGamma.jpeg'), units = "in",
     width = 3.25, height = 3.25, quality = 100, res = 300)
  # default margin
  mar.default <- c(5,4,4,2) + 0.1
  # set custom margin
  par(oma = c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  ppcomp(fit.gamma)
  dev.off()

  # histograms
  jpeg(file.path(wrkDir, 'histGamma.jpeg'), units = "in",
     width = 3.25, height = 3.25, quality = 100, res = 300)
  # default margin
  mar.default <- c(5,4,4,2) + 0.1
  # set custom margin
  par(oma = c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  denscomp(fit.gamma,
        xlim = c(minX, maxX),  main = "",
        xlab = expression(paste("Work (k"[B], "T)")))
  dev.off()

  # optimized plot
  # determine optimal bin size to fit density curve
  dd.old <- 1.0e+9
  for (cells in 11:29) {
    histOut <- hist(data, breaks = cells, freq = FALSE, plot = FALSE)
    hX <- (histOut$breaks[-1] +
             histOut$breaks[-length(histOut$breaks)]) * 0.5
    hY <- histOut$density
    d = dgamma(hX, shape = alpha, rate = rate)
    dd.new <- sum((d-hY)^2) / length(d)
    if (dd.new < dd.old) {
      dd.old <- dd.new
      binSize <- hX[2] - hX[1]
      maxd <- max(max(d), max(hY))
    }
  }

  xMin <- as.integer(as.integer(min(data) / binSize) * binSize - 1)
  xMax <- as.integer(as.integer(max(data) / binSize) * binSize + 1)
  cells <- as.integer((xMax-xMin) / binSize)
  # histograms
  jpeg(file.path(wrkDir, 'optimalGammaHist.jpeg'), units = "in",
       width = 3.5, height = 3.5, quality = 100, res = 300)
  # default margin
  mar.default <- c(5,4,4,2) + 0.1
  # set custom margin
  par(oma = c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  histOut <- hist(data, breaks = cells, prob = TRUE,
                  main = "",  ylab = "Density",
                  xlab = expression(paste("Work (", k[B], "T)")),
                  xlim = c(xMin, xMax), ylim = c(0, maxd))
  curve(dgamma(x, shape = alpha, rate = rate),
        col = "red", lwd = 2, add = TRUE, xaxt = "n", yaxt = "n")
  legend('topright',
         legend = c("gamma"),
         col = c("red"),
         lty = c(1), cex = 0.75,  lwd = 2,
         inset = c(0.01),
         box.lty = 0)
  #title(xlab = XLAB, line = 2)
  dev.off()
  # CDF
  x <- sort(data)
  x.length <- length(x)
  y <- pgamma(x, shape = alpha, rate = rate)
  jpeg(file.path(wrkDir, 'optimalCDF.jpeg'), units = "in",
       width = 3.5, height = 3.5, quality = 100, res = 300)
  # default margin
  mar.default <- c(5,4,4,2) + 0.1
  # set custom margin
  par(oma = c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  plot(x, y, type = 'l', lty=1, lwd=2, ylim = c(0, 1),
       ylab = "CDF", col = 'red',
       xlab = expression(paste("Work (", k[B], "T)")))
  points(x, ((1:x.length) - 0.5) / x.length)
  legend('bottomright',
         legend = c("gamma"),
         col = c("red"),
         lty = c(1), cex = 0.75,  lwd = 2,
         inset = c(0.01),
         box.lty = 0)
  # title(xlab = XLAB, line = 2)
  dev.off()


  return(output)
}

