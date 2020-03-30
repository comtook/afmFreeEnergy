#' @title Estimate Free Energy Assuming Underlying Norm Distribution
#'
#' @description This function tests underlying norm distribution of work and estimates free energy.
#' @param wrkDir Directory for image output.
#' @keywords norm
#' @import fitdistrplus
#' @export
#'

GaussianEstimator <- function(wrkDir, data) {

  # calculate free energy
  fit.norm <- fitdist(data, distr = "norm", method = "mle")
  mu <- fit.norm$estimate[1]   # mean
  sigma <- fit.norm$estimate[2]   # sigma
  mu.sd <- fit.norm$sd[1]
  sigma.sd <- fit.norm$sd[2]
  FE <- (mu - 0.5 * sigma * sigma)   # unit kT

  # estimate error
  FE.boot <- bootNorm(data, 10000)

  # Kolmogorov-Smirnov test
  kS <- ks.test(data, "pnorm", mean = mu, sd = sigma)

  output <- list(FE.boost = FE.boot$mean,
                 FE.boost.sd = FE.boot$sd,
                 FE.mle = FE,
                 mean = mu,
                 sd = sigma,
                 mean.sd = mu.sd,
                 sd.sd = sigma.sd,
                 KS_D = kS$statistic,
                 KS_P = kS$p.value)

  # standard plot by fitdistrplus
  minX <- min(data)
  maxX <- max(data)
  # cumulative density function
  jpeg(file.path(wrkDir, 'cdfNorm.jpeg'), units = "in",
     width = 3.25, height = 3.25, quality = 100, res = 300)
  # default margin
  mar.default <- c(5,4,4,2) + 0.1
  # set custom margin
  par(oma = c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")

  cdfcomp(fit.norm, lwd = 2, cex = 0.8,
        xlim = c(minX, maxX), ylim = c(0, 1), main = "",
        xlab = expression(paste("Work (k"[B], "T)")))

  dev.off()

  # Q-Q plot
  jpeg(file.path(wrkDir, 'qqNorm.jpeg'), units = "in",
     width = 3.25, height = 3.25,quality = 100, res = 300)
  # default margin
  mar.default <- c(5,4,4,2) + 0.1
  # set custom margin
  par(oma = c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  qqcomp(fit.norm)
  dev.off()

  # P-P plot
  jpeg(file.path(wrkDir, 'ppNorm.jpeg'), units = "in",
     width = 3.25, height = 3.25,quality = 100, res = 300)
  # default margin
  mar.default <- c(5,4,4,2) + 0.1
  # set custom margin
  par(oma = c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")

  ppcomp(fit.norm)
  dev.off()

  # histograms
  jpeg(file.path(wrkDir, 'histNorm.jpeg'), units = "in",
     width = 3.25, height = 3.25,quality = 100, res = 300)
  # default margin
  mar.default <- c(5,4,4,2) + 0.1
  # set custom margin
  par(oma = c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")

  denscomp(fit.norm,
        xlim = c(minX, maxX),  main = "",
        xlab = expression(paste("Work (k"[B], "T)")))

  dev.off()

  # determine optimal bin size to fit density curve
  dd.old <- 1.0e+9
  for (cells in 11:29) {
    histOut <- hist(data, breaks = cells, freq = FALSE, plot = FALSE)
    hX <- (histOut$breaks[-1] +
             histOut$breaks[-length(histOut$breaks)]) * 0.5
    hY <- histOut$density
    d = dnorm(hX, mean = mu, sd = sigma)
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
  jpeg(file.path(wrkDir, 'optimalNormHist.jpeg'), units = "in",
       width = 3.5, height = 3.5, quality = 100, res = 300)
  # default margin
  mar.default <- c(5,4,4,2) + 0.1
  # set custom margin
  par(oma = c(0, 0, 0, 0), mar = mar.default + c(-1, 0, -2, 0), las = "1")
  histOut <- hist(data, breaks = cells, prob = TRUE,
                  main = "",  ylab = "Density",
                  xlab = expression(paste("Work (", k[B], "T)")),
                  xlim = c(xMin, xMax), ylim = c(0, maxd))
  curve(dnorm(x, mean = mu, sd = sigma),
        col = "red", lwd = 2, add = TRUE, xaxt = "n", yaxt = "n")
  legend('topright',
         legend = c("norm"),
         col = c("red"),
         lty = c(1), cex = 0.75,  lwd = 2,
         inset = c(0.01),
         box.lty = 0)
  #title(xlab = XLAB, line = 2)
  dev.off()
  # CDF
  x <- sort(data)
  x.length <- length(x)
  y <- pnorm(x, mean = mu, sd = sigma)
  jpeg(file.path(wrkDir, 'optimalNormCDF.jpeg'), units = "in",
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
         legend = c("norm"),
         col = c("red"),
         lty = c(1), cex = 0.75,  lwd = 2,
         inset = c(0.01),
         box.lty = 0)
  # title(xlab = XLAB, line = 2)
  dev.off()

  return(output)

}

