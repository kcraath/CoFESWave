#' Title
#'
#' @param my.data
#' @param my.series
#' @param loess.span
#' @param dt
#' @param dj
#' @param lowerPeriod
#' @param upperPeriod
#' @param make.pval
#' @param method
#' @param params
#' @param n.sim
#' @param date.format
#' @param date.tz
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
CoFESWave.Transform <- function (my.data,
                                 my.series = 1, loess.span = 0.75,
                                 dt = 1, dj = 1/20,
                                 lowerPeriod = 2 * dt, upperPeriod = floor(nrow(my.data)/3) * dt,
                                 make.pval = TRUE,
                                 method = "white.noise", params = NULL, n.sim = 100,
                                 date.format = NULL, date.tz = NULL,
                                 verbose = TRUE)
{
  #### 0. Input Check ####

  ## 0.1. Verbose
  if (verbose == T) {
    out <- function(...) {
      cat(...)
    }
  }
  else {
    out <- function(...) {
    }
  }

  ## 0.1. Series

  ## 0.1.1. Preparation
  loess.data.frame = function(x, loess.span) {
    x.smoothed = x
    for (i in 1:ncol(x)) {
      day.index = 1:nrow(x)
      my.loess.x = loess(x[, i] ~ day.index, span = loess.span)
      x.loess = as.numeric(predict(my.loess.x, data.frame(x = 1:nrow(x))))
      x.smoothed[, i] = x.loess
    }
    return(x.smoothed)
  }
  ind = which(names(my.data) == my.series)
  x = data.frame(my.data[, ind])
  colnames(x) = my.series
  rownames(x) = rownames(my.data)

  ## 0.1.2. Test the series
  if (is.numeric(my.series)) {
    my.series = names(my.data)[my.series]
  }
  if (length(my.series) != 1) {
    stop("Please select (only) one series for analysis!\n")
  }
  if (is.element("date", my.series)) {
    stop("Please review your selection of series!\n")
  }
  if (!is.numeric(x[[my.series]])) {
    stop("Some values in your time series do not seem to be interpretable as numbers.\n")
  }
  if (sum(is.na(x[[my.series]])) > 0) {
    stop("Some values in your time series seem to be missing.\n")
  }
  if (sd(x[[my.series]]) == 0) {
    stop("Your time series seems to be constant, there is no need to search for periodicity.\n")
  }
  if (lowerPeriod > upperPeriod) {
    stop("Please choose lowerPeriod smaller than or (at most) equal to upperPeriod.\n")
  }
  if (loess.span != 0) {
    out("Smoothing the time series...\n")
    x.trend = loess.data.frame(x, loess.span)
    x = x - x.trend
    x = cbind(x, x.trend)
    colnames(x) = c(my.series, paste(my.series, ".trend",
                                     sep = ""))
  }
  if (is.element("date", names(my.data))) {
    x = cbind(date = my.data$date, x)
  }
  out("Starting wavelet transformation...\n")
  if (make.pval == T) {
    out("... and simulations... \n")
  }

  ### End of 0. Input Check ###



  #### 1. Functions ####

  ## 1.1. Wavelet Tranformation
  wt <- function (x, start = 1,
                  dt = 1, dj = 1/20,
                  lowerPeriod = 2 * dt, upperPeriod = floor(length(x) * dt/3),
                  make.pval = TRUE,
                  method = "white.noise",
                  params = NULL, n.sim = 100, save.sim = FALSE)
  {
    ## 1. Sub-functions: Morlet Wavelet Transform
    morlet.wavelet.transform = function(x) {
      x = (x - mean(x))/sd(x)
      xpad = c(x, rep(0, pad.length))
      fft.xpad = fft(xpad)
      wave = matrix(0, nrow = scales.length, ncol = N)
      wave = wave + (0+1i) * wave
      for (ind.scale in (1:scales.length)) {
        my.scale = scales[ind.scale]
        norm.factor = pi^(1/4) * sqrt(2 * my.scale/dt)
        expnt = -((my.scale * omega.k - omega0)^2/2) * (omega.k >
                                                          0)
        daughter = norm.factor * exp(expnt)
        daughter = daughter * (omega.k > 0)
        wave[ind.scale, ] = fft(fft.xpad * daughter, inverse = TRUE)/N
      }
      wave = wave[, 1:series.length]
      return(wave)
    }

    ## 2. Sub-functions: SurrogateData
    SurrogateData <- function (x, method = "white.noise",
                               params = list(AR = list(p = 1),
                                             ARIMA = list(p = 1, q = 1, include.mean = TRUE,
                                                          sd.fac = 1,trim = FALSE, trim.prop = 0.01)))
    {
      if (method == "white.noise")
        x.sur <- rnorm(length(x))
      if (method == "shuffle")
        x.sur <- sample(x, length(x))
      if (method == "Fourier.rand")
        x.sur <- FourierRand(x)
      if (method == "AR") {
        x.sur <- AR(x, params = params)
      }
      if (method == "ARIMA") {
        x.sur <- ARIMA(x, params = params)
      }
      return(invisible(x.sur))
    }

    ## 3. Sub-functions: COI
    COI <- function (start = start, dt = dt, nc = nc, nr = nr, Period = Period)
    {
      axis.1 = seq(from = start, by = dt, length.out = nc)
      axis.2 = log2(Period)
      omega0 = 6
      fourier.factor = (2 * pi)/omega0
      coi = fourier.factor * sqrt(2) * dt * c(1e-05, 1:((nc + 1)/2 -
                                                          1), rev((1:(nc/2 - 1))), 1e-05)
      coi.x = c(axis.1[c(1, 1)] - dt * 0.5, axis.1, axis.1[c(nc,
                                                             nc)] + dt * 0.5)
      logyint = axis.2[2] - axis.2[1]
      yl = c(log2(Period[nr]) + 0.5 * logyint, log2(Period[1]) -
               0.5 * logyint)
      yr = rev(yl)
      coi.y = c(yl, log2(coi), yr)
      output <- list(x = coi.x, y = coi.y, axis.1 = axis.1, axis.2 = axis.2)
      return(invisible(output))
    }

    ## 4. Preparation: Wavelet Tranformation
    series.length = length(x)
    pot2 = trunc(log2(series.length) + 0.5)
    pad.length = 2^(pot2 + 1) - series.length
    omega0 = 6
    fourier.factor = (2 * pi)/omega0
    min.scale = lowerPeriod/fourier.factor
    max.scale = upperPeriod/fourier.factor
    J = as.integer(log2(max.scale/min.scale)/dj)
    scales = min.scale * 2^((0:J) * dj)
    scales.length = length(scales)
    periods = fourier.factor * scales
    N = series.length + pad.length
    omega.k = 1:floor(N/2)
    omega.k = omega.k * (2 * pi)/(N * dt)
    omega.k = c(0, omega.k, -omega.k[floor((N - 1)/2):1])

    ## 5. Wavelet Tranformation
    Wave = morlet.wavelet.transform(x)
    Power = Mod(Wave)^2/matrix(rep(scales, series.length), nrow = scales.length)
    Phase = Arg(Wave)
    Ampl = Mod(Wave)/matrix(rep(sqrt(scales), series.length),
                            nrow = scales.length)
    Power.avg <- rowMeans(Power)
    Period <- periods
    Scale <- scales
    nr <- scales.length
    nc <- series.length

    ## 6. Simulate p-value
    Power.pval = NULL
    Power.avg.pval = NULL
    series.sim = NULL
    if (make.pval == T) {
      Power.pval = matrix(0, nrow = nr, ncol = nc)
      Power.avg.pval = rep(0, nr)
      if (save.sim == T) {
        series.sim = matrix(NA, nrow = nc, ncol = n.sim)
      }
      pb = txtProgressBar(min = 0, max = n.sim, style = 3)
      for (ind.sim in 1:n.sim) {
        x.sim = SurrogateData(x, method = method)
        if (save.sim == T) {
          series.sim[, ind.sim] = x.sim
        }
        WT.sim = morlet.wavelet.transform(x)
        Power.sim = Mod(WT.sim)^2/matrix(rep(scales, series.length), nrow = scales.length)
        Power.avg.sim = rowMeans(Power.sim)
        rm(WT.sim)
        Power.pval[Power.sim >= Power] = Power.pval[Power.sim >=
                                                      Power] + 1
        Power.avg.pval[Power.avg.sim >= Power.avg] = Power.avg.pval[Power.avg.sim >=
                                                                      Power.avg] + 1
        setTxtProgressBar(pb, ind.sim)
      }
      close(pb)
      Power.pval = Power.pval/n.sim
      Power.avg.pval = Power.avg.pval/n.sim
    }

    ## 7. COI
    coi = COI(start = start, dt = dt, nc = nc, nr = nr, Period = Period)

    ## 8. Output
    output = list(Wave = Wave, Phase = Phase, Ampl = Ampl, Power = Power,
                  Power.avg = Power.avg, Power.pval = Power.pval, Power.avg.pval = Power.avg.pval,
                  Period = Period, Scale = Scale, coi.1 = coi$x, coi.2 = coi$y,
                  nc = nc, nr = nr, axis.1 = coi$axis.1, axis.2 = coi$axis.2,
                  series.sim = series.sim)
    return(invisible(output))
  }
  ### End of 1. Functions ###



  #### 2. Main function ####
  my.wt = wt(x = x[[my.series]], start = 1, dt = dt, dj = dj,
             lowerPeriod = lowerPeriod, upperPeriod = upperPeriod,
             make.pval = make.pval, method = method, params = params,
             n.sim = n.sim, save.sim = F)
  Ridge = ridge(my.wt$Power)

  ### End of 2. Main function ###



  #### 3. Output ####
  output <- list(series = x, loess.span = loess.span,
                 dt = dt, dj = dj,
                 Wave = my.wt$Wave,
                 Phase = my.wt$Phase,
                 Ampl = my.wt$Ampl,
                 Power = my.wt$Power,
                 Power.avg = my.wt$Power.avg,
                 Power.pval = my.wt$Power.pval,
                 Power.avg.pval = my.wt$Power.avg.pval,
                 Ridge = Ridge,
                 Period = my.wt$Period,
                 Scale = my.wt$Scale,
                 nc = my.wt$nc, nr = my.wt$nr,
                 coi.1 = my.wt$coi.1, coi.2 = my.wt$coi.2,
                 axis.1 = my.wt$axis.1, axis.2 = my.wt$axis.2,
                 date.format = date.format, date.tz = date.tz)

  class(output) = "CoFESWave.Transform"
  # out("Class attributes are accessible through following names:\n")
  # out(names(output), "\n")
  return(invisible(output))
}
