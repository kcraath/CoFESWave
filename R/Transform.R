#' CoFESWave.Transform
#'
#' @title Computation of the wavelet power spectrum of a single time series
#'
#' @description The time series is selected from an input data frame by specifying either
#' its name or its column number. Optionally, the time series is detrended,
#' using \code{loess} with parameter \code{loess.span}. Internally, the series
#' will be further standardized before it undergoes wavelet transformation.
#'
#' The wavelet power spectrum is computed by applying the Morlet
#' wavelet. P-values to test the null hypothesis that a period
#' (within \code{lowerPeriod} and \code{upperPeriod})
#' is irrelevant at a certain time are calculated if desired;
#' this is accomplished with the help of a simulation algorithm.
#' There is a selection of models from which to
#' choose the alternative hypothesis. The selected model will be fitted to
#' the data and simulated according to estimated parameters in order
#' to provide surrogate time series.
#'
#' Wavelet transformation, as well as p-value computations, are
#' carried out by calling subroutine \code{wt}.
#'
#' The name and parts of the layout of subroutine \code{wt} were inspired by
#' a similar function developed by Huidong Tian and Bernard Cazelles
#' (archived R package \code{WaveletCo}). The basic concept of the simulation algorithm
#' and of ridge determination build on ideas developed by these authors. The major
#' part of the code for the computation of the cone of influence and the code for
#' Fourier-randomized surrogate time series has been adopted from Huidong Tian.
#'
#' Wavelet computation, the simulation algorithm and ridge determination build
#' heavily on the use of matrices in order to minimize computation time in R.
#'
#' This function provides a broad variety of final as well as intermediate results
#' which can be further analyzed in detail.
#'
#'
#' @param my.data data frame of time series (including header, and dates as row
#' names or as separate columnnamed \code{"date"} if available)
#' @param my.series name or column index indicating the series to be analyzed,
#' e.g. \code{1}, \code{2}, \code{"dji"}, \code{"ftse"}.
#'
#' Default: \code{1}.
#' @param loess.span parameter \code{alpha} in \code{loess} controlling the degree
#' of time series smoothing, if the time series is to be detrended; no detrending
#' if \code{loess.span = 0}.
#'
#' Default: \code{0.75}.
#' @param dt time resolution, i.e. sampling resolution in the time domain,
#' \code{1/dt} = number of observations per time unit.
#' For example: a natural choice of \code{dt} in case of hourly data is \code{dt = 1/24},
#' resulting in one time unit equaling one day. This is also the time unit in which periods are measured.
#' If \code{dt = 1}, the time interval between two consecutive observations will equal one time unit.
#'
#' Default: \code{1}.
#' @param dj frequency resolution, i.e. sampling resolution in the frequency domain,
#' \code{1/dj} = number of suboctaves (voices per octave).
#'
#' Default: \code{1/20}.
#' @param lowerPeriod lower Fourier period (measured in time units determined by \code{dt},
#' see the explanations concerning \code{dt}) for wavelet decomposition.\cr If \code{dt = 1},
#' the minimum admissible value is 2.
#'
#' Default: \code{2*dt}.
#' @param upperPeriod upper Fourier period (measured in time units determined by \code{dt},
#' see the explanations concerning \code{dt}) for wavelet decomposition.
#'
#' Default: \code{(floor of one third of time series length)*dt}.
#' @param make.pval Compute p-values? Logical.
#'
#' Default: \code{TRUE}.
#' @param method the method of generating surrogate time series; select from:
#' \tabular{rlll}{
#'   \tab \code{"white.noise"}  \tab : \tab white noise \cr
#'   \tab \code{"shuffle"}      \tab : \tab shuffling the given time series \cr
#'   \tab \code{"Fourier.rand"} \tab : \tab time series with a similar spectrum \cr
#'   \tab \code{"AR"}           \tab : \tab AR(p) \cr
#'   \tab \code{"ARIMA"}        \tab : \tab ARIMA(p,0,q)
#' }
#'
#' Default: \code{"white.noise"}.
#' @param params
#' a list of assignments between methods (AR, and ARIMA) and lists of parameter values
#' applying to surrogates. Default: \code{NULL}.
#'
#' Default includes two lists named \code{AR} and \code{ARIMA}:
#'   \itemize{
#'     \item \code{AR = list(...)}, a list containing one single element:
#'       \tabular{rlll}{
#'         \tab \code{p} \tab : \tab AR order.\cr
#'         \tab          \tab   \tab Default: \code{1}.
#'       }
#'     \item \code{ARIMA = list(...)}, a list of six elements:
#'       \tabular{rlll}{
#'         \tab \code{p}            \tab : \tab  AR order. \cr
#'         \tab                     \tab   \tab  Default: \code{1}.\cr
#'         \tab \code{q}            \tab : \tab  MA order. \cr
#'         \tab                     \tab   \tab  Default: \code{1}.\cr
#'         \tab \code{include.mean} \tab : \tab  Include a mean/intercept term? \cr
#'         \tab                     \tab   \tab  Default: \code{TRUE}.\cr
#'         \tab \code{sd.fac}       \tab : \tab  magnification factor to boost the \cr
#'         \tab                     \tab   \tab  residual standard deviation. \cr
#'         \tab                     \tab   \tab  Default: \code{1}.\cr
#'         \tab \code{trim}         \tab : \tab  Simulate trimmed data? \cr
#'         \tab                     \tab   \tab  Default: \code{FALSE}.\cr
#'         \tab \code{trim.prop}    \tab : \tab  high/low trimming proportion. \cr
#'         \tab                     \tab   \tab  Default: \code{0.01}.\cr
#'       }}
#' @param n.sim number of simulations.
#'
#' Default: \code{100}.
#' @param date.format optional, and for later reference: the format of calendar date
#' (if available in the input data frame) given as a character string, e.g. \code{"\%Y-\%m-\%d"},
#' or equivalently \code{"\%F"}; see \code{strptime} for a list of implemented date conversion specifications.
#' Explicit information given here will be overwritten by any later specification given in
#' e.g. \code{wt.image}. If unspecified, date formatting will be attempted according to \code{as.Date}.
#'
#' Default: \code{NULL}.
#' @param date.tz optional, and for later reference: a character string specifying the time zone of
#' calendar date (if available in the input data frame); see \code{strptime}.
#' Explicit information given here will be overwritten by any specification given in e.g. \code{wt.image}.
#'
#' If unspecified, \code{""} (the local time zone) will be used.
#' Default: \code{NULL}.
#' @param verbose Print verbose output on the screen? Logical.
#'
#' Default: \code{TRUE}.
#'
#' @details Wavelet transformation, as well as p-value computations, are
#' carried out by calling the internal function \code{wt}.
#' @return A list of class \code{"analyze.wavelet"} with elements of different dimensions.
#' %%%%%%%%%%%%%%%%%
#' The elements of matrix type (namely, \code{Wave}, \code{Phase}, \code{Ampl},
#' \code{Power}, \code{Power.pval}, \code{Ridge}) have the following structure:
#' \cr columns correspond to observations (observation epochs; "epoch" meaning point in time),
#' rows correspond to scales (Fourier periods) whose values are given in \code{Scale} (\code{Period}).
#' %%%%%%%%%%%%%%%%%
#' Here is a detailed list of all elements:
#'   \item{series}{a data frame with the following columns:
#'
#'       \tabular{rlll}{
#'         \tab date        \tab : \tab the calendar date \cr
#'         \tab             \tab   \tab (if available as column in \code{my.data}) \cr
#'         \tab <x>         \tab : \tab the series which has been analyzed \cr
#'         \tab             \tab   \tab (detrended, if \code{loess.span != 0}; \cr
#'                                       \tab             \tab   \tab  original name retained) \cr
#'         \tab <x>.trend   \tab : \tab the trend series (if \code{loess.span != 0})
#'       }
#'
#'     Row names are taken over from \code{my.data}, and so are dates if given as row names.
#'   }
#' %%%%%%%%%%%%%%%%%
#' \item{loess.span}{parameter \code{alpha} in \code{loess} controlling the degree of time series smoothing
#'   if the time series was detrended; no detrending if \code{loess.span = 0}}
#' %%%%%%%%%%%%%%%%%
#' \item{dt}{time resolution, i.e. sampling resolution in the time domain, \code{1/dt} = number of observations per time unit}
#' \item{dj}{frequency resolution, i.e. sampling resolution in the frequency domain, \code{1/dj} = number of suboctaves (voices per octave)}
#' %%%%%%%%%%%%%%%%%
#' \item{Wave}{complex wavelet transform of the series}
#' \item{Phase}{phases}
#' \item{Ampl}{amplitudes}
#' %%%%%%%%%%%%%%%%%
#' \item{Power}{wavelet power in the time/frequency domain}
#' \item{Power.avg}{average wavelet power in the frequency domain (averages over time)}
#' \item{Power.pval}{p-values of wavelet power}
#' \item{Power.avg.pval}{p-values of average wavelet power}
#' %%%%%%%%%%%%%%%%%
#' \item{Ridge}{wavelet power ridge, in the form of a matrix of \code{0}s and \code{1}s}
#' %%%%%%%%%%%%%%%%%
#' \item{Period}{the Fourier periods
#'   (measured in time units determined by \code{dt}, see the explanations concerning \code{dt})}
#' \item{Scale}{the scales (the Fourier periods divided by the Fourier factor)}
#' %%%%%%%%%%%%%%%%%
#' \item{nc}{number of columns = number of observations = number of observation epochs; "epoch" meaning point in time}
#' \item{nr}{number of rows = number of scales (Fourier periods)}
#' %%%%%%%%%%%%%%%%%
#' \item{coi.1, coi.2}{borders of the region where the wavelet transforms are not influenced by edge effects (cone of influence).
#'   The coordinates of the borders are expressed in terms of internal axes \code{axis.1} and \code{axis.2}.}
#' %%%%%%%%%%%%%%%%%
#' \item{axis.1}{tick levels corresponding to the time steps used for (cross-)wavelet transformation: \code{1, 1+dt, 1+2dt, ...}.
#'   The default time axis in plot functions provided by \code{WaveletComp} is determined by observation epochs, however; "epoch" meaning point in time. }
#' \item{axis.2}{tick levels corresponding to the log of Fourier periods: \code{log2(Period)}. This determines the period axis in plot functions provided by \code{WaveletComp}.}
#' %%%%%%%%%%%%%%%%%
#' \item{date.format}{the format of calendar date (if available)}
#' \item{date.tz}{the time zone of calendar date (if available)}
#'
#' @author CoFES. Credits are also due toAngi Roesch, Harald Schmidbauer, Huidong Tian, and Bernard Cazelles.
#'
#' @references
#' Aguiar-Conraria L., and Soares M.J., 2011. The Continuous Wavelet Transform: A Primer.
#' NIPE Working Paper Series 16/2011.
#'
#' Carmona R., Hwang W.-L., and Torresani B., 1998. Practical Time Frequency Analysis.
#' Gabor and Wavelet Transforms with an Implementation in S.Academic Press, San Diego.
#'
#' Cazelles B., Chavez M., Berteaux, D., Menard F., Vik J.O., Jenouvrier S., and
#' Stenseth N.C., 2008. Wavelet analysis of ecological time series. Oecologia 156, 287--304.
#'
#' Liu Y., Liang X.S., and Weisberg R.H., 2007. Rectification of the Bias in the
#' Wavelet Power Spectrum.  Journal of Atmospheric and Oceanic Technology 24, 2093--2102.
#'
#' Tian, H., and Cazelles, B., 2012. \code{WaveletCo}. Available
#' at \url{https://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013;
#' accessed July 26, 2013.
#'
#' Torrence C., and Compo G.P., 1998. A practical guide to wavelet analysis.
#' Bulletin of the American Meteorological Society 79 (1), 61--78.
#'
#' @seealso \code{\link{CoFESWave.reconstruct}}
#'
#' @examples \dontrun{
#' ## The following example is adopted from Liu et al., 2007:
#'
#' series.length <- 6*128*24
#' x1 <- periodic.series(start.period = 1*24, length = series.length)
#' x2 <- periodic.series(start.period = 8*24, length = series.length)
#' x3 <- periodic.series(start.period = 32*24, length = series.length)
#' x4 <- periodic.series(start.period = 128*24, length = series.length)
#'
#' x <- x1 + x2 + x3 + x4
#'
#' plot(x, type = "l", xlab = "index", ylab = "", xaxs = "i",
#'      main = "hourly series with periods of 1, 8, 32, 128 days")
#'
#' ## The following dates refer to the local time zone
#' ## (possibly allowing for daylight saving time):
#' my.date <- seq(as.POSIXct("2014-10-14 00:00:00", format = "\%F \%T"),
#'                by = "hour",
#'                length.out = series.length)
#' my.data <- data.frame(date = my.date, x = x)
#'
#' ## Computation of wavelet power:
#' ## a natural choice of 'dt' in the case of hourly data is 'dt = 1/24',
#' ## resulting in one time unit equaling one day.
#' ## This is also the time unit in which periods are measured.
#' ## There is an option to store the date format and time zone as additional
#' ## parameters within object 'my.wt' for later reference.
#'
#' my.wt <- analyze.wavelet(my.data, "x",
#'                          loess.span = 0,
#'                          dt = 1/24, dj = 1/20,
#'                          lowerPeriod = 1/4,
#'                          make.pval = TRUE, n.sim = 10,
#'                          date.format = "\%F \%T", date.tz = "")
#'
#' ## Plot of wavelet power spectrum (with equidistant color breakpoints):
#' wt.image(my.wt, color.key = "interval", main = "wavelet power spectrum",
#'          legend.params = list(lab = "wavelet power levels"),
#'          periodlab = "period (days)")
#'
#' ## Plot of average wavelet power:
#' wt.avg(my.wt, siglvl = 0.05, sigcol = "red",
#'        periodlab = "period (days)")
#'
#' ## Please see our guide booklet for further examples:
#' ## URL http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf.
#'
#' }
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
