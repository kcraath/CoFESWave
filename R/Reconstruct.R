#' CoFESWave.reconstruct
#'
#' @title Reconstruction of a (detrended) time series from output provided by an
#' object of class \code{"analyze.wavelet"} or \code{"analyze.coherency"}
#'
#' @description This function reconstructs a (detrended) time series analyzed
#' by wavelet transformation using either function
#' \code{analyze.wavelet} or function \code{analyze.coherency}, subject to
#' optional criteria concerning: minimum wavelet power, significance of wavelet power at a
#' given significance level, specification of (Fourier) periods or
#' period bands, exclusive use of the power ridge and/or the cone of influence.
#' An option is provided to prevent the reconstructed series from final rescaling
#' (applying the original (detrended) series' mean and standard deviation).
#'
#' (If the object provided as input is of class \code{"analyze.coherency"},
#' then the number or name of the time series must be specified.)
#'
#' Optional: plot of wavelets used for reconstruction, plot of
#' reconstructed series against original (detrended)
#' series. An option is given to individualize the time axis
#' by specifying tick marks and labels.
#'
#' Output includes the original (detrended) and the reconstructed
#' time series, along with reconstruction waves and parameters.
#'
#' @param WT an object of class \code{"analyze.wavelet"} or \code{"analyze.coherency"}}
#' @param my.series In case \code{class(WT) = "analyze.coherency"}: number (\code{1} or
#' \code{2}) or name of the series to be analyzed.
#'
#' Default: \code{1}.
#' @param lvl minimum level of wavelet power to be applied for the inclusion of reconstruction waves.
#'
#' Default: \code{0}.
#' @param only.ridge Select only the wavelet power ridge? Logical.
#'
#' Default: \code{FALSE}.
#' @param only.sig Use wavelet power significance in reconstruction? Logical.
#'
#' Default: \code{TRUE}.
#' @param siglvl level of wavelet power significance to be applied for the inclusion of reconstruction waves.
#'
#' Default: 0.05.
#' @param only.coi Exclude borders influenced by edge effects in reconstruction, i.e. include the cone of influence only? Logical.
#'
#' Default: \code{FALSE}.
#' @param sel.period a vector of numbers to select Fourier periods (or closest available periods)
#' and corresponding wavelets for the reconstruction.
#'
#' Default: \code{NULL}.
#' @param sel.lower a number to define a lower Fourier period (or the closest available) for the selection
#' of a band of wavelets for the reconstruction.\cr (Only effective if \code{sel.period = NULL}.)
#'
#' Default: \code{NULL}.
#' @param sel.upper a number to define an upper Fourier period (or the closest available) for the selection
#' of a band of wavelets for the reconstruction.\cr (Only effective if \code{sel.period = NULL}.)
#'
#' Default: \code{NULL}.
#' @param rescale Shall the reconstructed series finally be rescaled to attain the original (detrended) series' mean and standard deviation? Logical.
#'
#' Default: \code{TRUE}.
#' @param plot.waves Shall reconstruction waves be plotted? Logical.
#'
#' Default: \code{FALSE}.
#' @param plot.rec Shall the reconstructed series (together with the original (detrended) series) be plotted? Logical.
#'
#' Default: \code{TRUE}.
#' @param lty parameter for the plot of original vs. reconstructed series: line type, e.g. \code{1:2}.
#'
#' Default: \code{1}.
#' @param lwd parameter for the plot of original vs. reconstructed series: line width, e.g. \code{1:2}.
#'
#' Default: \code{1}.
#' @param col parameter for the plot of original vs. reconstructed series: color of lines.
#'
#' Default: \code{1:2}.
#' @param ylim numeric vector of length \code{2}, providing the range of vertical coordinates for the plot.
#'
#' Default: \code{NULL}.
#' @param show.legend Include legend into the plot of original vs. reconstructed series? Logical.
#'
#' Default: \code{TRUE}.
#' @param legend.coords coordinates to position the legend (as in function \code{legend}).
#'
#' Default: \code{"topleft"}.
#' @param legend.horiz Set the legend horizontally rather than vertically? Logical.
#'
#' Default: \code{FALSE}.
#' @param legend.text legend text.
#'
#' Default: \code{c("original (detrended)", "reconstructed")}.
#' @param label.time.axis Label the time axis? Logical.
#'
#' Default: \code{TRUE}.
#' @param show.date Show calendar dates? (Effective only if dates are available as row
#' names or by variable
#' \code{date} in the data frame which is analyzed.) Logical.
#'
#' Default: \code{FALSE}.
#' @param date.format the format of calendar date given as a character string,
#' e.g. \code{"\%Y-\%m-\%d"}, or equivalently \code{"\%F"};
#' see \code{strptime} for a list of implemented date conversion specifications.
#' Explicit information given here will overturn any specification
#' stored in \code{WT}. If unspecified, date formatting is attempted according to \code{as.Date}.
#'
#' Default: \code{NULL}.
#' @param date.tz a character string specifying the time zone of calendar date; see \code{strptime}.
#' Explicit information given here will overturn
#' any specification stored in \code{WT}. If unspecified,
#' \code{""} (the local time zone) is used.
#'
#' Default: \code{NULL}.
#' @param timelab Time axis label.
#'
#' Default: \code{"index"}; in case of a calendar axis: \code{"calendar date"}.
#' @param timetck length of tick marks on the time axis as a fraction of the smaller of the
#' width or height of the plotting region; see \code{par}.
#' If \code{timetck >= 0.5}, \code{timetck} is interpreted as a fraction of the length of
#' the time axis, so if \code{timetck = 1} (and \code{timetcl = NULL}), vertical grid lines
#' will be drawn. \cr Setting \code{timetck = NA} is to use \code{timetcl = -0.5} (which is the R
#' default setting of \code{tck} and \code{tcl}).
#'
#' Default here: \code{0.02}.
#' @param timetcl length of tick marks on the time axis as a fraction of the height of a line
#' of text; see \code{par}. With \code{timetcl = -0.5} (which is the R default setting of
#' \code{tcl}), ticks will be drawn outward.
#'
#' Default here: \code{0.5}.
#' @param spec.time.axisa list of tick mark and label specifications for individualized time axis labeling
#' (only effective if \code{label.time.axis = TRUE}):
#'
#'\itemize{
#'  \item[\code{at}:] locations of tick marks (when \code{NULL}, default plotting will be applied).
#'  Valid tick marks can be provided as numerical values or as dates. Dates are used only in the case
#'  \code{show.date = TRUE}, however, and date formats should conform to \code{as.Date} or the
#'  format given in \code{date.format}. \cr
#'  Default: \code{NULL}.
#'  \item[\code{labels}:] either a logical value specifying whether annotations at the tick marks
#'  are the tick marks themselves, or any vector of labels. If \code{labels} is non-logical,
#'  \code{at} should be of same length. \cr
#'  Default: \code{TRUE}.
#'  \item[\code{las}:] the style of axis labels, see \code{par}. \cr
#'  Default: \code{1} (always horizontal).
#'  \item[\code{hadj}:] adjustment of labels horizontal to the reading direction, see \code{axis}. \cr
#'  Default: \code{NA} (centering is used).
#'  \item[\code{padj}:] adjustment of labels perpendicular to the reading direction (this can be a vector
#'  of adjustments for each label),
#'  see \code{axis}. \cr
#'  Default: \code{NA} (centering is used).
#'}
#' Mismatches will result in a reset to default plotting.
#' @param main.waves an overall title for the plot of reconstruction waves.
#'
#' Default: \code{NULL}.
#' @param main.rec an overall title for the plot of original vs. reconstructed series.
#'
#' Default: \code{NULL}.
#' @param main an overall title for both plots.
#'
#' Default: \code{NULL}.
#' @param lwd.axis line width of axes.
#'
#' Default: \code{1}.
#' @param verbose Print verbose output on the screen? Logical.
#'
#' Default: \code{TRUE}.
#'
#' @value A list of class \code{reconstruct} with the following elements:
#' %%%%%%%%%%%%%%%%%
#'
#' \item{series}{a data frame building on \code{WT$series} with the following columns:
#'     \tabular{rlll}{
#'       \tab date             \tab : \tab the calendar date (if available as column\cr
#'                                                            \tab                  \tab   \tab in WT$series) \cr
#'       \tab <x>              \tab : \tab series <x>, with original name retained \cr
#'       \tab                  \tab : \tab (detrended, if \code{loess.span != 0}) \cr
#'       \tab <x>.trend        \tab : \tab the trend series (if \code{loess.span != 0}) \cr
#'       \tab <x>.r            \tab : \tab the reconstructed (detrended) series
#'     }
#'
#'   Row names are taken over from WT$series, and so are dates if given as row names.
#'   If \code{WT} is of class \code{analyze.coherency}, the second series in the coherency analysis is retained;
#'   if \code{loess.span != 0}, the second series is retained in the detrended version, and the trend is retained as well.
#' }
#' %%%%%%%%%%%%%%%%%
#' \item{rec.waves}{data frame of scaled waves used for reconstruction}
#' %%%%%%%%%%%%%%%%%
#' \item{loess.span}{parameter \code{alpha} in \code{loess} controlling the degree of time series smoothing,
#'   if the time series was detrended; no detrending if \code{loess.span = 0}.}
#' %%%%%%%%%%%%%%%%%
#' \item{lvl}{minimum level of wavelet power for waves (wave segments) to be included
#'   in the reconstruction}
#' \item{only.coi}{Was the influence of edge effects excluded? I.e. was the cone of influence used only?}
#' \item{only.sig}{Was wavelet power significance used in reconstruction?}
#' \item{siglvl}{level of wavelet power significance}
#' \item{only.ridge}{Was the wavelet power ridge used only?}
#' %%%%%%%%%%%%%%%%%
#' \item{rnum.used}{the vector of Fourier period numbers used for reconstruction }
#' %%%%%%%%%%%%%%%%%
#' \item{rescale}{Was the reconstructed series rescaled according to the mean and standard deviation
#'   taken from the original (detrended) series?}
#' %%%%%%%%%%%%%%%%%
#' \item{dt}{time resolution, i.e. sampling resolution in the time domain, \code{1/dt} = number of observations per time unit}
#' \item{dj}{frequency resolution, i.e. sampling resolution in the frequency domain, \code{1/dj} = number of suboctaves (voices per octave)}
#' %%%%%%%%%%%%%%%%%
#' \item{Period}{the Fourier periods
#'   (measured in time units determined by \code{dt}, see the explanations concerning \code{dt})}
#' \item{Scale}{the scales (the Fourier periods divided by the Fourier factor)}
#' %%%%%%%%%%%%%%%%%
#' \item{nc}{number of columns = number of observations = number of observation epochs; "epoch" meaning point in time}
#' \item{nr}{number of rows = number of scales (Fourier periods)}
#' %%%%%%%%%%%%%%%%%
#' \item{axis.1}{tick levels corresponding to the time steps used for (cross-)wavelet transformation: \code{1, 1+dt, 1+2dt, ...}.
#'   The default time axis in plot functions provided by \code{WaveletComp} is determined by observation epochs, however; "epoch" meaning point in time. }
#' \item{axis.2}{tick levels corresponding to the log of Fourier periods: \code{log2(Period)}. This determines the period axis in plot functions provided by \code{WaveletComp}.}
#' %%%%%%%%%%%%%%%%%
#' \item{date.format}{the format of calendar date (if available)}
#' \item{date.tz}{the time zone of calendar date (if available)}
#'
#'
#' @author CoFES. Credits are also due toAngi Roesch and Harald Schmidbauer.
#'
#' @references
#' Carmona R., Hwang W.-L., and Torresani B., 1998.
#' Practical Time Frequency Analysis. Gabor and Wavelet Transforms with an Implementation in S.
#' Academic Press, San Diego.
#'
#' Liu Y., Liang X.S., and Weisberg R.H., 2007.
#' Rectification of the Bias in the Wavelet Power Spectrum.
#' Journal of Atmospheric and Oceanic Technology 24, 2093--2102.
#'
#' Torrence C., and Compo G.P., 1998.
#' A practical guide to wavelet analysis.
#' Bulletin of the American Meteorological Society 79 (1), 61--78.
#'
#' @seealso \code{\link{CoFESWave.Transform}
#'
#' @examples \dontrun{
#' ## The following example is adopted from Liu et al., 2007:
#'
#' series.length = 6*128*24
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
#' my.data <- data.frame(x = x)
#'
#' ## Computation of wavelet power:
#' ## a natural choice of 'dt' in the case of hourly data is 'dt = 1/24',
#' ## resulting in one time unit equaling one day.
#' ## This is also the time unit in which periods are measured.
#' my.w <- analyze.wavelet(my.data, "x",
#'                         loess.span = 0,
#'                         dt = 1/24, dj = 1/20,
#'                         lowerPeriod = 1/4,
#'                         make.pval = TRUE, n.sim = 10)
#'
#' ## Plot of wavelet power spectrum (with equidistant color breakpoints):
#' wt.image(my.w, color.key = "interval",
#'          legend.params = list(lab = "wavelet power levels"),
#'          periodlab = "period (days)")
#'
#' ## Reconstruction of the time series,
#' ## including significant components only:
#' reconstruct(my.w)
#'
#' ## The same reconstruction, but showing wave components first:
#' reconstruct(my.w, plot.waves = TRUE)
#'
#' ## Reconstruction, including all components whether significant or not:
#' reconstruct(my.w, only.sig = FALSE)
#'
#' ## Reconstruction, including significant components,
#' ## but selected periods only (e.g. ignoring period 8):
#' reconstruct(my.w, sel.period = c(1,32,128))
#'
#' ## Reconstruction, including significant components,
#' ## but the ridge only:
#' reconstruct(my.w, only.ridge = TRUE)
#'
#' ## Alternate styles of the time axis:
#'
#' ## The plot with time elapsed in days, starting from 0 and proceeding
#' ## in steps of 50 days (50*24 hours),
#' ## instead of the (default) time index:
#' index.ticks  <- seq(1, series.length, by = 50*24)
#' index.labels <- (index.ticks-1)/24
#'
#' ## Insert your specification of time axis:
#' reconstruct(my.w, only.ridge = TRUE,
#'             timelab = "time elapsed (days)",
#'             spec.time.axis = list(at = index.ticks, labels = index.labels))
#'
#' ## See the periods involved:
#' my.rec <- reconstruct(my.w, only.ridge = TRUE)
#' print(my.rec$Period[my.rec$rnum.used])
#'
#' ## The original and reconstructed time series can be retrieved:
#' plot(my.rec$series$x, type = "l", xlab = "index", ylab = "")
#' lines(my.rec$series$x.r, col="red")
#' legend("topleft", legend = c("original","reconstructed"),
#'        lty = 1, col = c("black","red"))
#'
#' ## Please see also the examples in our guide booklet,
#' ## URL http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf.
#'
#' }
CoFESWave.reconstruct <- function (WT, my.series = 1,
                                   lvl = 0, only.ridge = FALSE, only.sig = TRUE,
                                   siglvl = 0.05, only.coi = FALSE,
                                   sel.period = NULL, sel.lower = NULL, sel.upper = NULL,
                                   rescale = TRUE, plot.waves = FALSE, plot.rec = TRUE,
                                   lty = 1, lwd = 1, col = 1:2, ylim = NULL, show.legend = TRUE,
                                   legend.coords = "topleft", legend.horiz = FALSE, legend.text = NULL,
                                   label.time.axis = TRUE, show.date = FALSE, date.format = NULL,
                                   date.tz = NULL, timelab = NULL, timetck = 0.02, timetcl = 0.5,
                                   spec.time.axis = list(at = NULL, labels = TRUE, las = 1,
                                                         hadj = NA, padj = NA), main.waves = NULL, main.rec = NULL,
                                   main = NULL, lwd.axis = 1, verbose = TRUE)
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

  ## 0.2. Series

  ## 0.2.1. Data Preparation
  series.data = WT$series
  if (class(WT) == "CoFESWave.Transform") {
    out("Your input object class is 'CoFESWave.Transform'...\n")
    my.series = ifelse(names(series.data)[1] == "date", names(series.data)[2],
                       names(series.data)[1])
    Wave = WT$Wave
    Power = WT$Power
    Power.pval = WT$Power.pval
    Ridge = WT$Ridge
  }
  if (class(WT) == "analyze.coherency") {
    out("Your input object class is 'analyze.coherency'...\n")
    if (is.numeric(my.series)) {
      if (!is.element(my.series, c(1, 2))) {
        stop("Please choose either series number 1 or 2!")
      }
      my.series = ifelse(names(series.data)[1] == "date",
                         names(series.data)[my.series + 1], names(series.data)[my.series])
    }
    ind = which(names(series.data) == my.series)
    which.series.num = ifelse(names(series.data)[1] == "date",
                              ind - 1, ind)
    if (!is.element(which.series.num, c(1, 2))) {
      stop("Your series name is not available, please check!")
    }
    if (which.series.num == 1) {
      Wave = WT$Wave.x
      Power = WT$Power.x
      Power.pval = WT$Power.x.pval
      Ridge = WT$Ridge.x
    }
    if (which.series.num == 2) {
      Wave = WT$Wave.y
      Power = WT$Power.y
      Power.pval = WT$Power.y.pval
      Ridge = WT$Ridge.y
    }
  }

  out(paste("Your time series '", my.series, "' will be reconstructed...",
            sep = ""), "\n")
  out("Starting the reconstruction process...\n")

  ## 0.2.2. Data Preparation continued
  nc = WT$nc
  nr = WT$nr
  dt = WT$dt
  dj = WT$dj
  Scale = WT$Scale
  Period = WT$Period
  loess.span = WT$loess.span

  rec.waves = matrix(0, nrow = nr, ncol = nc)
  for (s.ind in seq_len(nr)) {
    rec.waves[s.ind, ] = (Re(Wave[s.ind, ])/sqrt(Scale[s.ind])) *
      dj * sqrt(dt)/(pi^(-1/4) * 0.776)
  }
  rec.waves = rec.waves * (Power >= lvl)
  comment.lvl = paste("minimum power level: ", lvl, sep = "")

  if (only.ridge == T) {
    rec.waves = rec.waves * Ridge
    rec.waves[Ridge == 0] = NA
  }
  comment.ridge = paste("only ridge: ", only.ridge, sep = "")

  if (only.sig == T) {
    if (!is.null(Power.pval)) {
      rec.waves = rec.waves * (Power.pval < siglvl)
      rec.waves[Power.pval >= siglvl] = NA
    }
  }
  if (only.sig == F) {
    siglvl = NA
  }
  comment.sig = paste("significance level: ", siglvl, sep = "")

  if (only.coi == T) {
    for (i in (1:nc)) {
      for (s.ind in seq_len(nr)) {
        if (Scale[s.ind] > 2^WT$coi.2[i]) {
          rec.waves[s.ind, i] = NA
        }
      }
    }
  }
  comment.coi = paste("only coi: ", only.coi, sep = "")

  rnum.used = which(rowSums(rec.waves, na.rm = T) != 0)
  comment.periods = "period: all relevant"
  if (length(sel.period) != 0) {
    sel.rnum = numeric()
    for (i in (1:length(sel.period))) {
      sel.rnum = union(sel.rnum, which(abs(Period - sel.period[i]) ==
                                         min(abs(Period - sel.period[i]))))
    }
    rec.waves = rec.waves[sel.rnum, ]
    comment.periods = paste("period: ", paste(as.character(round(Period[sel.rnum],
                                                                 1)), collapse = ", "), sep = "")
    if (length(sel.rnum) == 1) {
      rec.waves = t(rec.waves)
    }
    rnum.used = intersect(rnum.used, sel.rnum)
  }

  if (length(sel.period) == 0 & ((length(sel.lower) != 0) |
                                 (length(sel.upper) != 0))) {
    if (length(sel.lower) == 0) {
      sel.lower = min(Period)
    }
    if (length(sel.upper) == 0) {
      sel.upper = max(Period)
    }
    if (sel.lower > sel.upper) {
      sel.lower.h = sel.lower
      sel.lower = sel.upper
      sel.upper = sel.lower.h
    }
    sel.rnum = which(((Period >= sel.lower) & (Period <=
                                                 sel.upper)))
    rec.waves = rec.waves[sel.rnum, ]
    sel.period.band = Period[sel.rnum]
    sel.period.range = as.character(round(sel.period.band,
                                          1))
    if (length(sel.rnum) > 1) {
      sel.period.range = paste(round(range(sel.period.band),
                                     1), collapse = " - ")
    }
    if (length(sel.rnum) == 1) {
      rec.waves = t(rec.waves)
    }
    comment.periods = paste("period: ", sel.period.range,
                            sep = "")
    rnum.used = intersect(rnum.used, sel.rnum)
  }

  ### End of 0. Input Check ###



  #### 1. Output (based on paramters) ####
  x.r = colSums(rec.waves, na.rm = T)
  x = series.data[[my.series]]

  ## 1.1. Output
  if (rescale == T) {
    x.r = (x.r - mean(x.r)) * sd(x)/sd(x.r) + mean(x)
  }
  if (is.null(date.format)) {
    date.format = WT$date.format
  }
  if (is.null(date.tz)) {
    date.tz = ifelse(is.null(WT$date.tz), "", WT$date.tz)
  }
  if (!is.list(spec.time.axis))
    spec.time.axis = list()
  if (is.null(spec.time.axis$at))
    spec.time.axis$at = NULL
  if (is.null(spec.time.axis$labels))
    spec.time.axis$labels = T
  if (is.null(spec.time.axis$las))
    spec.time.axis$las = 1
  if (is.null(spec.time.axis$hadj))
    spec.time.axis$hadj = NA
  if (is.null(spec.time.axis$padj))
    spec.time.axis$padj = NA

  time.axis.warning = F
  time.axis.warning.na = F
  chronology.warning = F
  if ((!is.null(spec.time.axis$at)) & (label.time.axis == F))
    warning("\nPlease set label.time.axis = TRUE to make time axis specification effective.",
            immediate. = TRUE)

  ## 1.2. label.time.axis
  if (label.time.axis == T) {
    time.axis.default = (is.null(spec.time.axis$at))
    if (show.date == T) {
      if (is.element("date", names(series.data))) {
        my.date = series.data$date
      }
      else {
        my.date = rownames(series.data)
      }
      if (is.null(date.format)) {
        chronology.warning = inherits(try(as.Date(my.date,
                                                  tz = date.tz), silent = T), "try-error")
        if (!chronology.warning) {
          my.date = as.Date(my.date, tz = date.tz)
          chronology.warning = ifelse(sum(is.na(my.date)) >
                                        0, TRUE, sum(diff(my.date, tz = date.tz) <
                                                       0) > 0)
        }
      }
      if (!is.null(date.format)) {
        chronology.warning = inherits(try(as.POSIXct(my.date,
                                                     format = date.format, tz = date.tz), silent = T),
                                      "try-error")
        if (!chronology.warning) {
          my.date = as.POSIXct(my.date, format = date.format,
                               tz = date.tz)
          chronology.warning = ifelse(sum(is.na(my.date)) >
                                        0, TRUE, sum(diff(my.date, tz = date.tz) <
                                                       0) > 0)
        }
      }
      if (chronology.warning) {
        show.date = F
        time.axis.default = T
        timelab = "index"
      }
    }
    if (!time.axis.default) {
      if (show.date == F) {
        time.axis.warning = (!is.numeric(spec.time.axis$at))
      }
      if (show.date == T) {
        if (is.null(date.format)) {
          time.axis.warning = inherits(try(as.Date(spec.time.axis$at,
                                                   tz = date.tz), silent = T), "try-error")
          if (!time.axis.warning) {
            time.axis.warning = (sum(!is.na(as.Date(spec.time.axis$at,
                                                    tz = date.tz))) == 0)
          }
        }
        if (!is.null(date.format)) {
          time.axis.warning = inherits(try(as.POSIXct(spec.time.axis$at,
                                                      format = date.format, tz = date.tz), silent = T),
                                       "try-error")
          if (!time.axis.warning) {
            time.axis.warning = (sum(!is.na(as.POSIXct(spec.time.axis$at,
                                                       format = date.format, tz = date.tz))) ==
                                   0)
          }
        }
      }
      if (!time.axis.warning) {
        if (is.logical(spec.time.axis$labels)) {
          time.axis.warning = (ifelse(length(spec.time.axis$labels) !=
                                        1, TRUE, is.na(spec.time.axis$labels)))
        }
        if (!is.logical(spec.time.axis$labels)) {
          time.axis.warning = (length(spec.time.axis$labels) !=
                                 length(spec.time.axis$at))
        }
      }
    }
    time.axis.default = (time.axis.default | time.axis.warning)
    if ((is.null(timelab) & time.axis.default) | (!is.null(timelab) &
                                                  time.axis.warning)) {
      timelab = ifelse(show.date, "calendar date", "index")
    }
  }

  ## 1.3. plot.waves
  if (plot.waves == T) {
    out("Reconstruction waves are being plotted...\n")
    if (is.null(main) == F) {
      main.waves = main
    }
    range.rec.waves = range(rec.waves, na.rm = T)
    matplot(1:nc, t(rec.waves), type = "l", ylim = range.rec.waves,
            main = main.waves, sub = paste(comment.lvl, ", ",
                                           comment.sig, ", ", comment.coi, ", ", comment.ridge,
                                           ", ", comment.periods, sep = ""), xaxs = "i",
            xaxt = "n", xlab = "", ylab = "")
    if (label.time.axis == T) {
      if (show.date == F) {
        if (time.axis.default) {
          A.1 = axis(1, lwd = lwd.axis, labels = NA,
                     tck = timetck, tcl = timetcl)
          mtext(A.1, side = 1, at = A.1, line = par()$mgp[2] -
                  0.5, font = par()$font.axis, cex = par()$cex.axis)
        }
        if (!time.axis.default) {
          time.tick = spec.time.axis$at
          time.tick.label = spec.time.axis$labels
          which.na = which(is.na(time.tick))
          if (length(which.na) > 0) {
            time.axis.warning.na = T
          }
          axis(1, lwd = lwd.axis, at = time.tick, labels = time.tick.label,
               tck = timetck, tcl = timetcl, las = spec.time.axis$las,
               hadj = spec.time.axis$hadj, padj = spec.time.axis$padj,
               mgp = par()$mgp - c(0, 0.5, 0), font = par()$font.axis,
               cex.axis = par()$cex.axis)
        }
        mtext(timelab, side = 1, line = par()$mgp[1] -
                1, font = par()$font.lab, cex = par()$cex.lab)
      }
      if (show.date == T) {
        par(new = TRUE)
        if (time.axis.default) {
          plot(my.date, seq(range.rec.waves[1], range.rec.waves[2],
                            length.out = nc), ylim = range.rec.waves,
               type = "n", xaxs = "i", yaxt = "n", xlab = "",
               ylab = "", lwd = lwd.axis, tck = timetck,
               tcl = timetcl, mgp = par()$mgp - c(0, 0.5,
                                                  0), font = par()$font.axis, cex.axis = par()$cex.axis)
        }
        if (!time.axis.default) {
          plot(my.date, seq(range.rec.waves[1], range.rec.waves[2],
                            length.out = nc), ylim = range.rec.waves,
               type = "n", xaxs = "i", xaxt = "n", yaxt = "n",
               xlab = "", ylab = "")
          if (is.null(date.format)) {
            time.tick = as.Date(spec.time.axis$at, tz = date.tz)
          }
          if (!is.null(date.format)) {
            time.tick = as.POSIXct(spec.time.axis$at,
                                   format = date.format, tz = date.tz)
          }
          time.tick.label = spec.time.axis$labels
          which.na = which(is.na(time.tick))
          if (length(which.na) > 0) {
            time.axis.warning.na = T
          }
          if (is.logical(time.tick.label)) {
            if (time.tick.label == T) {
              time.tick.label = time.tick
            }
          }
          axis(1, lwd = lwd.axis, at = time.tick, labels = time.tick.label,
               tck = timetck, tcl = timetcl, las = spec.time.axis$las,
               hadj = spec.time.axis$hadj, padj = spec.time.axis$padj,
               mgp = par()$mgp - c(0, 0.5, 0), font = par()$font.axis,
               cex.axis = par()$cex.axis)
        }
        mtext(timelab, side = 1, line = par()$mgp[1] -
                1, font = par()$font.lab, cex = par()$cex.lab)
      }
    }
  }

  ## 1.4. Output cont. (1)
  x.rec = cbind(series.data, x.r = x.r)
  my.rec.series = paste(my.series, ".r", sep = "")
  colnames(x.rec) = c(colnames(series.data), my.rec.series)
  rownames(x.rec) = rownames(series.data)
  if (plot.waves & plot.rec) {
    par(ask = T)
  }

  ## 1.5. plot.rec
  if (plot.rec == T) {
    out("Original (detrended) and reconstructed series are being plotted...\n")
    if (is.null(main) == F) {
      main.rec = main
    }
    range.rec = range(x.rec[, c(my.series, my.rec.series)],
                      na.rm = T)
    if (is.null(ylim)) {
      ylim = range.rec
    }
    matplot(1:nc, x.rec[, c(my.series, my.rec.series)], type = "l",
            ylim = ylim, lty = lty, lwd = lwd, col = col, main = main.rec,
            sub = paste(comment.lvl, ", ", comment.sig, ", ",
                        comment.coi, ", ", comment.ridge, ", ", comment.periods,
                        sep = ""), xaxs = "i", xaxt = "n", xlab = "",
            ylab = "")
    par(ask = F)
    if (show.legend == T) {
      if (is.null(legend.text)) {
        legend.text = c(paste("original", ifelse(loess.span !=
                                                   0, paste(" (detrended, span: ", loess.span,
                                                            ")", sep = ""), ""), sep = ""), "reconstructed")
      }
      legend(legend.coords, horiz = legend.horiz, legend = legend.text,
             lty = lty, lwd = lwd, col = col, text.font = par()$font.lab,
             cex = par()$cex.lab)
    }
    if (label.time.axis == T) {
      if (show.date == F) {
        if (time.axis.default) {
          A.1 = axis(1, lwd = lwd.axis, labels = NA,
                     tck = timetck, tcl = timetcl)
          mtext(A.1, side = 1, at = A.1, line = par()$mgp[2] -
                  0.5, font = par()$font.axis, cex = par()$cex.axis)
        }
        if (!time.axis.default) {
          time.tick = spec.time.axis$at
          time.tick.label = spec.time.axis$labels
          which.na = which(is.na(time.tick))
          if (length(which.na) > 0) {
            time.axis.warning.na = T
          }
          axis(1, lwd = lwd.axis, at = time.tick, labels = time.tick.label,
               tck = timetck, tcl = timetcl, las = spec.time.axis$las,
               hadj = spec.time.axis$hadj, padj = spec.time.axis$padj,
               mgp = par()$mgp - c(0, 0.5, 0), font = par()$font.axis,
               cex.axis = par()$cex.axis)
        }
        mtext(timelab, side = 1, line = par()$mgp[1] -
                1, font = par()$font.lab, cex = par()$cex.lab)
      }
      if (show.date == T) {
        par(new = TRUE)
        if (time.axis.default) {
          plot(my.date, seq(ylim[1], ylim[2], length.out = nc),
               ylim = ylim, type = "n", xaxs = "i", yaxt = "n",
               xlab = "", ylab = "", lwd = lwd.axis, tck = timetck,
               tcl = timetcl, mgp = par()$mgp - c(0, 0.5,
                                                  0), font = par()$font.axis, cex.axis = par()$cex.axis)
        }
        if (!time.axis.default) {
          plot(my.date, seq(ylim[1], ylim[2], length.out = nc),
               ylim = ylim, type = "n", xaxs = "i", xaxt = "n",
               yaxt = "n", xlab = "", ylab = "")
          if (is.null(date.format)) {
            time.tick = as.Date(spec.time.axis$at, tz = date.tz)
          }
          if (!is.null(date.format)) {
            time.tick = as.POSIXct(spec.time.axis$at,
                                   format = date.format, tz = date.tz)
          }
          time.tick.label = spec.time.axis$labels
          which.na = which(is.na(time.tick))
          if (length(which.na) > 0) {
            time.axis.warning.na = T
          }
          if (is.logical(time.tick.label)) {
            if (time.tick.label == T) {
              time.tick.label = time.tick
            }
          }
          axis(1, lwd = lwd.axis, at = time.tick, labels = time.tick.label,
               tck = timetck, tcl = timetcl, las = spec.time.axis$las,
               hadj = spec.time.axis$hadj, padj = spec.time.axis$padj,
               mgp = par()$mgp - c(0, 0.5, 0), font = par()$font.axis,
               cex.axis = par()$cex.axis)
        }
        mtext(timelab, side = 1, line = par()$mgp[1] -
                1, font = par()$font.lab, cex = par()$cex.lab)
      }
    }
  }

  ## 1.6. Output cont. (2)
  if (chronology.warning) {
    warning("\nPlease check your calendar dates, format and time zone: dates may not be in an unambiguous format or chronological. The default numerical axis was used instead.")
  }
  if (time.axis.warning == T) {
    warning("\nPlease check your time axis specifications. Default settings were used.")
  }
  if (time.axis.warning.na == T) {
    warning("\nNAs were produced with your time axis specifications.")
  }

  ### End of 1. Output (based on paramters) ###



  #### 2. Output ####
  output <- list(series = x.rec, rec.waves = rec.waves,
                 loess.span = loess.span, lvl = lvl,
                 only.coi = only.coi, only.sig = only.sig,
                 siglvl = siglvl, only.ridge = only.ridge,
                 rnum.used = rnum.used, rescale = rescale,
                 dt = dt, dj = dj, Period = Period, Scale = Scale,
                 nc = nc, nr = nr,
                 axis.1 = WT$axis.1, axis.2 = WT$axis.2,
                 date.format = date.format, date.tz = date.tz)
  class(output) = "CoFESWave.reconstruct"
  # out("Class attributes are accessible through following names:\n")
  # out(names(output), "\n")
  return(invisible(output))
}
