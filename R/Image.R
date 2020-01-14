#' CoFESWave.image
#'
#' @title Image plot of the wavelet power spectrum of a single time series
#'
#' @description This function plots the wavelet power spectrum of a single time series,
#' which is provided by an object of
#' class \code{"analyze.wavelet"}, or alternatively of class \code{"analyze.coherency"}.
#' (In the latter case, the series number or name must be specified.)
#' The vertical axis shows the Fourier periods. The horizontal axis shows time step counts, but can
#' be easily transformed into a calendar axis if dates are provided in either row names or as a variable
#' named \code{"date"} in the data frame at hand. Both axes can be relabeled.
#' In particular, an option is given to individualize the period and/or time axis
#' by specifying tick marks and labels.
#'
#' An option is given to raise wavelet power values to any (positive) exponent before
#' plotting in order to accentuate the contrast of the image.
#'
#' The color levels can be defined according to quantiles of values or
#' according to equidistant breakpoints (covering the interval from 0
#'                                       to maximum level), with the number of levels as a further
#' parameter. A user-defined maximum level can be applied.
#' In addition, there is an option to adopt an individual color palette.
#'
#' Further plot design options concern: plot of the cone of
#' influence, plot of wavelet power contour lines at a specified
#' level of significance, plot of power ridges.
#'
#' Finally, there is an option to insert and format a color legend (a
#' right-hand vertical color bar) and to set the plot title.  For
#' further processing of the plot, graphical parameters of plot
#' regions are provided as output.
#'
#' The name and parts of the layout were inspired by a similar function developed by
#' Huidong Tian and Bernard Cazelles (archived R package \code{WaveletCo}).
#' @param WT an object of class \code{"analyze.wavelet"} or \code{"analyze.coherency"}
#' @param my.series In case \code{class(WT) = "analyze.coherency"}:
#' number (\code{1} or \code{2}) or name of the series to be analyzed.
#'
#' Default: \code{1}.
#' @param exponent Exponent applied to values before plotting in order to accentuate
#' the contrast of the image; the exponent should be positive.
#'
#' Default: \code{1}.
#' @param plot.coi Plot cone of influence? Logical.
#'
#' Default: \code{TRUE}.
#' @param plot.contour Plot contour lines to border the area of wavelet power
#' significance? Logical.
#'
#' Default: \code{TRUE}.
#' @param siglvl level of wavelet power significance to be applied to the plot of
#' contour lines.
#'
#' Default: \code{0.1}.
#' @param col.contour color of contour lines. Default: \code{"white"}.
#' @param plot.ridge Plot the wavelet power ridge? Logical.
#'
#' Default: \code{TRUE}.
#' @param lvl minimum level of wavelet power for ridge to be plotted.
#'
#' Default: \code{0}.
#' @param col.ridge ridge color.
#'
#' Default: \code{"black"}.
#' @param color.key How to assign colors to power and coherence levels? Two options:
#'
#' \tabular{rlll}{
#'   \tab  \code{"interval"} or \code{"i"} \tab : \tab equidistant breakpoints \cr
#'   \tab                                  \tab   \tab (from \code{0} through maximum value) \cr
#'   \tab  \code{"quantile"} or \code{"q"} \tab : \tab quantiles
#' }
#'
#' Default: \code{"quantile"}.
#' @param n.levels Number of color levels.
#'
#' Default: \code{100}.
#' @param color.palette Definition of color levels. (The color palette will be assigned
#' to levels in reverse order!)
#'
#' Default: \code{"rainbow(n.levels, start = 0, end = .7)"}.
#' @param maximum.level Maximum plot level of wavelet power considered; only effective
#' in case of equidistant breakpoints (\code{color.key} equaling \code{"i"}).
#'
#' Default: \code{NULL} (referring to maximum level observed).
#' @param useRaster Use a bitmap raster instead of polygons to plot the image? Logical.
#'
#' Default: \code{TRUE}.
#' @param max.contour.segments limit on the number of segments in a single contour line,
#' positive integer.
#'
#' Default: \code{250000} (\code{options(...)} default settings: \code{25000}).
#' @param plot.legend Plot color legend (a vertical bar of colors and breakpoints)? Logical.
#'
#' Default: \code{TRUE}.
#' @param legend.params a list of parameters for the plot of the color legend; parameter
#' values can be set selectively (style in parts adopted from \code{image.plot} in the R
#' package \code{fields} by Douglas Nychka):
#'
#'   \itemize{
#'     \item[\code{width}:] width of legend bar. \cr
#'     Default: \code{1.2}.
#'     \item[\code{shrink}:] a vertical shrinkage factor. \cr
#'     Default: \code{0.9}.
#'     \item[\code{mar}:] right margin of legend bar. \cr
#'     Default: \code{5.1}.
#'     \item[\code{n.ticks}:] number of ticks for labels. \cr
#'     Default: \code{6}.
#'     \item[\code{label.digits}:] digits of labels. \cr
#'     Default: \code{1}.
#'     \item[\code{label.format}:] format of labels. \cr
#'     Default: \code{"f"}.
#'     \item[\code{lab}:] axis label. \cr
#'     Default: \code{NULL}.
#'     \item[\code{lab.line}:] line (in user coordinate units) where to put the axis label. \cr
#'     Default: \code{2.5}.
#'   }
#' @param label.time.axis Label the time axis? Logical.
#'
#' Default: \code{TRUE}.
#' @param show.date Show calendar dates? (Effective only if dates are available as row names
#' or by variable \code{date} in the data frame which is analyzed.) Logical.
#'
#' Default: \code{FALSE}.
#' @param date.format the format of calendar date given as a character string, e.g.
#' \code{"\%Y-\%m-\%d"}, or equivalently \code{"\%F"}; see \code{strptime} for a list of
#' implemented date conversion specifications. Explicit information given here will overturn any
#' specification stored in \code{WT}. If unspecified, date formatting is attempted
#' according to \code{as.Date}.
#'
#' Default: \code{NULL}.
#' @param date.tz a character string specifying the time zone of calendar date;
#' see \code{strptime}. Explicit information given here will overturn
#' any specification stored in \code{WT}. If unspecified, \code{""} (the local time zone) is used.
#'
#' Default: \code{NULL}.
#' @param timelab Time axis label.
#'
#' Default: \code{"index"}; in case of a calendar axis: \code{"calendar date"}.
#' @param timetck length of tick marks on the time axis as a fraction of the smaller of the
#' width or height of the plotting region; see \code{par}. If \code{timetck >= 0.5}, \code{timetck}
#' is interpreted as a fraction of the length of the time axis, so if \code{timetck = 1}
#' (and \code{timetcl = NULL}), vertical grid lines will be drawn. \cr Setting
#' \code{timetck = NA} is to use \code{timetcl = -0.5} (which is the R default setting
#' of \code{tck} and \code{tcl}).
#'
#' Default here: \code{0.02}.
#' @param timetcl length of tick marks on the time axis as a fraction of the height of a
#' line of text; see \code{par}. With \code{timetcl = -0.5} (which is the R default
#' setting of \code{tcl}), ticks will be drawn outward.
#'
#' Default here: \code{0.5}.
#' @param spec.time.axis a list of tick mark and label specifications for individualized
#' time axis labeling (only effective if \code{label.time.axis = TRUE}):
#'
#'   \itemize{
#'     \item[\code{at}:] locations of tick marks (when \code{NULL}, default plotting will be applied).
#'     Valid tick marks can be provided as numerical values or as dates. Dates are used only in the case \code{show.date = TRUE}, however,
#'     and date formats should conform to \code{as.Date} or the format given in \code{date.format}. \cr
#'     Default: \code{NULL}.
#'     \item[\code{labels}:] either a logical value specifying whether annotations at the tick marks are the tick marks themselves,
#'     or any vector of labels. If \code{labels} is non-logical, \code{at} should be of same length. \cr
#'     Default: \code{TRUE}.
#'     \item[\code{las}:] the style of axis labels, see \code{par}. \cr
#'     Default: \code{1} (always horizontal).
#'     \item[\code{hadj}:] adjustment of labels horizontal to the reading direction, see \code{axis}. \cr
#'     Default: \code{NA} (centering is used).
#'     \item[\code{padj}:] adjustment of labels perpendicular to the reading direction (this can be a vector of adjustments for each label),
#'     see \code{axis}. \cr
#'     Default: \code{NA} (centering is used).
#'   }
#' Mismatches will result in a reset to default plotting.
#' @param label.period.axis Label the (Fourier) period axis? Logical.
#'
#' Default: \code{TRUE}.
#' @param periodlab (Fourier) period axis label.
#'
#' Default: \code{"period"}.
#' @param periodtck length of tick marks on the period axis as a fraction of the smaller
#' of the width or height of the plotting region; see \code{par}. If \code{periodtck >= 0.5},
#' \code{periodtck} is interpreted as a fraction of the length of the period axis, so if
#' \code{periodtck = 1} (and \code{periodtcl = NULL}), horizontal grid lines will be drawn.
#' \cr Setting \code{periodtck = NA} is to use \code{periodtcl = -0.5} (which is the R default
#' setting of \code{tck} and \code{tcl}).
#'
#' Default here: \code{0.02}.
#' @param periodtcl length of tick marks on the period axis as a fraction of the height of
#' a line of text; see \code{par}. With \code{periodtcl = -0.5} (which is the R default
#' setting of \code{tcl}) ticks will be drawn outward.
#'
#' Default here: \code{0.5}.
#' @param spec.period.axis a list of tick mark and label specifications for individualized
#' period axis labeling (only effective if \code{label.period.axis = TRUE}):
#'
#'   \itemize{
#'     \item[\code{at}:] locations of tick marks (when \code{NULL}, default plotting will be applied). Valid tick marks can be provided as
#'     numerical and positive values only. \cr
#'     Default: \code{NULL}.
#'     \item[\code{labels}:] either a logical value specifying whether annotations at the tick marks are the tick marks themselves,
#'     or any vector of labels. If \code{labels} is non-logical, \code{at} should be of same length. \cr
#'     Default: \code{TRUE}.
#'     \item[\code{las}:] the style of axis labels, see \code{par}. \cr
#'     Default: \code{1} (always horizontal).
#'     \item[\code{hadj}:] adjustment of labels horizontal to the reading direction, see \code{axis}. \cr
#'     Default: \code{NA} (centering is used).
#'     \item[\code{padj}:] adjustment of labels perpendicular to the reading direction (this can be a vector of adjustments for each label),
#'     see \code{axis}. \cr
#'     Default: \code{NA} (centering is used).
#'   }
#' Mismatches will result in a reset to default plotting.
#' @param main an overall title for the plot.
#'
#' Default: \code{NULL}.
#' @param lwd line width of contour lines and ridge.
#'
#' Default: \code{2}.
#' @param lwd.axis line width of axes (image and legend bar).
#'
#' Default: \code{1}.
#' @param graphics.reset Reset graphical parameters? Logical.
#'
#' Default: \code{TRUE}.
#' @param verbose Print verbose output on the screen? Logical.
#'
#' Default: \code{FALSE}.
#'
#' @return A list of class \code{graphical parameters} with the following elements:
#' \item{op}{original graphical parameters}
#' \item{image.plt}{image plot region}
#' \item{legend.plt}{legend plot region}
#'
#' @references
#' Aguiar-Conraria L., and Soares M.J., 2011.
#' The Continuous Wavelet Transform: A Primer.
#' NIPE Working Paper Series 16/2011.
#'
#' Carmona R., Hwang W.-L., and Torresani B., 1998.
#' Practical Time Frequency Analysis. Gabor and Wavelet Transforms with an Implementation in S.
#' Academic Press, San Diego.
#'
#' Cazelles B., Chavez M., Berteaux, D., Menard F., Vik J.O., Jenouvrier S., and Stenseth N.C., 2008.
#' Wavelet analysis of ecological time series.
#' Oecologia 156, 287--304.
#'
#' Liu Y., Liang X.S., and Weisberg R.H., 2007.
#' Rectification of the Bias in the Wavelet Power Spectrum.
#' Journal of Atmospheric and Oceanic Technology 24, 2093--2102.
#'
#' Tian, H., and Cazelles, B., 2012. \code{WaveletCo}.
#' Available at \url{https://cran.r-project.org/src/contrib/Archive/WaveletCo/}, archived April 2013; accessed July 26, 2013.
#'
#' Torrence C., and Compo G.P., 1998.
#' A practical guide to wavelet analysis.
#' Bulletin of the American Meteorological Society 79 (1), 61--78.
#'
#' @author CoFES. Credits are also due toAngi Roesch, Harald Schmidbauer, Huidong Tian, and Bernard Cazelles.
#'
#' @seealso \code{\link{CoFESWave.Transform}}, \code{\link{CoFESWave.reconstruct}}
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
#' my.wt <- analyze.wavelet(my.data, "x",
#'                          loess.span = 0,
#'                          dt = 1/24, dj = 1/20,
#'                          lowerPeriod = 1/4,
#'                          make.pval = TRUE, n.sim = 10)
#'
#' ## Plot of wavelet power spectrum
#' ## with color breakpoints referring to quantiles:
#' wt.image(my.wt, main = "wavelet power spectrum",
#'          legend.params = list(lab = "wavelet power levels (quantiles)",
#'                               lab.line = 3.5,
#'                               label.digits = 2),
#'          periodlab = "period (days)")
#' ## Note:
#' ## The default time axis shows an index of given points in time,
#' ## which is the count of hours in our example.
#'
#' ## The same plot, but with equidistant color breakpoints:
#' wt.image(my.wt, color.key = "i", main = "wavelet power spectrum",
#'          legend.params = list(lab = "wavelet power levels (equidistant)"),
#'          periodlab = "period (days)")
#'
#' ## Alternative styles of the time axis:
#'
#' ## The plot with time elapsed in days, starting from 0 and proceeding
#' ## in steps of 50 days (50*24 hours),
#' ## instead of the (default) time index:
#' index.ticks  <- seq(1, series.length, by = 50*24)
#' index.labels <- (index.ticks-1)/24
#' ## Insert your specification of the time axis:
#' wt.image(my.wt, color.key = "i", main = "wavelet power spectrum",
#'          legend.params = list(lab = "wavelet power levels (equidistant)"),
#'          periodlab = "period (days)", timelab = "time elapsed (days)",
#'          spec.time.axis = list(at = index.ticks, labels = index.labels))
#'
#' ## The plot with (automatically produced) calendar axis:
#' wt.image(my.wt, color.key = "i", main = "wavelet power spectrum",
#'          legend.params = list(lab = "wavelet power levels (equidistant)"),
#'          periodlab = "period (days)",
#'          show.date = TRUE, date.format = "\%F \%T")
#'
#' ## Individualizing your calendar axis (works with 'show.date = TRUE')...
#' ## How to obtain, for example, monthly date ticks and labels:
#'
#' ## The sequence of tick positions:
#' monthly.ticks <- seq(as.POSIXct("2014-11-01 00:00:00", format = "\%F \%T"),
#'                      as.POSIXct("2016-11-01 00:00:00", format = "\%F \%T"),
#'                      by = "month")
#' ## Observe that the following specification may produce an error:
#' ## 'seq(as.Date("2014-11-01"), as.Date("2016-11-01"), by = "month")'
#' ## Time of the day is missing here!
#'
#' ## The sequence of labels (e.g. information on month and year only):
#' monthly.labels <- strftime(monthly.ticks, format = "\%b \%Y")
#'
#' ## Insert your specification of the time axis:
#' wt.image(my.wt, color.key = "i", main = "wavelet power spectrum",
#'          legend.params = list(lab = "wavelet power levels (equidistant)"),
#'          periodlab = "period (days)",
#'          show.date = TRUE, date.format = "\%F \%T",
#'          spec.time.axis = list(at = monthly.ticks, labels = monthly.labels,
#'                                las = 2))
#' ## Note:
#' ## The monthly ticks specify the midpoints of the colored cells and match
#' ## the location of corresponding (default) time index ticks.
#'
#' ## Furthermore, the plot with an individualized period axis:
#' wt.image(my.wt, color.key = "i", main = "wavelet power spectrum",
#'          legend.params = list(lab = "wavelet power levels (equidistant)"),
#'          periodlab = "period (days)",
#'          show.date = TRUE, date.format = "\%F \%T",
#'          spec.time.axis = list(at = monthly.ticks, labels = monthly.labels,
#'                                las = 2),
#'          spec.period.axis = list(at = c(1,8,32,128)))
#'
#' ## Switching the time axis from index to time elapsed in hours
#' ## (starting from 0, and proceeding in steps of 500 hours),
#' ## and the period axis from days to hours:
#' index.ticks  <- seq(1, series.length, by = 500)
#' index.labels <- index.ticks - 1
#' wt.image(my.wt, color.key = "i", main = "wavelet power spectrum",
#'          legend.params = list(lab = "wavelet power levels (equidistant)"),
#'          timelab = "time elapsed (hours)", periodlab = "period (hours)",
#'          spec.time.axis = list(at = index.ticks, labels = index.labels),
#'          spec.period.axis = list(at = c(1,8,32,128), labels = c(1,8,32,128)*24))
#'
#' ## A plot with different colors:
#' wt.image(my.wt, main = "wavelet power spectrum",
#'          legend.params = list(lab = "wavelet power levels (quantiles)",
#'                               lab.line = 3.5, label.digits = 2),
#'          color.palette = "gray((1:n.levels)/n.levels)", col.ridge = "yellow",
#'          periodlab = "period (days)")
#'
#' ## In the case of monthly (or quarterly) data, the time axis should be
#' ## labeled at equally spaced time points. An example:
#'
#' monthyear <- seq(as.Date("2014-01-01"), as.Date("2018-01-01"),
#'                  by = "month")
#' monthyear <- strftime(monthyear, format = "\%b \%Y")
#'
#' xx <- periodic.series(start.period = 6, length = length(monthyear))
#' xx <- xx + 0.2*rnorm(length(monthyear))
#'
#' plot(xx, type = "l", xlab = "index", ylab = "", xaxs = "i",
#'      main = "monthly series with period of 6 months")
#'
#' monthly.data <- data.frame(date = monthyear, xx = xx)
#'
#' my.wt <- analyze.wavelet(monthly.data, "xx", loess.span = 0,
#'                          dt = 1, dj = 1/250,
#'                          make.pval = TRUE, n.sim = 250)
#' ## Note:
#' ## The natural choice of 'dt' in this example is 'dt = 1',
#' ## resulting in periods measured in months.
#' ## (Setting 'dt = 1/12' would result in periods measured in years.)
#'
#' ## The default wavelet power plot then shows the monthly:
#' wt.image(my.wt, main = "wavelet power spectrum",
#'          periodlab = "period (months)")
#'
#' ## The following plot shows the elapsed time, measured in months:
#' wt.image(my.wt, main = "wavelet power spectrum",
#'          periodlab = "period (months)", timelab = "time elapsed (months)",
#'          spec.time.axis = list(at = 1:length(monthyear),
#'                                labels = (1:length(monthyear))-1))
#'
#' ## In case you prefer the monthyear labels themselves:
#' wt.image(my.wt,  main = "wavelet power spectrum",
#'          periodlab = "period (months)", timelab = "month and year",
#'          spec.time.axis = list(at = 1:length(monthyear), labels = monthyear))
#'
#' ## You may sometimes wish to enhance your plot with additional information.
#' ## There is an option to add further objects to the image plot region,
#' ## by setting 'graphics.reset = FALSE'
#' ## (but recall previous par settings after plotting):
#'
#' op <- par(no.readonly = TRUE)
#' wt.image(my.wt, main = "wavelet power spectrum",
#'          periodlab = "period (months)",
#'          spec.period.axis = list(at = c(2,4,6,8,12)),
#'          spec.time.axis = list(at = 1:length(monthyear),
#'                                labels = substr(monthyear,1,3)),
#'          graphics.reset = FALSE)
#' abline(h = log2(6), lty = 3)
#' abline(v = seq(1, length(monthyear), by = 12), lty = 3)
#' mtext(2014:2018, side = 1,
#'       at = seq(1, length(monthyear), by = 12), line = 2)
#' par(op)
#'
#' ## For further axis plotting options:
#' ## Please see the examples in our guide booklet,
#' ## URL http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf.
#'
#' }
CoFESWave.image <- function (WT, my.series = 1, exponent = 1, plot.coi = TRUE, plot.contour = TRUE,
                             siglvl = 0.1, col.contour = "white", plot.ridge = TRUE, lvl = 0,
                             col.ridge = "black", color.key = "quantile", n.levels = 100,
                             color.palette = "rainbow(n.levels, start = 0, end = .7)",
                             maximum.level = NULL, useRaster = TRUE, max.contour.segments = 250000,
                             plot.legend = TRUE,
                             legend.params = list(width = 1.2, shrink = 0.9,
                                                  mar = 5.1, n.ticks = 6, label.digits = 1, label.format = "f",
                                                  lab = NULL, lab.line = 2.5),
                             label.time.axis = TRUE,
                             show.date = FALSE, date.format = NULL, date.tz = NULL, timelab = NULL,
                             timetck = 0.02, timetcl = 0.5,
                             spec.time.axis = list(at = NULL,
                                                   labels = TRUE, las = 1, hadj = NA, padj = NA),
                             label.period.axis = TRUE,
                             periodlab = NULL, periodtck = 0.02, periodtcl = 0.5,
                             spec.period.axis = list(at = NULL,
                                                     labels = TRUE, las = 1, hadj = NA, padj = NA),
                             main = NULL,
                             lwd = 2, lwd.axis = 1, graphics.reset = TRUE, verbose = FALSE)
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

  ## 0.2. Options
  default.options = options()
  options(max.contour.segments = as.integer(max.contour.segments))

  ## 0.3. Series
  axis.1 <- WT$axis.1
  axis.2 <- WT$axis.2
  series.data = WT$series
  if (exponent <= 0) {
    stop("Please use a positive exponent, or return to default setting (exponent=1)!")
  }
  if (class(WT) == "CoFESWave.Transform") {
    out("Your input object class is 'CoFESWave.Transform'...\n")
    my.series = ifelse(names(series.data)[1] == "date", names(series.data)[2],
                       names(series.data)[1])
    Power = WT$Power^exponent
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
      Power = WT$Power.x^exponent
      Power.pval = WT$Power.x.pval
      Ridge = WT$Ridge.x
    }
    if (which.series.num == 2) {
      Power = WT$Power.y^exponent
      Power.pval = WT$Power.y.pval
      Ridge = WT$Ridge.y
    }
  }
  out(paste("A wavelet power image of your time series '",
            my.series, "' will be plotted...", sep = ""), "\n")

  ### End of 0. Input Check ###



  #### 1. Plot Preparation ####

  ## 1.1. Preparation (based on parameters)
  topl_colors <- colorRampPalette(c("#64d8cb","#26a69a", "#90a4ae","#5f5fc4", "#283593"))
  if ((!is.null(maximum.level)) & (is.element(color.key, c("quantile",
                                                           "q")))) {
    warning("\nPlease set color.key = 'i' to make your maximum level specification effective.",
            immediate. = TRUE)
  }
  if (is.element(color.key, c("interval", "i"))) {
    if (is.null(maximum.level)) {
      maximum.level = max(Power)
    }
    if (maximum.level < max(Power)) {
      stop(paste("... plot can't be produced! Your choice of maximum plot level is smaller than the maximum level observed! Please choose maximum.level larger than ",
                 max(Power), " or return to default setting (maximum.level = NULL)!",
                 sep = ""))
    }
  }
  if (is.element(color.key, c("interval", "i"))) {
    wavelet.levels = seq(from = 0, to = maximum.level, length.out = n.levels +
                           1)
  }
  if (is.element(color.key, c("quantile", "q"))) {
    wavelet.levels = quantile(Power, probs = seq(from = 0,
                                                 to = 1, length.out = n.levels + 1))
  }
  key.cols = rev(eval(parse(text = color.palette)))
  if (!is.list(legend.params))
    legend.params = list()
  if (is.null(legend.params$width))
    legend.params$width = 1.2
  if (is.null(legend.params$shrink))
    legend.params$shrink = 0.9
  if (is.null(legend.params$mar))
    legend.params$mar = ifelse(is.null(legend.params$lab),
                               5.1, 6.1)
  if (is.null(legend.params$n.ticks))
    legend.params$n.ticks = 6
  if (is.null(legend.params$label.digits))
    legend.params$label.digits = 1
  if (is.null(legend.params$label.format))
    legend.params$label.format = "f"
  if (is.null(legend.params$lab.line))
    legend.params$lab.line = 2.5
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
  if (!is.list(spec.period.axis))
    spec.period.axis = list()
  if (is.null(spec.period.axis$at))
    spec.period.axis$at = NULL
  if (is.null(spec.period.axis$labels))
    spec.period.axis$labels = T
  if (is.null(spec.period.axis$las))
    spec.period.axis$las = 1
  if (is.null(spec.period.axis$hadj))
    spec.period.axis$hadj = NA
  if (is.null(spec.period.axis$padj))
    spec.period.axis$padj = NA
  period.axis.warning = F
  period.axis.warning.na = F
  if ((!is.null(spec.period.axis$at)) & (label.period.axis ==
                                         F))
    warning("\nPlease set label.period.axis = TRUE to make period axis specification effective.",
            immediate. = TRUE)
  op = par(no.readonly = TRUE)
  image.plt = par()$plt
  legend.plt = NULL

  ## 1.2. plot.legend
  if (plot.legend == T) {
    legend.plt = par()$plt
    char.size = par()$cin[1]/par()$din[1]
    hoffset = char.size * par()$mar[4]
    legend.width = char.size * legend.params$width
    legend.mar = char.size * legend.params$mar
    legend.plt[2] = 1 - legend.mar
    legend.plt[1] = legend.plt[2] - legend.width
    vmar = (legend.plt[4] - legend.plt[3]) * ((1 - legend.params$shrink)/2)
    legend.plt[4] = legend.plt[4] - vmar
    legend.plt[3] = legend.plt[3] + vmar
    image.plt[2] = min(image.plt[2], legend.plt[1] - hoffset)
    par(plt = legend.plt)
    key.marks = round(seq(from = 0, to = 1, length.out = legend.params$n.ticks) *
                        n.levels)
    key.labels = formatC(as.numeric(wavelet.levels), digits = legend.params$label.digits,
                         format = legend.params$label.format)[key.marks +
                                                                1]
    image(1, seq(from = 0, to = n.levels), matrix(wavelet.levels,
                                                  nrow = 1), col = key.cols, breaks = wavelet.levels,
          useRaster = T, xaxt = "n", yaxt = "n", xlab = "",
          ylab = "")
    axis(4, lwd = lwd.axis, at = key.marks, labels = NA,
         tck = 0.02, tcl = (par()$usr[2] - par()$usr[1]) *
           legend.params$width - 0.04)
    mtext(key.labels, side = 4, at = key.marks, line = 0.5,
          las = 2, font = par()$font.axis, cex = par()$cex.axis)
    text(x = par()$usr[2] + (1.5 + legend.params$lab.line) *
           par()$cxy[1], y = n.levels/2, labels = legend.params$lab,
         xpd = NA, srt = 270, font = par()$font.lab, cex = par()$cex.lab)
    box(lwd = lwd.axis)
    par(new = TRUE, plt = image.plt)
  }

  ### End of 1. Plot Preparation ###



  #### 2. Plot ####
  image(axis.1, axis.2, t(Power), col = key.cols, breaks = wavelet.levels,
        useRaster = useRaster, ylab = "", xlab = "", axes = FALSE,
        main = main)

  ## 2.1. Further Improvement
  if ((plot.contour == T) & (is.null(Power.pval) == F)) {
    contour(axis.1, axis.2, t(Power.pval) < siglvl, levels = 1,
            lwd = lwd, add = TRUE, col = col.contour, drawlabels = FALSE)
  }
  if (plot.ridge == T) {
    Ridge = Ridge * (Power >= lvl)
    contour(axis.1, axis.2, t(Ridge), levels = 1, lwd = lwd,
            add = TRUE, col = col.ridge, drawlabels = FALSE)
  }
  if (plot.coi == T) {
    polygon(WT$coi.1, WT$coi.2, border = NA, col = rgb(1,
                                                       1, 1, 0.5))
  }
  box(lwd = lwd.axis)

  ## 2.2. label.period.axis
  if (label.period.axis == T) {
    period.axis.default = (is.null(spec.period.axis$at))
    if (!period.axis.default) {
      if (is.numeric(spec.period.axis$at)) {
        period.axis.warning = ((sum(spec.period.axis$at <=
                                      0, na.rm = T) > 0) | (sum(!is.na(spec.period.axis$at)) !=
                                                              sum(is.finite(spec.period.axis$at))))
      }
      else {
        period.axis.warning = TRUE
      }
      if (is.logical(spec.period.axis$labels)) {
        period.axis.warning = (period.axis.warning |
                                 ifelse(length(spec.period.axis$labels) != 1,
                                        TRUE, is.na(spec.period.axis$labels)))
      }
      if (!is.logical(spec.period.axis$labels)) {
        period.axis.warning = (period.axis.warning |
                                 (length(spec.period.axis$labels) != length(spec.period.axis$at)))
      }
    }
    period.axis.default = (period.axis.default | period.axis.warning)
    if ((is.null(periodlab)) | (!is.null(periodlab) & period.axis.warning)) {
      periodlab = "period"
    }
    if (period.axis.default) {
      period.tick = unique(trunc(axis.2))
      period.tick[period.tick < log2(WT$Period[1])] = NA
      period.tick = na.omit(period.tick)
      period.tick.label = 2^(period.tick)
      axis(2, lwd = lwd.axis, at = period.tick, labels = NA,
           tck = periodtck, tcl = periodtcl)
      axis(4, lwd = lwd.axis, at = period.tick, labels = NA,
           tck = periodtck, tcl = periodtcl)
      mtext(period.tick.label, side = 2, at = period.tick,
            las = 1, line = par()$mgp[2] - 0.5, font = par()$font.axis,
            cex = par()$cex.axis)
    }
    if (!period.axis.default) {
      period.tick = log2(spec.period.axis$at)
      period.tick[(period.tick < log2(WT$Period[1]))] = NA
      period.tick.label = spec.period.axis$labels
      which.na = which(is.na(period.tick))
      if (length(which.na) > 0) {
        period.axis.warning.na = T
      }
      if (is.logical(period.tick.label)) {
        if (period.tick.label == T) {
          period.tick.label = 2^(period.tick)
        }
      }
      axis(2, lwd = lwd.axis, at = period.tick, labels = period.tick.label,
           tck = periodtck, tcl = periodtcl, las = spec.period.axis$las,
           hadj = spec.period.axis$hadj, padj = spec.period.axis$padj,
           mgp = par()$mgp - c(0, 0.5, 0), font = par()$font.axis,
           cex.axis = par()$cex.axis)
      axis(4, lwd = lwd.axis, at = period.tick, labels = NA,
           tck = periodtck, tcl = periodtcl)
    }
    mtext(periodlab, side = 2, line = par()$mgp[1] - 0.5,
          font = par()$font.lab, cex = par()$cex.lab)
  }

  ## 2.3. label.time.axis
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
  }

  ## 2.4. Further Improvement cont. (1)
  half.increment.2 = (axis.2[2] - axis.2[1])/2
  my.ylim = c(min(axis.2) - half.increment.2, max(axis.2) +
                half.increment.2)
  index.condition = ((show.date == F) | (label.time.axis ==
                                           F))
  if ((WT$dt != 1) & (index.condition)) {
    par(new = TRUE)
    plot(1:WT$nc, seq(min(axis.2), max(axis.2), length.out = WT$nc),
         xlim = c(0.5, WT$nc + 0.5), ylim = my.ylim, type = "n",
         xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "",
         ylab = "")
  }

  ## 2.5. label.time.axis cont.
  if (label.time.axis == T) {
    if ((is.null(timelab) & time.axis.default) | (!is.null(timelab) &
                                                  time.axis.warning)) {
      timelab = ifelse(show.date, "calendar date", "index")
    }
    if (show.date == F) {
      if (time.axis.default) {
        A.1 = axis(1, lwd = lwd.axis, labels = NA, tck = timetck,
                   tcl = timetcl)
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
      mtext(timelab, side = 1, line = par()$mgp[1] - 1,
            font = par()$font.lab, cex = par()$cex.lab)
    }
    if (show.date == T) {
      half.increment = difftime(my.date[2], my.date[1],
                                tz = date.tz)/2
      if (is.null(date.format)) {
        my.xlim = as.Date(range(my.date) - c(half.increment,
                                             -half.increment), tz = date.tz)
      }
      else {
        my.xlim = as.POSIXct(range(my.date) - c(half.increment,
                                                -half.increment), format = date.format, tz = date.tz)
      }
      par(new = TRUE)
      if (time.axis.default) {
        plot(my.date, seq(min(axis.2), max(axis.2), length.out = WT$nc),
             xlim = my.xlim, ylim = my.ylim, type = "n",
             xaxs = "i", yaxs = "i", yaxt = "n", xlab = "",
             ylab = "", lwd = lwd.axis, tck = timetck, tcl = timetcl,
             mgp = par()$mgp - c(0, 0.5, 0), font = par()$font.axis,
             cex.axis = par()$cex.axis)
      }
      if (!time.axis.default) {
        plot(my.date, seq(min(axis.2), max(axis.2), length.out = WT$nc),
             xlim = my.xlim, ylim = my.ylim, type = "n",
             xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
             xlab = "", ylab = "")
        if (is.null(date.format)) {
          time.tick = as.Date(spec.time.axis$at, tz = date.tz)
        }
        if (!is.null(date.format)) {
          time.tick = as.POSIXct(spec.time.axis$at, format = date.format,
                                 tz = date.tz)
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
      mtext(timelab, side = 1, line = par()$mgp[1] - 1,
            font = par()$font.lab, cex = par()$cex.lab)
    }
  }

  ## 2.6. Further Improvement cont. (2)
  options(default.options)
  if (graphics.reset == T) {
    par(op)
  }
  if (period.axis.warning == T) {
    warning("\nPlease check your period axis specifications. Default settings were used.")
  }
  if (period.axis.warning.na == T) {
    warning("\nNAs were produced with your period axis specifications.")
  }
  if (chronology.warning) {
    warning("\nPlease check your calendar dates, format and time zone: dates may not be in an unambiguous format or chronological. The default numerical axis was used instead.")
  }
  if (time.axis.warning == T) {
    warning("\nPlease check your time axis specifications. Default settings were used.")
  }
  if (time.axis.warning.na == T) {
    warning("\nNAs were produced with your time axis specifications.")
  }

  ### End of 2. Plot ###



  #### 3. Output ####
  output = list(op = op,
                image.plt = image.plt,
                legend.plt = legend.plt)
  class(output) = "graphical parameters"
  # out("Class attributes are accessible through following names:\n")
  # out(names(output), "\n")
  return(invisible(output))
}
