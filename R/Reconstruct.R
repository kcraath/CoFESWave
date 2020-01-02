#' Title
#'
#' @param WT
#' @param my.series
#' @param lvl
#' @param only.ridge
#' @param only.sig
#' @param siglvl
#' @param only.coi
#' @param sel.period
#' @param sel.lower
#' @param sel.upper
#' @param rescale
#' @param plot.waves
#' @param plot.rec
#' @param lty
#' @param lwd
#' @param col
#' @param ylim
#' @param show.legend
#' @param legend.coords
#' @param legend.horiz
#' @param legend.text
#' @param label.time.axis
#' @param show.date
#' @param date.format
#' @param date.tz
#' @param timelab
#' @param timetck
#' @param timetcl
#' @param spec.time.axis
#' @param main.waves
#' @param main.rec
#' @param main
#' @param lwd.axis
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
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
