#' plot_CoFESWaveWCO
#'
#' @description To plot (absolute value) of wavelet coherency or
#' multiple or partial wavelet coherencies of several series, coi (and levels of sigificance)
#' @param Coherency.Object output of function CoFEScoherency or CoFESmpcoherency
#' @param pE a constant to enhance the quality of picture
#' @param horizons.label Label the horizontal axis with labels (instead of numbers)? Logical.
#'
#' Default: \code{TRUE}.
#' @param low.FP lower periods used in the frequency band
#' @param up.FP upper periods used in the frequency band
#'
#' @author CoFES. Credits are also due to L. Aguiar-Conraria and M.J. Soares.
#'
#' @examples
#' ## more work
plot_CoFESWaveWCO <- function(Coherency.Object,
                              pE=5,
                              horizons.label = TRUE,
                              low.FP = 32, up.FP = 128)
{ ######################
  # Coherency plot
  ######################


  # make test and debug easier.
  # Coherency.Object = testpco
  # horizons.label = TRUE
  # pE=5
  # low.FP = 32; up.FP = 128



  #### 0. Data Preparation ####

  ## 0.0. extract items for plot
  name_class <- class(Coherency.Object)
  date0 <- Coherency.Object$date
  periods <- Coherency.Object$periods
  coi <- Coherency.Object$coi
  phaseDif <- Coherency.Object$phase.dif

  if (name_class == "CoFEScoherency") {
    data.toplot <- Coherency.Object$rwco
    string_title_co <- "Wavelet Coherency"
  }
  if (name_class == "CoFESmpcoherency") {
    data.toplot <- abs(Coherency.Object$mpwco)
    string_title_co <- "Partial Wavelet Coherency"
  }


  ## 0.1. process extracted items
  periods<-log2(periods)
  coi<-log2(coi)

  C<-(t(data.toplot))^pE
  range.C <- range(C)

  min.lim <- round(range.C[1])
  max.lim <- round(range.C[2])
  times <- 1:ncol(data.toplot)


  ## 0.2. Date
  signal_date <- (!is.null(date0) & (length(date0)!=0))
  if (signal_date) {
    date = date0
  } else {
    date = 1:tail(times,1)
  }


  ## 0.3. basic definition
  n.len <- length(Coherency.Object$coi) # we can output number of obs in coherecy functions, for now I use this.
  n.break <- 10
  len_break <- round(n.len / n.break)

  if (horizons.label) {
    scale_labraw <- c("0", "Daily", "Weekly", "Biweekly", "Monthly","Bimonthly", "Quarterly",
                      "Half_year", "Annual", "Bi_annual", "Four_year", "Eight_year", "Sixteen_year")
  } else {
    scale_labraw <- c(0, 2^c(1:20))
  }


  ## 0.4. plot: scale ticker
  levels    <- seq(from = 1, to = floor(log(n.len,2)), 1)
  scale_at0 <- c(1, 2^levels)
  scale_at  <- log2(scale_at0)
  scale_lab <- scale_labraw[1:length(scale_at0)]

  ## 0.5. plot: date ticker
  date_at  <- c(1, len_break*seq(1,(n.break-1),1), n.len)
  date_lab <- date[date_at]
  if (signal_date) { date_lab <- format(date_lab, "%Y-%m") }

  ## 0.6. plot: frequency ticker
  freq_at  <- c(-pi,-(pi/2),0,pi/2,pi)
  freq_lab <- c(expression(-pi),expression(-pi/2),0,expression(pi/2),expression(pi))

  ## 0.7. plot: frequency plot title
  name.band <- paste0(scale_labraw[1+log2(low.FP)], " ~ ", scale_labraw[1+log2(up.FP)]," freq. band")
  ### End of 0. Data Preparation ###

  ## 0.8. plot: color palette
  topl_colors <- colorRampPalette(c("#64d8cb","#26a69a","#5f5fc4", "#283593"))



  #### 1. plot ####

  ## 1.0. basic settings
  oldpar <- par(oma=rep(2.5, 4),
                mfrow = c(1,1),
                mar=c(1,1,1,1))

  # 1.1. Plot of coherency
  layout(matrix(c(1,1,1,1),1,1))
  par(mar=c(6.1,3.1,3.1,2.1))
  fields::image.plot(times,periods,C, zlim = c(min.lim,max.lim),
                     axes = FALSE, main=string_title_co,
                     xlab = " ", ylab = "Periods", col = topl_colors(124))
  polygon(times,coi,border="#5f5fc4", lwd=3)
  axis(side=1, at=date_at,  lab=date_lab,  las=2)
  axis(side=2, at=scale_at, lab=scale_lab, las=1)
  box()

  # 1.2. Plot of phase differences
  plot(phaseDif,type="l",col="#5f5fc4",lwd=1.0,ylim=c(-pi,pi),
       xlab=" ",ylab="Phase Diff",main=name.band,las=3,axes=FALSE,
       panel.first = rect(xleft = rep(date_at[1],4) - max(date_at),
                          ybottom = c(-pi,-pi/2,0,pi/2),
                          xright = rep(tail(date_at,1),4)+ max(date_at),
                          ytop= c(-pi/2,0,pi/2,pi),
                          col=c('gray50', '#8E93CB','#C9CAE6','#E2F6F5','#9EE5DF'),
                          border=NA))
  axis(side=1, at=date_at, lab=date_lab, las=2)
  axis(side=2, at=freq_at, lab=freq_lab, las=1)
  grid()
  box()

  par(oldpar)
  ### End of 1. Coherency and Phase Differences ###


  ### End of the function ###
}
