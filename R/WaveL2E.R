#' WaveL2E
#'
#' @title Wavelet Analysis with L2E and L2E_chi^2 methods
#'
#' @description This function is based on L2E and L2E_chi^2 thresholding techniques.
#'
#' @param x a series
#' @param date a date series
#' @param block number of observations within a block
#'
#' Default: \code{1}.
#' @param base_plot Plot wavelet power spectrum of x? Logical
#'
#' Default: \code{TRUE}.
#' @param L2E Apply L2E thresholding method? Logical
#'
#' Default: \code{TRUE}.
#' @param Chi_square Apply L2E_Chi_square thresholding method? Logical
#'
#' Default: \code{TRUE}.
#'
#' @return A list of class \code{"WaveL2E"} with elements of different dimensions.
#' %%%%%%%%%%%%%%%%%
#' Here is a detailed list of all elements:
#'#' %%%%%%%%%%%%%%%%%
#' \item{sig}{tmp}
#' \item{w}{tmp}
#' \item{dis0}{tmp}
#' \item{thresh}{threshold of \code{x} with \eqn{L2E} thresholding}
#' \item{qthresh}{threshold of \code{x} with \eqn{L2E(\chi^2)} thresholding}
#' \item{original}{wavelet analysis of \code{x} without thresholding}
#' \item{Ana_Wave}{wavelet analysis of \code{x} with \eqn{L2E} thresholding}
#' \item{Emp_WaveL2E}{tmp}
#' \item{Emp_WaveL2E_MAD}{tmp}
#' \item{recon_L2E}{reconstructed series of \code{x} with \eqn{L2E} thresholding}
#' \item{PTV_L2E}{percentage of total volume of \code{x} with \eqn{L2E} thresholding}
#' \item{PSL_L2E}{percentage of significance area of \code{x} with \eqn{L2E} thresholding}
#' \item{recon_Chi_square}{reconstructed series of \code{x} with \eqn{L2E(\chi^2)} thresholding}
#' \item{PTV_Chi_square}{percentage of total volume of \code{x} with \eqn{L2E(\chi^2)} thresholding}
#' \item{PSL_Chi_square}{percentage of significance area of \code{x} with \eqn{L2E(\chi^2)} thresholding}
#' \item{date}{the corrsponding date series}
#'
#' @author CoFES.
#' @examples
#' ## WaveL2E(x)
WaveL2E <- function(x, date = NULL, block = 1, base_plot = TRUE,
                    L2E = TRUE,
                    Chi_square = TRUE)
{ #Maybe we should add in a static option (if static - length of Data, or we are starting at index =2 so (length - 1))
  #For this one w should also just have a table for the variance and the weight - only one value each on for each block

  # library(latex2exp)
  topl_colors <- colorRampPalette(c("#64d8cb","#26a69a", "#90a4ae","#5f5fc4", "#283593"))

  #date = NULL; block = 1; base_plot = TRUE; L2E = TRUE; Chi_square = TRUE; # make test and debug easier.
  # x = Data[,1]; date = date_Sector;
  # block = 50; base_plot = TRUE;
  # L2E = TRUE; Chi_square = TRUE; # make test and debug easier.



  #### 0. Data Preparation ####

  ## 0.0. Date
  Signal_date = !is.null(date)
  if (Signal_date) {
    Data_raw <- data.frame(date =  date, x = x )
  } else {
    Data_raw <- data.frame(x = x)
  }


  ## 0.1. Base Wavelet Transformation
  ## 0.1.1. transform
  my.w_original <- CoFESWave.Transform(Data_raw, "x",
                                                loess.span = 0.0, make.pval = F, verbose = F ,
                                                n.sim = 10000, upperPeriod = length(x))
  recon_org0 <- CoFESWave.reconstruct(my.w_original, plot.waves = FALSE, lwd = c(1,2),
                                         legend.coords = "topright",plot.rec = FALSE, verbose = F)

  ## 0.1.2. plot before WaveL2E
  if (base_plot) {
    CoFESWave.image(my.w_original, periodlab = " ",timelab = "  " , main = " ",
                          legend.params = list(lab = "wavelet power levels", mar = 5.1, cex = 6, n.ticks = 10),
                          color.key = "quantile", lwd = 2, plot.ridge = FALSE, color.palette = "topl_colors(n.levels)",
                          show.date = Signal_date) #"topl_colors(n.levels)"
    title(latex2exp::TeX(" CWT Power Spectrum"), cex.main = 1.3,
          xlab = " ", ylab = " ", cex.lab = 1.3)
  }

  ## 0.1.3. base transformation for L2E and Chi-square
  Power <- my.w_original$Power
  Time <- ncol(Power)
  nrow_power <- nrow(Power)
  my.w1 <- my.w_original # for L2E
  my.w2 <- my.w_original # for Chi-square


  ## 0.2. Block
  winD <- c(0, rev(seq(from = Time, to = 1, by = -block)))


  ## 0.3. Values returned
  dis <- NULL
  sig <- NULL
  w <- NULL
  thresh <- NULL
  qthresh <- NULL
  signal_recovered <- NULL

  ### End of 0. Data Preparation ###


  #### 1. Function for WaveL2E ####
  coefs = function(X) {
    crit = function(x,X=X,d=d) {
      sig = exp(x[1]);
      # The sigma parameter to be solved for when minimizing the L2E criterion
      w = exp(x[2]);
      # The weight parameter to be solved for when min the L2E criterion
      fX = exp(log(dnorm(X,0,sig)) %*% rep(1,d));
      # The Noise is N(0,sig)
      w^2/(2*sqrt(pi)*sig)^d - 2*w*mean( fX ) }
    # The L2E criterion

    if(is.matrix(X)) { d=ncol(X) } else { d=1; X=cbind(X,drop=F) }

    x0  = log(c(sd(X),0.8));
    # Starting value
    n=nrow(X);
    # The amount of rows to sort in the zeroing process
    # Added in constaints of sigma and w (w [0,1), sigma [0, inf))
    ans = nlminb( x0, crit, X=X, d=d, lower = c(-Inf, -Inf), upper = c(1000,-1e-5));
    # Eventhough transformations have been implemented, we still need constraints for large block sizes
    # When block size is too large revert to smaller interpretable block size or length of series (static implementation)
    # Optimization using PORT routines
    sig = exp(ans$par[1]);
    # Optimal Sigma
    w   = exp(ans$par[2]);
    # Optimal weight
    #print(c(sig,w))
    # Show sigma and weight

    # specify which coefficients should be zeroed

    dis0 = sqrt(X^2 %*% rep(1,d))
    # create positive values
    thresh = sort(dis0)[ w*n ]
    # identify the threshold value based on the optimal weight
    izero = seq(n)[dis0<thresh]
    # index of the values that will be zeroed based on the threshold

    #Assuming Chi-Square Distribution
    index = qchisq(0.95, df = (length(dis0)-1))  #by definition df
    # Identify the critical value for the Chi-Squared distribution
    qthresh = sort(dis0)[w*index]
    if (is.na(qthresh)) {
      qthresh = max(dis0)
    } else {
      qthresh = qthresh
    }
    # identify the new threshold value based on the weight
    qizero = seq(n)[dis0<qthresh]
    # index of the values that will be zeroed based on the Chi-Squared threshold

    return( list( sig=sig, w=w, izero=izero, thresh=thresh,
                  dis0=dis0, qthresh = qthresh, qizero = qizero) )
  }
  ### End of 1. Function for WaveL2E ###


  #### 2. Main function ####

  ## 2.1. Apply L2E
  ind_block <- 1:(length(winD)-1)

  signal_izero  <- matrix(data = FALSE, nrow = nrow_power, ncol = Time)
  signal_qizero <- matrix(data = FALSE, nrow = nrow_power, ncol = Time)

  for(i in ind_block){
    # X <- as.matrix(Power[,winD[i]:(winD[i+1]-1)])
    X <- as.matrix(Power[,(winD[i]+1):winD[i+1]])
    print(paste("loop", i, ", index:", winD[i]+1, "to", winD[i+1], ", length of", winD[i+1] - winD[i]))

    ans = coefs(X);

    #### Changes ###########################
    # Changed these to assign the same value to the whole block - so we have a comparable time series for each of these ###
    sig[(winD[i]+1):winD[i+1]] <- (ans$sig)
    w[(winD[i]+1):winD[i+1]] <- (ans$w)
    thresh[(winD[i]+1):winD[i+1]] <- ans$thresh
    qthresh[(winD[i]+1):winD[i+1]] <- ans$qthresh
    ######################################
    dis <- cbind(dis, ans$dis0)
    signal_izero[ans$izero,(winD[i]+1):winD[i+1]] <- TRUE
    signal_qizero[ans$qizero,(winD[i]+1):winD[i+1]] <- TRUE

  }

  ## 2.2. Threshold
  if (L2E) {
    my.w1$Power[signal_izero] <- 0
    my.w1$Wave[signal_izero] <- 0
  }

  if (Chi_square) {
    my.w2$Power[signal_qizero] <- 0
    my.w2$Wave[signal_qizero] <- 0
  }
  ### End of 2. Main function ###


  #### 3. Conditional Output ####

  ## 3.1. L2E = TRUE
  if (L2E) {
    # 1. Wavelet Plot
    CoFESWave.image(my.w1, periodlab = " ",timelab = "  " , main = " ",
                          legend.params = list(lab = "wavelet power levels", mar = 5.1, cex = 6, n.ticks = 10),
                          color.key = "quantile", lwd = 2, plot.ridge = FALSE, color.palette = "topl_colors(n.levels)",
                          show.date = Signal_date)
    title(latex2exp::TeX(" CWT Power Spectrum Thresholded - L_2E"), cex.main = 1.3,
          xlab = " ", ylab = " ", cex.lab = 1.3)

    # 2. Reconstructed Time-Seires
    recon_org1 <- CoFESWave.reconstruct(my.w1, plot.waves = FALSE, lwd = c(1,2),
                                           legend.coords = "topright",plot.rec = FALSE, verbose = F)

    # 3. Percentage of Total Volume
    PTV_1 <- 100*mean(as.matrix(my.w1$Power)/as.matrix(my.w_original$Power))

    # 4. Percentage of Significance Level
    PSL_1 <- 100*base::sum(as.matrix(my.w1$Power))/base::sum(as.matrix(my.w_original$Power))
  }

  ## 3.2. Chi_square = TRUE
  if (Chi_square) {
    # 1. Wavelet Plot
    CoFESWave.image(my.w2, periodlab = " ",timelab = "  " , main = " ",
                          legend.params = list(lab = "wavelet power levels", mar = 5.1, cex = 6, n.ticks = 10),
                          color.key = "quantile", lwd = 2, plot.ridge = FALSE, color.palette = "topl_colors(n.levels)",
                          show.date = Signal_date)
    title(latex2exp::TeX(" CWT Power Spectrum Thresholded - L_2E_{$\\chi$^2}"), cex.main = 1.3,
          xlab = " ", ylab = " ", cex.lab = 1.3)

    # 2. Reconstructed Time-Seires
    recon_org2 <- CoFESWave.reconstruct(my.w2, plot.waves = FALSE, lwd = c(1,2),
                                           legend.coords = "topright",plot.rec = FALSE, verbose = F)

    # 3. Percentage of Total Volume
    PTV_2 <- 100*mean(as.matrix(my.w2$Power)/as.matrix(my.w_original$Power))

    # 4. Percentage of Significance Level
    PSL_2 <- 100*base::sum(as.matrix(my.w2$Power))/base::sum(as.matrix(my.w_original$Power))
  }

  ########## Change ##############
  # Empirical WaveL2E + MAD
  # This is temporary we need to fix the code so we do evaluate block0 not disgard
  sig <- zerofill(sig)
  w <- zerofill(w)
  signal_recovered_mad <- (x - mad(w)*rnorm(n = 1, mean = 0, sd = mad(sig)))/(1-mad(w))
  signal_recovered_mad <- (x - mad(w,na.rm = TRUE)*rnorm(n = 1, mean = 0, sd = mad(sig,na.rm = TRUE)))/(1-mad(w,na.rm = TRUE))

  # The reason there results might not be as accurate could be a scaling issue
  # see https://github.com/cran/WaveletComp/blob/master/R/reconstruct.R - line 99
  signal_recovered <- (x - (w)*rnorm(n = 1, mean = 0, sd = (sig)))/(1-(w))
  signal_recovered <- (x - (w)*rnorm(n = 1, mean = 0, sd = (sig)))/(1-(w))
  #x.r <- (x.r-mean(x.r))*sd(x)/sd(x.r) + mean(x)
  signal_recovered <- (signal_recovered - mean(signal_recovered))*sd(x)/sd(signal_recovered) + mean(x)
  signal_recovered <- (signal_recovered - mean(signal_recovered,na.rm = TRUE))*sd(x,na.rm = TRUE)/sd(signal_recovered,na.rm = TRUE) + mean(x,na.rm = TRUE)

  ##################################

  ## 3.3. Plot of series: Data
  if (L2E & Chi_square) {
    df.plot1 <- recon_org1$series
    df.plot2 <- recon_org2$series
    df.plot3 <- signal_recovered

    df.plot <- cbind(df.plot1, df.plot2[,ncol(df.plot2)])

    if (Signal_date) {
      df.xts <- xts(df.plot[,2:4], order.by = df.plot$date)
    } else {
      # df.xts <- xts(df.plot, order.by = index(df.plot))
      df.xts <- ts(df.plot)
    }

    colnames(df.xts) <- c("original series", "reconstruction: L_2E", "reconstruction: L_2E Chi_square")
    title_series <- latex2exp::TeX(" CWT Series - L_2E, L_2E_{$\\chi$^2}")

  } else if (L2E & !Chi_square) {
    df.plot <- recon_org1$series

    if (Signal_date) {
      df.xts <- xts(df.plot[,2:3], order.by = df.plot$date)
    } else {
      df.xts <- xts(df.plot)
    }
    colnames(df.xts) <- c("original series", "reconstruction: L_2E")
    title_series <- latex2exp::TeX(" CWT Series - L_2E")

  } else if (!L2E & Chi_square) {
    df.plot <- recon_org2$series

    if (Signal_date) {
      df.xts <- xts(df.plot[,2:4], order.by = df.plot$date)
    } else {
      df.xts <- xts(df.plot)
    }
    colnames(df.xts) <- c("original series", "reconstruction: L_2E Chi_square", "EWaveL_2E")
    title_series <- latex2exp::TeX(" CWT Series - L_2E_{$\\chi$^2} and EWaveL_2E")

  } else {
    df.plot <- my.w_original$series

    if (Signal_date) {
      df.xts <- xts(df.plot[,2], order.by = df.plot$date)
    } else {
      df.xts <- xts(df.plot)
    }
    colnames(df.xts) <- c("original series")
    title_series <- latex2exp::TeX(" CWT Series")

  }

  ## 3.4. Plot of series
  if (Signal_date) {
    plot.series <- plot(df.xts,
         # col = topl_colors(2),
         col = c("#64d8cb","#26a69a", "#90a4ae", "#5f5fc4"),
         lwd = c(1,1.5,2.5), lty = c(1,6,3),
         grid.col = NA, yaxis.right = FALSE,
         legend.loc = "topleft", auto.legend=TRUE,
         main=title_series)
    # addLegend(legend.loc = "topleft", col = c("#64d8cb","#26a69a", "#90a4ae"),
    #           legend.names = c("test1", "test2", "test3"))
    matplot(sig, type = "l", main = "Noise SD")
    matplot(w, type = "l", main = "Threshold Weight")
    print(plot.series)
  } else {
    matplot(df.xts, type = "l",
            col = c("#64d8cb","#26a69a", "#90a4ae", "#5f5fc4"),
            lwd = c(1,1.5,2.5), lty = c(1,6,3),
            main=title_series,
            ylab = '')
    legend("topleft", fill = c("#64d8cb","#26a69a", "#90a4ae", "#5f5fc4"),
           legend = colnames(df.xts),
           box.lty=0, cex=0.7, pt.cex = 1)
    matplot(sig, type = "l", main = "Noise SD")
    matplot(w, type = "l", main = "Threshold Weight")
  }
  ### End of 3. Conditional Output ###



  #### 4. Output ####
  output<-list(sig=sig, w=w, dis0=dis,
               thresh=thresh, qthresh = qthresh,
               original = my.w_original,
               Ana_Wave = my.w1,
               Emp_WaveL2E = signal_recovered,
               Emp_WaveL2E_MAD = signal_recovered_mad)

  if (L2E) {
    output <- c(output,
                list(recon_L2E = recon_org1,
                     PTV_L2E = PTV_1,
                     PSL_L2E = PSL_1))
  }

  if (Chi_square) {
    output <- c(output,
                list(recon_Chi_square = recon_org2,
                     PTV_Chi_square = PTV_2,
                     PSL_Chi_square = PSL_2))
  }

  if (!is.null(date)) {
    output <- c(output,
                list(date = date))
  }

  class(output) <- "WaveL2E"
  return(invisible(output))
  ### End of 4. Output ###


  ### End of the function ###
}
