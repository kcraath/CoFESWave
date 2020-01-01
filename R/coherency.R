#' Title
#'
#' @param x
#' @param y
#' @param dt
#' @param dj
#' @param low.period
#' @param up.period
#' @param pad
#' @param sigma
#' @param wt.type
#' @param wt.size
#' @param ws.type
#' @param ws.size
#' @param n.sur
#' @param p
#' @param q
#' @param Phase_diff
#' @param low.fp
#' @param up.fp
#' @param date
#'
#' @return
#' @export
#'
#' @examples
CoFEScoherency <- function(x,y,
                           dt=1,dj=0.25,
                           low.period=2*dt,up.period=length(x)*dt,
                           pad=0,sigma=1,
                           wt.type=0,wt.size=5,
                           ws.type=0,ws.size=5,
                           n.sur=0,
                           p=0,q=0,
                           Phase_diff = FALSE,
                           low.fp = 32,up.fp = 128,
                           date = NULL)
{
  # make test and debug easier.
  # x = rCGW[,1];y = rXLE[,1];
  # dt=1;dj=0.25;
  # low.period=1;up.period=512;
  # pad=0;sigma=1;
  # wt.type=0;wt.size=5;
  # ws.type=0;ws.size=5;
  # n.sur=0;
  # p=0;q=0;
  # Phase_diff = TRUE;
  # low.fp = 32;up.fp = 128;

  #### 0. Data Preparation ####

  ## 0.1. Test input
  if (missing(x)||missing(y)) {stop("You must input two vectors")}


  if (!(is.numeric(x) && is.vector(x)))
  {
    stop(sprintf("argument %s must be a vector", sQuote("x")))
  }

  if (!(is.numeric(y) && is.vector(y)))
  {
    stop(sprintf("argument %s must be a vector", sQuote("y")))
  }

  nTimes <- length(x)
  if(nTimes!=length(y)){stop("vectores must be of same length")}


  ## 0.2. Computation of quantities that depend on the wavelet
  if (sigma==1){
    center.frequency=6.0    #  No need to compute: This corresponds to Morlet with omega_0=6.
  } else {
    tol<- 5e-8   #  May be changed, if we wish
    ck<- tol/(sqrt(2)*pi^(0.25))
    center.frequency<- sqrt(2)*(sqrt(log(sqrt(sigma))-log(ck)))/sigma
    center.frequency=ceiling(center.frequency)   #  choose an integer center.frequency
  }
  #  Computation of Fourier factor
  fourier.factor <- (2*pi)/center.frequency

  #  Computation of radius in time
  sigma.time <- sigma/sqrt(2)

  ## 0.3. Padding
  #  Computation of extra.length (number of zeros used to pad x with; it depends on pad)
  if ( (pad > 0) && (log2(pad)%%1 == 0) ) {
    #  Zero padding to selected size
    pot2 <- log2(pad)
    extra.length <- 2^pot2-nTimes
    if (extra.length <= 0){
      print("pad smaller than size of series; next power of 2 used")
      #  Zero padding to size=next power of 2+2
      pot2 <- ceiling(log2(nTimes))
      extra.length <- 2^(pot2+2)-nTimes
    }
  } else {
    #  Zero padding to size=next power of 2+2
    pot2 <- ceiling(log2(nTimes))
    extra.length <- 2^(pot2+2)-nTimes
  }

  ## 0.4. Computation of scales and periods
  s0 <- low.period/fourier.factor   #  Convert low.period to minimum scale

  if (up.period > nTimes*dt)
  { print("up.period is too long; it will be adapted")
    up.period<- nTimes*dt
  }
  up.scale <- up.period/fourier.factor       #  Convert up.period to maximum scale
  J <- as.integer(log2(up.scale/s0)/dj)    #  Index of maximum scale
  scales <- s0*2^((0:J)*dj)   #  Vector of scales
  nScales <- length(scales)   #  Number of scales
  periods <- fourier.factor*scales   #  Conversion of scales to periods

  ## 0.5. Computation of coi
  coiS <- fourier.factor/sigma.time   #
  coi <-  coiS*dt*c(1e-8,1:floor((nTimes-1)/2),floor((nTimes-2)/2):1,1e-8)

  ## 0.6. Computation of angular frequencies
  N <- nTimes+extra.length
  wk <- (1:floor(N/2))
  wk <-  wk*((2*pi)/(N*dt))
  wk <-  c(0., wk, -wk[(floor((N-1)/2)):1])

  ### End of 0. Data Preparation ###



  #### 1. Functions ####

  ## 1.1. Gabor wavelet tranformation
  gabor.wtransform <- function(x)
  {
    x <- (x-mean(x))/sd(x)
    xn <- c(x,rep(0,extra.length))   #  Pad x with zeros

    ftxn <- fft(xn) #  Computation of FFT of xn
    #  Computation of Wavelet Transform of x   #
    wave <- matrix(0,nScales,N)   #  Matrix to accomodate WT
    wave <- wave + 1i*wave;       #  Make it complex
    for (iScales in (1:nScales)) {
      #  Do the computation for each scale
      #  (at all times simultaneously)
      scaleP <- scales[iScales]   #  Particular scale
      norm <- pi^(1/4)*sqrt(2*sigma*scaleP/dt)
      expnt <- -( ((scaleP*wk-center.frequency)*sigma)^2/2) *(wk > 0)
      daughter <- norm*exp(expnt)
      daughter<- daughter*(wk>0)
      wave[iScales,] <- (fft(ftxn*daughter,inv=TRUE))/length(xn)
    }
    #  Truncate WT to correct size
    wave<- wave[,1:nTimes]
    return(wave)
  }


  ## 1.2. Window function
  ## Windows used in smoothing: necessary for coherency
  window  <-  function(wtype=0,n=5) {
    if  (!(n  ==  round(n))||  n<=0){
      stop("Size  of  window  must  be  a  positive  integer")
    }


    if  (wtype==0||  wtype=='bar'){
      w <- 2*(0:((n-1)/2))/(n-1)

      if  (n%%2){
        w <- c(w,w[((n-1)/2):1])
      } else {
        w <- c(w,w[(n/2):1]) }
    }

    # Case  wtype=1  --  Blackmann  window
    if  (wtype==1||  wtype=='bla')
    {w  <-  .42  -  .5*cos(2*pi*(0:(n-1))/(n-1) )  +  .08*cos(4*pi*(0:(n-1))/(n-1))
    }

    # Case  wtype=2  --  Rectangular  window  (boxcar)
    if  (wtype==2||  wtype=='rec')
    {  w  <-  rep(1,n)
    }

    # Case  wtype=3  --  Hamming  window
    if  (wtype==3||wtype=='ham')
    {
      if  (n  >  1) {
        w  <-  .54  -  .46*cos( 2*pi*( 0:(n-1) ) / (n-1) )
      } else {
        w  <-  .08
      }
    }

    # Case  wtype=4  --  Hanning  window
    if  (wtype==4||wtype=='han') {
      w  <-  .5*(1  -  cos(( 2*pi*(1:n))/(n+1)))
    }

    # Case  wtype=5  --  Triangular  window
    if  (wtype==5||wtype=='tri') {
      if  (n%%2) {
        w  <-  2*(1:((n+1)/2))/(n+1)
        w  <-  c(w,  w[((n-1)/2):1])
      }  else  {
        w  <- (2*(1:((n+1)/2))-1)/n
        w  <-  c(w,  w[(n/2):1])
      }
    }
    w <- w/sum(w)    # Normalize
    return(w)
  }


  ## 1.3. coherency
  ## computes complex wavelet coherency)
  gw.coherency <- function(x,y)
  {
    # Computation of Wavelet Transforms of x and y
    # Uses function gabor.wtransform
    wt.x <- gabor.wtransform(x)
    wt.y <- gabor.wtransform(y)
    # Computation of Wxy (cross), Wxx and Wyy
    Wxy <- wt.x*Conj(wt.y)
    Wxx <- (abs(wt.x))^2
    Wyy <- (abs(wt.y))^2
    sWxy<- Wxy  # To acomodate sWxy (smoothed cross)

    # Window for smoothing in scale
    ws.size <-as.integer(ws.size/(2*dj))
    FILS <- window(ws.type,ws.size)
    lFILS<- length(FILS)
    # Smoothing in scale direction
    for (iTime in (1:nTimes)) {
      aux <-  convolve(Wxx[,iTime],rev(FILS),type='open')
      Wxx[,iTime]<- aux[floor((lFILS+2)/2):( floor((lFILS+2)/2)+ nScales-1)]

      aux <-  convolve(as.vector(Wyy[,iTime]),rev(FILS),type='open')
      Wyy[,iTime]<- aux[floor((lFILS+2)/2):(floor((lFILS+2)/2)+nScales-1)]

      aux <-  convolve(as.vector(sWxy[,iTime]),rev(FILS),type='open')
      sWxy[,iTime]<- aux[floor((lFILS+2)/2):(floor((lFILS+2)/2)+nScales-1)]
    }
    #  Smoothing in time direction
    for (iScale in (1:nScales)) {
      wTSS <-as.integer(periods[iScale]*wt.size/dt)
      # Size of window is adapted to scale
      # Window for smoothing in time
      FILT <- window(wt.type,wTSS)
      lFILT<-length(FILT)

      aux <- convolve(as.vector(Wxx[iScale,]),rev(FILT),type='open')
      Wxx[iScale,] <- aux[floor((lFILT+2)/2):(floor((lFILT+2)/2)-1+nTimes)]

      aux <- convolve(as.vector(Wyy[iScale,]),rev(FILT),type='open')
      Wyy[iScale,] <- aux[floor((lFILT+2)/2):(floor((lFILT+2)/2)-1+nTimes)]

      aux <- convolve(as.vector(sWxy[iScale,]),rev(FILT),type='open')
      sWxy[iScale,]<- aux[floor((lFILT+2)/2):(floor((lFILT+2)/2)-1+nTimes)]
    }

    # Computation of  coherency
    wco <- sWxy/(sqrt(Wxx*Wyy))
    wco[which(Wxx<1e-10)]<- 0
    wco[which(Wyy<1e-10)]<- 0  # To avoid meaningless numbers

    # when Wxx or Wyy are nearly zero
    # Coherency in absolute value
    rwco<-abs(wco)
    output<-list(wco=wco,rwco=rwco,cross=Wxy,scross=sWxy)
    return(output)
  }
  ### End of 1. Functions ###



  #### 2. Main function ####

  ## 2.1. Main
  WCO <- gw.coherency(x,y)

  ## 2.2. details
  wco<-WCO$wco
  rwco<-WCO$rwco
  cross<-WCO$cross
  scross<-WCO$scross

  ### End of 2. Main function ###



  #### 3. Conditional Output ####

  ## 3.1. n.sur > 0
  # Computation of pvCO

  if (n.sur>0) {
    nRows<- nrow(wco)
    nCols<- ncol(wco)
    pvCo<- matrix(0,nRows,nCols)
    for  (iSur in 1:n.sur) {
      # Surrogate x and y
      if(p==0&q==0) {
        xSur <- rnorm(nTimes)
        ySur <- rnorm(nTimes)
      } else {
        fitx<-stats::arima(x,order=c(p,0,q))
        fity<-stats::arima(y,order=c(p,0,q))
        coefsx<-coef(fitx)
        coefsy<-coef(fity)

        if (p==0) {
          ar.coefsx<-c(0)
          ma.coefsx<-coefsx[1:q]
          ar.coefsy<-c(0)
          ma.coefsy<-coefsy[1:q]
        } else {
          ar.coefsx <- coefsx[1:p]
          ar.coefsy <- coefsy[1:p]
          if (q==0) {
            ma.coefsx<-c(0)
            ma.coefsy<-c(0)
          } else {
            ma.coefsx<-coefsx[(p+1):(p+q)]
            ma.coefsy<-coefsy[(p+1):(p+q)]
          }
        }
        modelx <- list(ar=ar.coefsx,ma=ma.coefsx)
        modely <- list(ar=ar.coefsy,ma=ma.coefsy)
        xSur <- stats::arima.sim(modelx,nTimes)
        ySur <- stats::arima.sim(modely,nTimes)
      }


      # Coherency of surrogates
      wcoaux <- gw.coherency(xSur,ySur)
      rwcoSur<-wcoaux$rwco
      pvCo[which(rwcoSur>=rwco)] <- pvCo[which(rwcoSur>=rwco)]+1
    }
    pvCo <- pvCo/n.sur
  }


  ## 3.2. Phrase_diff = TRUE
  if (Phase_diff) {
    cross.power<-abs(cross)
    average.cross<- as.vector(t(rowMeans(cross.power))) # Average Cross-Wavelet Power

    rwco<-WCO$rwco      # real coherency
    # periods<-WCO$periods

    if (missing(low.fp)) { low.fp <- periods[1] }
    if (missing(up.fp))  { up.fp  <- periods[length(periods)]}

    if (low.fp<periods[1]) {
      print("low.fp too low; we use first period")
      low.fp <- periods[1]
    }
    if (up.fp>periods[length(periods)]) {
      print("up.fp too high; we use last period")
      up.fp <- periods[length(periods)]
    }

    selPeriods<-((periods >= low.fp) & (periods <= up.fp))
    # Indexes of selected periods

    average.coer <- as.vector(t(rowMeans(rwco)))  # Average Wavelet Coherency
    # (i.e. time-average over all times)
    # Computation of phase.dif (output)
    phase.dif <- Arg(scross)  # Phase-Difference (use cross if you
    # prefer the other formula for phase difference
    # (see comment p.16, NIPE WP 16/2011)

    phase.dif <-colMeans(phase.dif[selPeriods,])
    # Mean Phase-Difference for selected periods

    # Computation of the instantaneous time-lag
    mean.period<- mean(periods[selPeriods])  # Mean of selected periods

    time.lag <- (phase.dif*mean.period)/(2*pi)  # Instantaneous Time-Lag
    # between  x and y for
    # selected periods
  }
  ### End of 3. Conditional Output ###



  #### 4. Output ####
  output <- list(wco=wco,
                 rwco=rwco,
                 cross=cross,
                 scross=scross,
                 periods=periods,
                 scales=scales,
                 coi=coi)

  if (n.sur>0){
    output<-c(output,
              list(pv=pvCo) )
  }

  if (Phase_diff){
    output<-c(output,
              list(phase.dif=phase.dif,
                   time.lag=time.lag,
                   average.coer=average.coer,
                   average.cross=average.cross) )
  }

  if (!is.null(date)) {
    output <- c(output,
                list(date = date))
  }

  class(output) <- "CoFEScoherency"
  return(invisible(output))
  ### End of 4. Output ###


  ### End of the function ###
}
