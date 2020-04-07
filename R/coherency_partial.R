#' CoFESmpcoherency
#'
#' @title Multiple or Partial Wavelet Coherencies of several series
#' @aliases CoFESmpcoherency
#' @description This function computes multiple(wmco) and/or partial wavelet coherencies (wpco)
#' of columns of  matrix X. These coherencies are computed using a Gabor Wavelet
#' Transform (different Gabor wavelets can be used).
#' The user has the options to use different windows for smoothing.
#' The user can also compute the corresponding p_values
#' (if we take n.sur>0).
#' The surrogates are constructed by fitting an ARMA(p,q)
#' model to our series and building new samples by drawing
#' errors from a  Gaussian distribution.
#'
#' The major part of the code has been adopted from L. Aguiar-Conraria and M.J. Soares.
#' @param X matrix with m columns (vectors corresponding to time-series)
#' whose partial or multiple coherencies we want to compute.
#' The first column X[,1] has a special role;
#' e.g. in case of multiple coherency, we want to compte R1.(2...m)
#' @param dt sampling rate
#'
#' Default: \code{1}.
#' @param dj frequency resolution
#'  (i.e. 1/dj = number of voices per octave)
#'
#' Default: \code{0.25}.
#' @param low.period lower period of the decomposition
#'
#' Default: \code{2*dt}.
#' @param up.period upper period of the decomposition
#'
#' Default: \code{length(x)*dt}.
#' @param pad an integer (power of 2) defining the total length of
#'            the vectors x and x after zero padding;
#'            if pad is not a power of 2, pad with zeros to
#'            total length: 2^(next power of 2 + 2)
#' @param sigma the sigma parameter for the Gabor wavelet
#'
#' Default: \code{1 - Morlet Wavelet}.
#' @param coher.type type of coherency we want to compute
#' \tabular{rlll}{
#'   \tab \code{0} or "part"  \tab : \tab partial (default) \cr
#'   \tab \code{1} or "mult"  \tab : \tab multiple
#' }
#'
#' @param index.p index of series for partial coherency
#'
#' Default: \code{2}.
#' @param wt.type type of window for smoothing in time direction
#' \tabular{rlll}{
#'   \tab \code{0}  \tab : \tab Bartlett \cr
#'   \tab \code{1}  \tab : \tab Blackman \cr
#'   \tab \code{2}  \tab : \tab Rectangular \cr
#'   \tab \code{3}  \tab : \tab Hamming \cr
#'   \tab \code{4}  \tab : \tab Hanning \cr
#'   \tab \code{5}  \tab : \tab Triangular
#' }
#'
#' Default: \code{Bartlett}.
#' @param wt.size (half) size of window in time
#'
#' Default: \code{5}.
#' @param ws.type type of window for smoothing in scale direction
#' \tabular{rlll}{
#'   \tab \code{0}  \tab : \tab Bartlett \cr
#'   \tab \code{1}  \tab : \tab Blackman \cr
#'   \tab \code{2}  \tab : \tab Rectangular \cr
#'   \tab \code{3}  \tab : \tab Hamming \cr
#'   \tab \code{4}  \tab : \tab Hanning \cr
#'   \tab \code{5}  \tab : \tab Triangular
#' }
#'
#' Default: \code{Bartlett}.
#' @param ws.size (half) size of window in scale
#'
#' Default: \code{5}.
#' @param n.sur integer, number of surrogate series, if we want
#'             to compute p-values for the Wavelet Coherency
#'
#' Default: \code{0 - no computation}.
#' @param p non-negative integers, orders of the ARMA(P,Q) model
#'
#' Default: \code{0}.
#' @param q non-negative integers, orders of the ARMA(P,Q) model
#'
#' Default: \code{0}.
#' @param Phase_diff Calculate the phase-difference? Logical.
#'
#' Default: \code{TRUE}.
#' @param low.fp lower periods used in the computation of phase-difference
#'
#' Default: \code{32}.
#' @param up.fp upper periods used in the computation of phase-difference
#'
#' Default: \code{128}.
#' @param date a date series
#'
#' @return A list of class \code{"CoFESmpcoherency"} with elements of different dimensions.
#' %%%%%%%%%%%%%%%%%
#' Here is a detailed list of all elements:
#'#' %%%%%%%%%%%%%%%%%
#' \item{wmpco}{Wavelet Multiple or Partial Coherency Matrix}
#' \item{periods}{the vector of Fourier periods (in time units)
#'                that corresponds to the used scales}
#' \item{scales}{the vector of scales, given by
#'                  s0*2^(j*dj);j=0,...,J1,
#'            where J1+1 is  number of scales and s0 is
#'            minimum scale}
#' \item{coi}{the "cone-of-influence", which is a vector of
#'            n_x=length(x) points that contains the limit of the
#'            region where the wavelet transforms are influenced
#'            by edge effects}
#' \item{low.fp}{lower periods where phase-diff is computed}
#' \item{up.fp}{upper periods where phase-diff is computed}
#' \item{pvMP}{matrix of p-values for the multiple/partial wavelet coherency}
#' \item{phase.dif}{mean phase difference (in the selected periods)}
#' \item{time.lag}{time-lag (in the selected periods)}
#' \item{average.wpco}{mean partial coherency (we take absolute value of partial coherency)}
#' \item{date}{a date series}
#' %%%%%%%%%%%%%%%%%
#'
#' @author CoFES. Credits are also due to L. Aguiar-Conraria and M.J. Soares.
#'
#' @references
#' Aguiar-Conraria, L. and Soares, M.J. (2010)
#'
#' Soares, M.J. (2010),
#' "Multiple and Partial Wavelet Coherencies"
#'  (private notes) available at http:
#'
#' Torrence, C. and Compo, T.C.,
#' "A Prectical Guide to Wavelet Analysis" (1998),
#' Bulletin of the American Meteorological Society, 79, 605.618.
#'
#' @examples \dontrun{
#' ## The following example is adopted from Raath et al., 2020:
#'
#' ## Please first run the examples in the WaveL2E and CoFEScoherency functions.
#'
#' # ----- Computation of Partial Coherency of Water , Energy (controlling for SPY) ---- #
#'
#' index.p=2
#' Data_thresh <- cbind(periodic_waveL2E_2$recon_L2E$series$x.r,
#'                      periodic_waveL2E$recon_L2E$series$x.r,
#'                      periodic_waveL2E_3$recon_L2E$series$x.r)
#' X<-matrix(Data_thresh,T,3)
#' W<-CoFESmpcoherency(X, coher.type='part',index.p=index.p,
#'                     low.period=low.period,up.period=up.period,
#'                     low.fp = lowFP1,up.fp = upFP1, Phase_diff = TRUE, date = date1)
#'
#' plot_CoFESWaveWCO(W,pE=5)
#' }
#'
CoFESmpcoherency <- function(X,
                             dt=1,dj=0.25,
                             low.period=2*dt,up.period=nrow(X)*dt,
                             pad=0,sigma=1,coher.type=0,index.p=2,
                             wt.type=0,wt.size=5,ws.type=0,ws.size=5,
                             n.sur=0,
                             p=0,q=0,
                             Phase_diff = FALSE,
                             low.fp = 32,up.fp = 128,
                             date = NULL)
{
  # make test and debug easier.
  # X = as.matrix(Data);
  # dt=1;dj=0.25;
  # low.period=2*dt;up.period=nrow(X)*dt;
  # pad=0;sigma=1;coher.type=0;index.p=2;
  # wt.type=0;wt.size=5;ws.type=0;ws.size=5;
  # n.sur=0;
  # p=0;q=0;
  # Phase_diff = TRUE;
  # low.fp = 32;up.fp = 128

  #### 0. Data Preparation ####

  ## 0.1. Test input
  if (missing(X)) { stop("Must input matrix X") }
  if (!(is.numeric(X) && is.matrix(X))) { stop(sprintf("argument %s must be a matrix", sQuote("X"))) }

  nTimes<-nrow(X)   #  Time length


  ## 0.2. Computation of quantities that depend on the wavelet
  if (sigma==1) {
    center.frequency=6.0    #  No need to compute: This corresponds to Morlet with omega_0=6.
  } else {
    tol<- 5e-8   #  May be changed, if we wish
    ck<- tol/(sqrt(2)*pi^(0.25))
    center.frequency<- sqrt(2)*(sqrt(log(sqrt(sigma))-log(ck)))/sigma
    center.frequency=ceiling(center.frequency)   #  choose an integer center.frequency
  }

  #  Computation of Fourier factor
  fourier.factor<- (2*pi)/center.frequency

  #  Computation of radius in time
  sigma.time<-sigma/sqrt(2)

  ## 0.3. Padding
  #  Computation of extra.length (number of zeros used to pad x with; it depends on pad)
  if ( (pad > 0) && (log2(pad)%%1 ==0) ) {
    # Zero padding to selected size
    pot2 <- log2(pad)
    extra.length <- 2^pot2-nTimes
    if (extra.length <=0) {
      print("pad smaller than size of series; next power of 2 used")
      # Zero padding to size=next power of 2+2
      pot2 <- ceiling(log2(nTimes))
      extra.length <- 2^(pot2+2)-nTimes
    }
  } else {
    # Zero padding to size=next power of 2+2
    pot2 <- ceiling(log2(nTimes))
    extra.length <- 2^(pot2+2)-nTimes
  }

  ## 0.4. Computation of scales and periods
  s0 <- low.period/fourier.factor   #  Convert low.period to minimum scale

  if (up.period>nTimes*dt) {
    print("up.period is too long; it will be adapted")
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
  wk<-  c(0., wk, -wk[(floor((N-1)/2)):1])

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


  ## 1.3. smooth
  ## to smooth the spectra
  smooth <- function(X) {
    #  Window for smoothing in scale
    ws.size <-as.integer(ws.size/(2*dj))
    FILS <- window(ws.type,ws.size)
    lFILS<- length(FILS)

    # Smoothing in scale direction
    for (iTime in 1:nTimes) {
      aux <-  convolve(X[,iTime],rev(FILS),type='open')
      X[,iTime]<- aux[floor((lFILS+2)/2):(floor((lFILS+2)/2)-1+nScales)]
    }

    # Smoothing in time direction
    for (iScale in 1:nScales) {
      wTSS <-as.integer(periods[iScale]*wt.size/dt)
      #  Size of window is adapted to scale
      #  Window for smoothing in time
      FILT <- window(wt.type,wTSS)
      lFILT<-length(FILT)

      aux <- convolve(X[iScale,],rev(FILT),type='open')
      X[iScale,]<- aux[floor((lFILT+2)/2):(floor((lFILT+2)/2)-1+nTimes)]
    }

    return(X)
  }

  detC<-function(A){
    deter<- prod(eigen(A, only.values=TRUE)$values)
    return(deter)
  }

  ## 1.4. mp coherency
  ## compute multiple and partial coherencies
  mpw.coherency <- function(X) {
    nSeries <- ncol(X)      #  Number of series
    S <- vector('list',nSeries*nSeries)
    dim(S) <- c(nSeries,nSeries)   #  Array  to accomodate all
    #  the smoothed spectra matrices
    #  e.g  S[[1,2]]=S(Wx1x2)
    SWT<-vector('list',nSeries)

    #  Computation of all smoothed cross spectra (will be saved in SWT)
    for (ii in 1:nSeries) {
      SWT[[ii]]<-gabor.wtransform(X[,ii])  #  Save WTs of all series
    }

    for (iiRow in 1:(nSeries)) {
      for (iiCol in (iiRow):nSeries) {
        SXY <- SWT[[iiRow]]*Conj(SWT[[iiCol]])
        SXY<-smooth(SXY)
        S[[iiRow,iiCol]] <- SXY

        if (iiRow!=iiCol) {
          S[[iiCol,iiRow]]<-Conj(SXY)
        }
      }
    }
    MPWCO <- matrix(0, nScales,nTimes)   #  Initialize MPWCO
    S11<-S[2:nSeries,2:nSeries]

    # PARTIAL COHERENCY
    # (rho_1j.(2...j-1 j+1  ...m))
    if (coher.type=='part' || coher.type==0) {
      MatNum <- matrix(0,nSeries-1,nSeries-1)
      MatDen1 <- matrix(0,nSeries-1,nSeries-1)
      MatDen2 <- matrix(0,nSeries-1,nSeries-1)

      vrows <- 1:nSeries
      J=index.p;   #  This is just for simplicity
      vrows<-vrows[-J]
      SJ1 <- S[vrows,2:nSeries]

      SJJ <- S[vrows,vrows]

      for (nscale in 1:nScales)
      {
        for (ntime in 1:nTimes)
        {
          for (iiRow in 1:(nSeries-1))
          {
            for (iiCol in 1:(nSeries-1))
            {
              MatNum[iiRow,iiCol]<- SJ1[[iiRow,iiCol]][nscale,ntime]

              MatDen1[iiRow,iiCol]<- S11[[iiRow,iiCol]][nscale,ntime]

              MatDen2[iiRow,iiCol]<- SJJ[[iiRow,iiCol]][nscale,ntime]
            }
          }
          Num <-detC(MatNum)
          Den<- abs(detC(MatDen1)*detC(MatDen2))
          MPWCO[nscale,ntime] <- ((-1)^J)*(Num/sqrt(Den))
        }
      }
    }

    # MULTIPLE COHERENCY
    # R1^2.(2...m)
    if(coher.type=='mult' || coher.type == 1) {
      MatNum <- matrix(0,nSeries,nSeries)
      MatDen <-matrix(0,nSeries-1,nSeries-1)

      for (nscale in 1:nScales)
      {
        for (ntime in 1:nTimes)
        {
          for (iiRow in 1:(nSeries))
          {
            for (iiCol in 1:(nSeries))
            {
              MatNum[iiRow,iiCol]<-S[[iiRow,iiCol]][nscale,ntime]
            }
          }
          for (iiRow in 1:(nSeries-1))
          {
            for (iiCol in 1:(nSeries-1))
            {
              MatDen[iiRow,iiCol]<- S11[[iiRow,iiCol]][nscale,ntime]
            }
          }
          MPWCO[nscale,ntime]=sqrt(abs(1-abs(detC(MatNum)/(S[[1,1]][nscale,ntime]*detC(MatDen)))))

        }
      }
    }

    return(MPWCO)
  }
  ### End of 1. Functions ###



  #### 2. Main function ####

  ## 2.1. Main
  MPWCO <- mpw.coherency(X)
  output <- list(mpwco=MPWCO,scales=scales,periods=periods,coi=coi)

  ### End of 2. Main function ###



  #### 3. Conditional Output ####

  ## 3.1. n.sur > 0
  # Computation of pvCO

  if (n.sur>0) {
    nSeries <- ncol(X)
    pvMP <- matrix(0,nScales,nTimes)   #  Initialize matrix pvMP
    XSur <- matrix(0,nrow(X),ncol(X))

    for  (iSur in 1:n.sur){
      for (jcol in 1:nSeries) {
        if(p==0&q==0) {
          XSur[,jcol] <-rnorm(nTimes)
        } else {
          fit <- stats::arima(X[,jcol],order=c(p,0,q))
          coefs <- coef(fit)
          if (p==0) {
            ar.coefs<-c(0)
            ma.coefs<-coefs[1:q]
          } else {
            ar.coefs<-coefs[1:p]
            if (q==0) {
              ma.coefs<-c(0)
            } else {
              ma.coefs<-coefs[(p+1):(p+q)]}
            }
          model<-list(ar=ar.coefs,ma=ma.coefs)
          surrogate <-stats::arima.sim(model,nTimes)
          XSur[,jcol]<- surrogate
          }
      }
      MPWCOSur<-mpw.coherency(XSur)
      pvMP[which(abs(MPWCOSur)>=abs(MPWCO))] <- pvMP[which(abs(MPWCOSur)>=abs(MPWCO))]+1
    }
    pvMP<-pvMP/n.sur
  }


  ## 3.2. Phrase_diff = TRUE
  if (Phase_diff) {
    wpco <- MPWCO
    rwpco <- abs(wpco) # real coherency

    # periods<-MPWCO$periods

    if(missing(low.fp)) { low.fp<- periods[1] }
    if (missing(up.fp)) { up.fp=periods[length(periods)] }

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

    if (n.sur>0) {
      min.prob <- apply(pvMP[selPeriods,],2,min)
      sig1 <- as.numeric(min.prob<0.05)
      isig1 <- which(sig1==1)
    }
    # code added

    average.wpco <- as.vector(t(rowMeans(rwpco)))  # Average Wavelet Coherency
    # (i.e. time-average over all times)
    # Computation of phase.dif (output)
    phase.dif <- Arg(wpco)  # Phase-Difference
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
  output<-list(mpwco=MPWCO,
               scales=scales,
               periods=periods,
               coi=coi,
               low.fp = low.fp,
               up.fp = up.fp)

  if (n.sur>0) {
    output <- c(output,
                list(pvMP=pvMP))
  }

  if (Phase_diff) {
    output <- c(output,
                list(phase.dif=phase.dif,
                     time.lag=time.lag,
                     average.wpco=average.wpco,
                     isig1 = isig1))
  }

  if (!is.null(date)) {
    output <- c(output,
                list(date = date))
  }

  class(output) <- "CoFESmpcoherency"
  return(invisible(output))
  ### End of 4. Output ###


  ### End of the function ###
}
