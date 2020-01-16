############################
##
## Methods Paper Figures
##
############################

# All the libraries
## You should make your repo public so other people can download this package.
## Or open this package and click on the "Build" button near the top-right corner of RStudio,
## then click on the "Install and Restart" button.
remove.packages("WaveCoFES")
remove.packages("CoFESWave")
devtools::install_github("kcraath/CoFESWave", auth_token = "ae2a03370845f12f616bed5b214903c3e56e68ae")
library(CoFESWave)


#install.packages("Metrics")
library("Metrics")
#install.packages('wavethresh')
library("wavethresh")
#install.packages('wmtsa')
library('wmtsa')
#Details (math included): https://www.rdocumentation.org/packages/wmtsa/versions/2.0-3/topics/wavShrink

#ebayesthresh
#https://cran.r-project.org/web/packages/EbayesThresh/EbayesThresh.pdf
#http://www.stats.ox.ac.uk/~silverma/ebayes/ebayesw.pdf
#install.packages("EbayesThresh")
library(EbayesThresh)
#install.packages("waveslim")
library(waveslim)

# An R Package of time series tools and utilities; Rmetrics - Financial Time Series Objects
#https://www.rdocumentation.org/packages/timeSeries
#install.packages("timeSeries")
require(timeSeries)

# An R package with a collection of econometric functions for performance and risk analysis
#https://www.rdocumentation.org/packages/PerformanceAnalytics
#install.packages("PerformanceAnalytics")
require(PerformanceAnalytics)

# R package which includes Quantitative Financial Modelling Frameworks.
#https://www.rdocumentation.org/packages/quantmod
#install.packages("quantmod")
require(quantmod)

# An R package for Wavelet analysis and reconstruction of time series,
# cross-wavelets and phase-difference (with filtering options),
# significance with simulation algorithms.
# https://www.rdocumentation.org/packages/WaveletComp/versions/1.0
#install.packages("WaveletComp")
require(WaveletComp)

# devtools: Tools to Make Developing R Packages Easier
#https://www.rdocumentation.org/packages/devtools
#install.packages("devtools")
require(devtools) # using devtools to download from github

# R package for time series analysis using the Wavelet Scalogram
# from https://github.com/rbensua/wavScalogram
#install_github("rbensua/wavScalogram")
require(wavScalogram)

# biwavelet: Conduct Univariate and Bivariate Wavelet Analyses
# https://www.rdocumentation.org/packages/biwavelet
#install.packages("biwavelet")
require(biwavelet)

#Latex expersions
# library(latex2exp)

#Latex tables
library(knitr)

#Interactive Charts
library(highcharter)

##################

###########################

# Load the WaveL2E Function

###########################

# Figure 1 - Set up #

# Identify the tickers of interest
tickers <- c("CGW","XLE", "SPY")

# Download these tickers from Yahoo for the dates in the presentation
getSymbols(tickers,src="yahoo", from = "2007-06-01",to = "2018-01-26") # original paper dates.

# Merge all the Price series into one dataframe
AllPrices <- do.call(merge, lapply(tickers, function(x) get(x)))

#Import Data from .RData file - this was the data from google
# If you try and download from google now this is the warning: "Error: 'getSymbols.google' is defunct.
# Google Finance stopped providing data in March, 2018."
# New Issue: https://github.com/joshuaulrich/quantmod/issues/221
#load("Water-Energy.Rdata")

#Some of these series have (NA) missing values for dates when others
# do not have missiong vaulesin the series so we interpolate for these values
AllPrices$CGW.Close <- interpNA(AllPrices$CGW.Close)
AllPrices$XLE.Close <- interpNA(AllPrices$XLE.Close)
AllPrices$SPY.Close <- interpNA(AllPrices$SPY.Close)

#Set up the correct data frame for prices
CGW <- as.data.frame((AllPrices$CGW.Close))
XLE <- as.data.frame((AllPrices$XLE.Close))
SPY <- as.data.frame((AllPrices$SPY.Close))

#Retrieve specific dates for this time frame (prices)
date1 <- index(AllPrices)

# Set Up Data Frame Prices

Data <- cbind(CGW[,1], XLE[,1], SPY[,1])

# WaveL2E Analysis of Prices

periodic_waveL2E <- WaveL2E(Data[,1], date = date1, block = 1)
periodic_waveL2E_2 <- WaveL2E(Data[,2], date= date1, block = 1)
periodic_waveL2E_3 <- WaveL2E(Data[,3], date= date1, block = 1)

#Set up the correct data frame for returns
rCGW <- as.data.frame(returns(AllPrices$CGW.Close))
rXLE <- as.data.frame(returns(AllPrices$XLE.Close))
rSPY <- as.data.frame(returns(AllPrices$SPY.Close))

#Retrieve specific dates for this time frame
date1 <- index(AllPrices)

# Set Up Data Frame Returns

Data <- cbind(rCGW[-1,1], rXLE[-1,1], rSPY[-1,1])

# WaveL2E Analysis of Returns

periodic_waveL2E_R <- WaveL2E(Data[,1], date = date1[-1], block = 1)
periodic_waveL2E_2_R <- WaveL2E(Data[,2], date= date1[-1], block = 1)
periodic_waveL2E_3_R <- WaveL2E(Data[,3], date= date1[-1], block = 1)

######## Figure 1 ##########

CoFESWave.image(periodic_waveL2E$original, periodlab = " ",timelab = "  " , main = " ",
         legend.params = list(lab = "wavelet power levels", mar = 5.1, cex = 4, n.ticks = 5),
         color.key = "quantile", lwd = 2, plot.ridge = FALSE, color.palette = "topl_colors(n.levels)")

CoFESWave.image(periodic_waveL2E_2$original, periodlab = " ",timelab = "  " , main = " ",
         legend.params = list(lab = "wavelet power levels", mar = 5.1, cex = 4, n.ticks = 5),
         color.key = "quantile", lwd = 2, plot.ridge = FALSE, color.palette = "topl_colors(n.levels)")

CoFESWave.image(periodic_waveL2E_R$original, periodlab = " ",timelab = "  " , main = " ",
         legend.params = list(lab = "wavelet power levels", mar = 5.1, cex = 4, n.ticks = 5),
         color.key = "quantile", lwd = 2, plot.ridge = FALSE, color.palette = "topl_colors(n.levels)")

CoFESWave.image(periodic_waveL2E_2_R$original, periodlab = " ",timelab = "  " , main = " ",
         legend.params = list(lab = "wavelet power levels", mar = 5.1, cex = 4, n.ticks = 5),
         color.key = "quantile", lwd = 2, plot.ridge = FALSE, color.palette = "topl_colors(n.levels)")

####### Figure 2 & Table 1 #############

#Sinusoids (periodic)
# Quarterly Investment Horizon
x1 <- periodic.series(start.period = 64, length = 2400)
# Biweekly Investment Horizon
x2 <- periodic.series(start.period = 15, length = 2400)
# Biannual Investment Horizon
x3 <- periodic.series(start.period = 125, length = 2400)
#Add up all the contributions
periodic_func <- x3 + x1 + x2 + rnorm(2400, mean = 0, sd = 0.1)
periodic_func_2 <- x3 + x1 + x2
Data <-cbind(periodic_func, periodic_func_2)
X_1 <- Data[,1]

# Static Implementation
figure_2 <- WaveL2E(X_1, block = length(X_1))

# Table 1 #

joint_P = rbind(cbind(figure_2$PTV_L2E, figure_2$PSL_L2E),cbind(figure_2$PTV_Chi_square, figure_2$PSL_Chi_square))
colnames(joint_P) <- c("PTV", "PSL")
rownames(joint_P) <- c("L2E","L2E Chi-Sq")
joint_P

####### Figure 3 & 4 ###############

## Not in R

####### Figure 5 & Table 1 + 2 ##################
#Dynamic Implementation
figure_5 <- WaveL2E(X_1, block = 1)

# Table 1 #

joint_P = rbind(cbind(figure_5$PTV_L2E, figure_5$PSL_L2E),cbind(figure_5$PTV_Chi_square, figure_5$PSL_Chi_square))
colnames(joint_P) <- c("PTV", "PSL")
rownames(joint_P) <- c("L2E","L2E Chi-Sq")
joint_P

# Table 2 #

#Comaprison Caluations

#Universal Hard and Soft Thershold
# https://www.rdocumentation.org/packages/wmtsa/versions/2.0-3/topics/wavShrink
hard_uni_1 <- wavShrink(X_1, wavelet="s8", shrink.fun="hard",
                        thresh.fun="universal", threshold=NULL,thresh.scale=1,
                        xform="dwt", noise.variance=-1.0, reflect=TRUE)
soft_uni_1 <- wavShrink(X_1, wavelet="s8", shrink.fun="soft",
                        thresh.fun="universal", threshold=NULL,thresh.scale=1,
                        xform="dwt", noise.variance=-1.0, reflect=TRUE)

#Ebayes Implementation - Empirical Bayes thresholding on the levels of a wavelet transform
# https://github.com/stephenslab/EbayesThresh/blob/master/R/ebayesthresh.wavelet.R
#Wavelet Transform
bayes.dwt <- dwt(X_1,  wf="la8")
#Emperical Bayesian Threshold
bayes_imp <- ebayesthresh.wavelet(bayes.dwt)
#nverse Transform
bayes_smooth <- idwt(bayes_imp)

# Comparison Results
rmse_Compare <- as.data.frame(cbind(round(rmse(Data[,2],hard_uni_1), 4),
                                    round(rmse(Data[,2],soft_uni_1), 4),
                                    round(rmse(Data[,2],bayes_smooth), 4),
                                    round(rmse(Data[,2],figure_5$recon_L2E$series$x.r), 4),
                                    round(rmse(Data[,2],figure_5$recon_Chi_square$series$x.r), 4),
                                    round(rmse(Data[,2],figure_2$recon_L2E$series$x.r), 4),
                                    round(rmse(Data[,2],figure_2$recon_Chi_square$series$x.r), 4)))
colnames(rmse_Compare) <- c("Uni Hard Thresh", "Uni Soft Thresh", "EBayes Thresh",
                            "WaveL2E Thresh (D)", "WaveL2E_C2 Thresh (D)", "WaveL2E Thresh (S)", "WaveL2E_C2 Thresh (S)")
rmse_Compare

############## Figure 6 ##############

#simulation from CWT: A Primer Paper  (adapted)
#10 year cycle (first periodic component) the second periodic component 1 year cycle that changes to a 3 yr cycle between
#10-20 years and 5 years cycle between 50-60 years adn with variance that changes 5 times (every 12 years)
t <- seq(from = 1/12, to = 100, by = 1/12) #60 years monthly data # seq(from = 1/10, to = 24000, by = 1)#seq(from = 1/12, to = 300, by = 1/12) #
#how many breaks in variance do we want
inc <- 5
set <- floor(length(t)/inc)
maxV <- 0.5
minV <- 0
eps_t <- NULL
ran_sig <- runif(inc, min = minV, max = maxV)
set.seed(24)
for(i in 1:inc){
  eps_t <- c(eps_t, rnorm(set, mean = 0, sd = ran_sig[i]))
}

eps_t <- c(rnorm(length(t)/5, mean = 0, sd = 0.5), rnorm(length(t)/5, mean = 0, sd = 0.05),rnorm(length(t)/5, mean = 0, sd = 0.25), rnorm(length(t)/5, mean = 0, sd = 0.15), rnorm(length(t)/5, mean = 0, sd = 0.05)) # N(0,1) noise

plot(density(eps_t))

Y_t <- NULL
Y_t2 <- NULL
for(i in 1:length(t)){
  if(t[i]>= 10 && t[i]<= 20){
    Y_t[i] <- cos(2*pi*t[i]/10) + cos(2*pi*t[i]/3) + eps_t[i]
    Y_t2[i] <- cos(2*pi*t[i]/10) + cos(2*pi*t[i]/3)
  }
  if(t[i]>= 50 && t[i]<= 60){
    Y_t[i] <- cos(2*pi*t[i]/10) + cos(2*pi*t[i]/5) + eps_t[i]
    Y_t2[i] <- cos(2*pi*t[i]/10) + cos(2*pi*t[i]/5)
  }
  else{
    Y_t[i] <- cos(2*pi*t[i]/10) + cos(2*pi*t[i]/1) + eps_t[i]
    Y_t2[i] <- cos(2*pi*t[i]/10) + cos(2*pi*t[i]/1)
  }
}

#############

Data <-cbind(Y_t, Y_t2)

X_1 <- Data[,1]


# Comaprison Caluations

# Thresholding Methods
#Limitation is that the time - series have to be some number with 2^power
#Error: Sample size is not divisible by 2^J

#Universal Hard and Soft Thershold
# https://www.rdocumentation.org/packages/wmtsa/versions/2.0-3/topics/wavShrink
hard_uni_1 <- wavShrink(X_1, wavelet="s8", shrink.fun="hard",
                        thresh.fun="universal", threshold=NULL,thresh.scale=1,
                        xform="dwt", noise.variance=-1.0, reflect=TRUE)
soft_uni_1 <- wavShrink(X_1, wavelet="s8", shrink.fun="soft",
                        thresh.fun="universal", threshold=NULL,thresh.scale=1,
                        xform="dwt", noise.variance=-1.0, reflect=TRUE)

#Ebayes Implementation - Empirical Bayes thresholding on the levels of a wavelet transform
# https://github.com/stephenslab/EbayesThresh/blob/master/R/ebayesthresh.wavelet.R
#Wavelet Transform
bayes.dwt <- dwt(X_1,  wf="la8")
#Emperical Bayesian Threshold
bayes_imp <- ebayesthresh.wavelet(bayes.dwt)
#nverse Transform
bayes_smooth <- idwt(bayes_imp)

# Static Implementation
block = length(X_1)
figure_6 <- WaveL2E(X_1, block = block)

# Comparison Results
signal_recovered <- figure_6$Emp_WaveL2E
rmse_Compare <- as.data.frame(cbind(round(rmse(Data[,2],hard_uni_1), 4),
                                    round(rmse(Data[,2],soft_uni_1), 4),
                                    round(rmse(Data[,2],bayes_smooth), 4),
                                    round(rmse(Data[,2],figure_6$recon_L2E$series$x.r), 4),
                                    round(rmse(Data[,2],figure_6$recon_Chi_square$series$x.r), 4),
                                    round(rmse(Data[,2],signal_recovered), 4)))
colnames(rmse_Compare) <- c("Uni Hard Thresh", "Uni Soft Thresh", "EBayes Thresh",
                            "WaveL2E Thresh", "WaveL2E_C2 Thresh", "EWaveL2E (w, sigma)")
rmse_Compare
joint_P = rbind(cbind(figure_6$PTV_L2E, figure_6$PSL_L2E),cbind(figure_6$PTV_Chi_square, figure_6$PSL_Chi_square))
colnames(joint_P) <- c("PTV", "PSL")
rownames(joint_P) <- c("L2E","L2E Chi-Sq")
joint_P


###### Dynamic implementation - block size 1
block = 1
figure_6_2 <- WaveL2E(X_1, block = block)
# Comparison Results
signal_recovered <- figure_6_2$Emp_WaveL2E
rmse_Compare <- as.data.frame(cbind(round(rmse(Data[,2],hard_uni_1), 4),
                                    round(rmse(Data[,2],soft_uni_1), 4),
                                    round(rmse(Data[,2],bayes_smooth), 4),
                                    round(rmse(Data[,2],figure_6_2$recon_L2E$series$x.r), 4),
                                    round(rmse(Data[,2],figure_6_2$recon_Chi_square$series$x.r), 4),
                                    round(rmse(Data[,2],signal_recovered), 4)))
colnames(rmse_Compare) <- c("Uni Hard Thresh", "Uni Soft Thresh", "EBayes Thresh",
                            "WaveL2E Thresh", "WaveL2E_C2 Thresh", "Backsolve (w_t, sigma_t)")
rmse_Compare

joint_P = rbind(cbind(figure_6_2$PTV_L2E, figure_6_2$PSL_L2E),cbind(figure_6_2$PTV_Chi_square, figure_6_2$PSL_Chi_square))
colnames(joint_P) <- c("PTV", "PSL")
rownames(joint_P) <- c("L2E","L2E Chi-Sq")
joint_P

################ Figure 7 + Table 5 & 6 ###################

all_signals <- c("hisine","losine", "linchirp", "twochirp", "quadchirp","mishmash1", "mishmash2", "mishmash3")
z <- lapply(all_signals, make.signal)
topl_colors <- grDevices:::colorRampPalette(c("#64d8cb","#26a69a", "#90a4ae","#5f5fc4", "#283593"))
ifultools::stackPlot(x=seq(1024),y=z, ylab= list(text = all_signals, col = topl_colors(length(all_signals))), col = topl_colors(length(all_signals)))

signal_test <- function(signal_to_noise, series_length, block){

  set.seed(123)
  rmse_Compare <- NULL
  joint_P_Compare <- NULL
  row_names <- NULL
  all_signals <- c("hisine","losine", "linchirp", "twochirp", "quadchirp","mishmash1", "mishmash2", "mishmash3")
  for(j in 1: length(all_signals)){
    signal_noise <- make.signal(all_signals[j], n=series_length, snr=signal_to_noise)
    signal_noN <- make.signal(all_signals[j], n=series_length, snr=Inf)
    Data <- cbind(signal_noise@data, signal_noN@data)

    plot(signal_noise)
    plot(signal_noN)

    X_1 <- Data[,1]

    periodic_waveL2E <- WaveL2E(X_1, block = block)

    # Comaprison Caluations

    # Thresholding Methods
    #Limitation is that the time - series have to be some number with 2^power
    #Error: Sample size is not divisible by 2^J

    #Universal Hard and Soft Thershold
    # https://www.rdocumentation.org/packages/wmtsa/versions/2.0-3/topics/wavShrink
    hard_uni_1 <- wavShrink(X_1, wavelet="s8", shrink.fun="hard",
                            thresh.fun="universal", threshold=NULL,thresh.scale=1,
                            xform="dwt", noise.variance=0, reflect=TRUE)
    soft_uni_1 <- wavShrink(X_1, wavelet="s8", shrink.fun="soft",
                            thresh.fun="universal", threshold=NULL,thresh.scale=1,
                            xform="dwt", noise.variance=0, reflect=TRUE)

    #Ebayes Implementation - Empirical Bayes thresholding on the levels of a wavelet transform
    # https://github.com/stephenslab/EbayesThresh/blob/master/R/ebayesthresh.wavelet.R
    #Wavelet Transform
    bayes.dwt <- dwt(X_1,  wf="la8")
    #Emperical Bayesian Threshold
    bayes_imp <- ebayesthresh.wavelet(bayes.dwt)
    #nverse Transform
    bayes_smooth <- idwt(bayes_imp)

    # Comparison Results
    # Recover the series using L2E back calculation
    signal_recovered <- periodic_waveL2E$Emp_WaveL2E_MAD
    rmse_Compare_new <- as.data.frame(cbind(round(rmse(Data[,2],hard_uni_1), 4),
                                            round(rmse(Data[,2],soft_uni_1), 4),
                                            round(rmse(Data[,2],bayes_smooth), 4),
                                            round(rmse(Data[,2],periodic_waveL2E$recon_L2E$series$x.r), 4),
                                            round(rmse(Data[,2],periodic_waveL2E$recon_Chi_square$series$x.r), 4),
                                            round(rmse(Data[,2],signal_recovered), 4)))
    colnames(rmse_Compare_new) <- c("Uni Hard Thresh", "Uni Soft Thresh", "EBayes Thresh",
                                    "WaveL2E Thresh", "WaveL2E_C2 Thresh", "EWaveL2E")
    rmse_Compare <- rbind(rmse_Compare, rmse_Compare_new)
    joint_P_new <- rbind(cbind(periodic_waveL2E$PTV_L2E, periodic_waveL2E$PSL_L2E),cbind(periodic_waveL2E$PTV_Chi_square, periodic_waveL2E$PSL_Chi_square))
    colnames(joint_P_new) <- c("PTV", "PSL")
    row_names <- rbind(row_names, rbind(paste("L2E", all_signals[j]),paste("L2E Chi-Sq",all_signals[j])))
    joint_P_Compare <- rbind(joint_P_Compare, joint_P_new)
  }

  rownames(rmse_Compare) <- all_signals
  rownames(joint_P_Compare) <- row_names

  return(list(rmse <- rmse_Compare, percentages <- joint_P_Compare))
}

test_table_2 <- signal_test(signal_to_noise = 2, series_length = 1024, block = 1)
test_table_2[[1]]
test_table_5 <- signal_test(signal_to_noise = 5, series_length = 1024, block = 1)
test_table_5[[1]]

############### Figure 8 ####################

### Load all the transformation and plotting functions for coherence and partial coherence

# plot_functions
# functions_transform

Data <- cbind(CGW[,1], XLE[,1], SPY[,1])

# ----- Create matrix with columns ETF columns  ----- #
T <- length(CGW[,1])
X<-matrix(Data,T,3)

# -----  Choice of wavelet parameters  ------ #
low.period<-1
up.period<-512

# -----   Choice of smoothing parameters for coherency ---- #
wt.type ='ham'
wt.size =3
ws.type ='ham'
ws.size =3

# Computation of coherency
# WCO<-CoFEScoherency(Data[,1],Data[,2],low.period=low.period,up.period=up.period,low.fp = 32,up.fp = 128)

# --- Lower and upper  periods (mid regime)
lowFP1<-32
upFP1<-128

WCO<-CoFEScoherency(Data[,1],Data[,2],low.period=low.period,up.period=up.period,
                    low.fp = lowFP1,up.fp = upFP1, Phase_diff = TRUE, date = date1)

# -------- PLOTS  ----------
plot_CoFESWaveWCO(WCO,pE=5, horizons.label = FALSE)


# Computation of coherency of (x, y)
WCO_recon<-CoFEScoherency(periodic_waveL2E_2$recon_L2E$series$x.r,
                          periodic_waveL2E$recon_L2E$series$x.r,
                          low.fp = lowFP1,up.fp = upFP1, Phase_diff = TRUE, date = date1)

plot_CoFESWaveWCO(WCO_recon,pE=5, horizons.label = FALSE)


# ----- Computation of Partial Coherency of Water , Energy (controlling for SPY) ---- #

index.p=2
Data_thresh <- cbind(periodic_waveL2E_2$recon_L2E$series$x.r,periodic_waveL2E$recon_L2E$series$x.r, periodic_waveL2E_3$recon_L2E$series$x.r)
X<-matrix(Data_thresh,T,3)
W<-CoFESmpcoherency(X, coher.type='part',index.p=index.p,
                 low.period=low.period,up.period=up.period,
                 low.fp = lowFP1,up.fp = upFP1, Phase_diff = TRUE, date = date1)

plot_CoFESWaveWCO(W,pE=5)

########## Figure 9 ############

# 22 Days - Fiscal Month

periodic_waveL2E_22 <- WaveL2E(XLE[,1], date= date1, block = 22)
periodic_waveL2E_22_R <- WaveL2E(rXLE[-1,1], date= date1[-1], block = 22)

plot(periodic_waveL2E_22$sig,type="l",col=topl_colors(3)[3],lwd=2.0,
     xlab=" ",ylab="Variance",main="Noise Variance",las=3, xaxt='n')
axis(side=1, at=c(1,120*seq(1,22,1)),lab=date1[c(1,120*seq(1,22,1))], las = 2)
grid(120, NA, lwd = 1)

plot(periodic_waveL2E_22$w,type="l",col=topl_colors(3)[3],lwd=2.0,
     xlab=" ",ylab="Weight Values",main="Threshold Weights",las=3, xaxt='n')
axis(side=1, at=c(1,120*seq(1,22,1)),lab=date1[c(1,120*seq(1,22,1))], las = 2)
grid(120, NA, lwd = 1)

df.xts <- (cbind(periodic_waveL2E_22$recon_L2E$series$x, (periodic_waveL2E_22$recon_L2E$series$x.r)))
df.xts <- as.xts(df.xts, order.by = date1)
colnames(df.xts) <- (c("Original Series", "WaveL2E"))
plot(df.xts["2011-07-01/2012-07-01"], col = topl_colors(2),
     lwd = c(1,2), lty = c(1,2), grid.col = NA, yaxis.right = FALSE,
     legend.loc = "bottomright", auto.legend=TRUE,
     main=" ", major.ticks = "months", ylab = "Prices (logscale)")

df.xts <- (cbind(periodic_waveL2E_22_R$recon_L2E$series$x, (periodic_waveL2E_22_R$Emp_WaveL2E),periodic_waveL2E_22_R$recon_L2E$series$x.r))
df.xts <- as.xts(df.xts, order.by = date1[-1])
colnames(df.xts) <- (c("Original Series", "EWaveL2E", "WaveL2E"))
plot(df.xts["2011-07-01/2012-07-01"], col = topl_colors(3),
     lwd = c(1,2,2), lty = c(1,2,6), grid.col = NA, yaxis.right = FALSE,
     legend.loc = "bottomright", auto.legend=TRUE,
     main=" ", major.ticks = "months", ylab = "Prices (logscale)")

