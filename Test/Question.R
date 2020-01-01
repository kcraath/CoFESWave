############################
##
## Example Script
##
############################

# All the libraries

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


# Sector ETF's Behavior

## # A tibble: 10 × 2
##    ticker                 sector
##     <chr>                  <chr>
## 1     XLY Consumer Discretionary
## 2     XLP       Consumer Staples
## 3     XLE                 Energy
## 4     XLF             Financials
## 5     XLV            Health Care
## 6     XLI            Industrials
## 7     XLB              Materials
## 8     XLK Information Technology
## 9     XLU              Utilities
## 10    SPY                  Index


# Identify the tickers of interest
tickers <- c("XLY", "XLP","XLE", "XLF", "XLV", "XLI", "XLB", "XLK", "XLU", "SPY")

# Download these tickers from Yahoo for the dates in the presentation
getSymbols(tickers,src="yahoo", from = "1999-01-01",to = "2019-01-01")

# Merge all the Price series into one dataframe
AllPrices <- do.call(merge, lapply(tickers, function(x) get(x)))

AllPrices$XLY.Close <- interpNA(AllPrices$XLY.Close)
AllPrices$XLP.Close <- interpNA(AllPrices$XLP.Close)
AllPrices$XLE.Close <- interpNA(AllPrices$XLE.Close)
AllPrices$XLF.Close <- interpNA(AllPrices$XLF.Close)
AllPrices$XLV.Close <- interpNA(AllPrices$XLV.Close)
AllPrices$XLI.Close <- interpNA(AllPrices$XLI.Close)
AllPrices$XLB.Close <- interpNA(AllPrices$XLB.Close)
AllPrices$XLK.Close <- interpNA(AllPrices$XLK.Close)
AllPrices$XLU.Close <- interpNA(AllPrices$XLU.Close)
AllPrices$SPY.Close <- interpNA(AllPrices$SPY.Close)

rXLY <- as.data.frame((AllPrices$XLY.Close))
rXLP <- as.data.frame((AllPrices$XLP.Close))
rXLE <- as.data.frame((AllPrices$XLE.Close))
rXLF <- as.data.frame((AllPrices$XLF.Close))
rXLV <- as.data.frame((AllPrices$XLV.Close))
rXLI <- as.data.frame((AllPrices$XLI.Close))
rXLB <- as.data.frame((AllPrices$XLB.Close))
rXLK <- as.data.frame((AllPrices$XLK.Close))
rXLU <- as.data.frame((AllPrices$XLU.Close))
rSPY <- as.data.frame((AllPrices$SPY.Close))

date_Sector <- index(AllPrices)

#----- Set Up Data Frame -----#


Data <-cbind(rXLY[,1], rXLP[,1], rXLE[,1], rXLF[,1], rXLV[,1], rXLI[,1], rXLB[,1], rXLK[,1], rXLU[,1], rSPY[,1])

###### Plot Series
topl_colors <- colorRampPalette(c("#64d8cb","#26a69a","#5f5fc4", "#283593"))
ts.plot(Data, col = topl_colors(10))

# Options 1 - without SVD ##########

WaveL2E(Data[,1], date = date_Sector, block = 50)



# Run WaveL2E
S1 <- WaveL2E_dates(Data[,1], date_Sector)
S2 <- WaveL2E_dates(Data[,2], date_Sector)
S3 <- WaveL2E_block(Data[,3], block = 50)
S4 <- WaveL2E_block(Data[,4], block = 100)
S5 <- WaveL2E(Data[,5])
S6 <- WaveL2E(Data[,6])
S7 <- WaveL2E(Data[,7])
S8 <- WaveL2E(Data[,8])
S9 <- WaveL2E(Data[,9])
S10 <- WaveL2E(Data[,10])

# Options 2 - use the other series by doing a SVD ##########

newData <- svd(Data)
Data <- Data%*%(newData$v)

# Run WaveL2E
S1 <- WaveL2E_dates(Data[,1], date_Sector)
S2 <- WaveL2E_dates(Data[,2], date_Sector)
S3 <- WaveL2E_dates(Data[,3], date_Sector)
S4 <- WaveL2E_dates(Data[,4], date_Sector)
S5 <- WaveL2E_dates(Data[,5], date_Sector)
S6 <- WaveL2E_dates(Data[,6], date_Sector)
S7 <- WaveL2E_dates(Data[,7], date_Sector)
S8 <- WaveL2E_dates(Data[,8], date_Sector)
S9 <- WaveL2E_dates(Data[,9], date_Sector)
S10 <- WaveL2E_dates(Data[,10], date_Sector)

###### PTV & PSL ####

Nexus_density <- rbind(cbind(S1$PTV, S2$PSL),cbind(S2$PTV, S2$PSL), cbind(S3$PTV, S3$PSL), cbind(S4$PTV, S4$PSL),
                       cbind(S5$PTV, S5$PSL), cbind(S6$PTV, S6$PSL), cbind(S7$PTV, S7$PSL), cbind(S8$PTV, S8$PSL), cbind(S9$PTV, S9$PSL), cbind(S10$PTV, S10$PSL)) 

colnames(Nexus_density) <- c("PTV", "PSL")
rownames(Nexus_density) <- c("XLY L2E","XLY L2E Chi-Sq", "XLP L2E", "XLP L2E Chi-sq","XLE L2E", "XLE L2E Chi-sq", "XLF L2E", "XLF L2E Chi-sq", 
                             "XLV L2E", "XLV L2E Chi-sq", "XLI L2E", "XLI L2E Chi-sq", "XLB L2E", "XLB L2E Chi-sq", "XLK L2E", "XLK L2E Chi-sq", "XLU L2E", "XLU L2E Chi-sq", "SPY L2E", "SPY L2E Chi-sq") 
Nexus_density

###########

# plot the level density (intra-scale Density)
plot(density(log(S1$dis0[1,])))
#plot the time block density
plot(density(log(S1$dis0[,100])), main = "Inter-scale Density")
abline(v=c(log(S1$thresh[100]), log(S1$qthresh[100])), col = c("red", "blue"))
#Plot the variance
plot.ts(log(S1$original$series$x[100:115]), main = "Time Series +/- Variance", ylab = "Observation Values", xlab = "Observation Index")
lines(log(S1$original$series$x[100:115]) + 8*(S1$sig[100:115]), col = "red", lty = 3) 
lines(log(S1$original$series$x[100:115]) - 8*(S1$sig[100:115]), col = "red", lty = 3) 

plot.ts(log(S1$sig), main = "Time Series +/- Variance", ylab = "Observation Values", xlab = "Observation Index")


mean((S1$sig) - mad(S1$original$Power[,1])/0.675)

###### Simulation Example ########

#simulation from CWT: A Primer Paper  (adapted)
#10 year cycle (first periodic component) the second periodic component 1 year cycle that changes to a 3 yr cycle between 
#10-20 years and 5 years cycle between 50-60 years adn with variance that changes 5 times (every 12 years)
t <- seq(from = 1/12, to = 60, by = 1/12) #60 years monthly data # seq(from = 1/10, to = 24000, by = 1) #
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

###### interscale AND intra-scale implementation

block = 50 #has to be > 2
periodic_waveL2E <- WaveL2E_block(Data[,1], block = block) #WaveL2E(Data[,1]) # 
noise <- NULL
# Comparison Results
# Recover the series using L2E back calculation
for(i in 1:(floor(length(Data[,1])/block)-1)){
  noise <- c(noise, periodic_waveL2E$w[i]*rnorm(n = block, mean = 0, sd = (periodic_waveL2E$sig[i]))/(1-periodic_waveL2E$w[i]))
}
signal_recovered <- (X_1[1:length(noise)] - noise)

plot.ts(signal_recovered[100:200], col = "#5f5fc4")
lines(Data[100:200,2])
rmse(Data[100:200,2], signal_recovered[100:200])

rmse_Compare <- as.data.frame(cbind(round(rmse(Data[1:length(noise),2],hard_uni_1[1:length(noise)]), 4), 
                                    round(rmse(Data[1:length(noise),2],soft_uni_1[1:length(noise)]), 4), 
                                    round(rmse(Data[1:length(noise),2],bayes_smooth[1:length(noise)]), 4), 
                                    round(rmse(Data[1:length(noise),2],periodic_waveL2E$recon1$series$x.r[1:length(noise)]), 4), 
                                    round(rmse(Data[1:length(noise),2],periodic_waveL2E$recon2$series$x.r[1:length(noise)]), 4), 
                                    round(rmse(Data[1:length(noise),2],signal_recovered), 4)))
colnames(rmse_Compare) <- c("Uni Hard Thresh", "Uni Soft Thresh", "EBayes Thresh", 
                            "WaveL2E Thresh", "WaveL2E_C2 Thresh", "Backsolve (w, sigma)")
rmse_Compare

par(mfrow=c(1,1))

plot(X_1[100:300], type = "l", lty = 1, main = "Accuracy (RMSE) of Threshold Methods (t=100,..., 300)", col = rainbow10equal[10], ylab = " ", xlab = "Time",  ylim =c(-6,4))

lines(Data[100:300,2], type = "l", lty = 1, main = "Original Series", col = "red", ylab = " ", xlab = "Time")

lines(signal_recovered[100:300], type = "l",lty = 2 , col = "blue" , main = "Recovered", ylab = " ")

lines(hard_uni_1[100:300], type = "l",lty = 3 , col = rainbow10equal[1], main = "Hard Universal", ylab = " ")

lines(soft_uni_1[100:300], type = "l",lty = 3,  col = rainbow10equal[2], main = "Soft Universal", ylab = "")

lines(bayes_smooth[100:300], type = "l",lty = 2, col = rainbow10equal[3], main = "Bayes Thresh",ylab = "")

lines(periodic_waveL2E$recon1$series$x.r[100:300], col = rainbow10equal[5], 
      main = c("WavL2E Thresh to Zeroes"), ylab = "", xlab = "index", lty = 2)

lines(periodic_waveL2E$recon2$series$x.r[100:300], col = rainbow10equal[6], 
      main = c("WavL2E Thresh using Chi-Square"), ylab = "", xlab = "Index", lty = 2)

abline(v = c(44, 140, 188) , col = topl_colors(3), lty = 3)

legend("bottomleft", c("Original series", "Original w/out Noise", "Recovered (L2E) Original Series", paste("Bayes Thresh", paste(rmse_Compare$`EBayes Thresh`)), paste("Cleaned Up (1)", paste(rmse_Compare$`WaveL2E Thresh`)), paste("Cleaned Up (2)", paste(rmse_Compare$`WaveL2E_C2 Thresh`))), col = c(rainbow10equal[10],"red", "blue", rainbow10equal[3], rainbow10equal[5:6]), lty = c(1,1,2,2,2,2), bty='n')



#### Clustering
# Identify the tickers of interest
tickers <- c("CGW","XLE", "SPY")

# Download these tickers from Yahoo for the dates in the presentation
getSymbols(tickers,src="yahoo", from = "2007-06-01",to = "2018-01-26")

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

# Wavelet Analysis

# Most of these simulations are run with a very low n-number so that results can be viewed quicker. 
# Increase the n-number for more reliable results

#Set up the correct data frame
rCGW <- as.data.frame((AllPrices$CGW.Close))
rXLE <- as.data.frame((AllPrices$XLE.Close))
rSPY <- as.data.frame((AllPrices$SPY.Close))

#Retrieve specific dates for this time frame
date1 <- index(AllPrices)

# Set Up Data Frame #

Data <- cbind(rCGW[,1], rXLE[,1], rSPY[,1])

periodic_waveL2E <- WaveL2E_dates(Data[,1], date = date1)
periodic_waveL2E_2 <- WaveL2E_dates(Data[,2], date= date1) #WaveL2E_block(Data[,2], block = block)

MV1 <- periodic_waveL2E$dis0
MV2 <- periodic_waveL2E_2$dis0

#3D Plots of Scale Specific Combinations between series

library(SDMTools)
new <- NULL
new[1] <- 1

#scale_data <- c(0,2,4,8,16,32,64,128) #needs to be some sequence of 2^J

levels <- seq(from = 1, to = floor(log(length(MV1[,1]),2)), 1)
scale_data <- c(0, 2^levels) 

# Zero padding to size=next power of 2+2
#	pot2 <- ceiling(log2(nTimes))
#	extra.length <- 2^(pot2+2)-nTimes  


for(i in 2:length(scale_data)){
  new[i] <- which(floor(periodic_waveL2E_2$original$Scale) > scale_data[i])[1]
}

#Cluster implementation

MV1_new <- NULL
MV2_new <- NULL
for(i in 2:length(scale_data)){
  MV1_new <- rbind(MV1_new, cbind(MV1[new[i-1]:new[i],], scale_data[i]))
  MV2_new <- rbind(MV2_new, cbind(MV2[new[i-1]:new[i],], scale_data[i]))
  
  data_length <- (length(MV1[1,])+1)
  data_chek1 <- MV1_new[MV1_new[,data_length] == scale_data[i],-data_length]
  data_chek2 <- MV2_new[MV2_new[,data_length] == scale_data[i],-data_length]
  df <- cbind((apply(data_chek1, 2, FUN = sum)), (apply(data_chek2,2, FUN = sum)))
  clusters <- kmeans(df[,1:2], 2, nstart = 10)
  clus_data <- cbind(df, as.factor(clusters$cluster))
  
  periodic_waveL2E$Ana_Wave$Wave[MV1_new[,data_length] == scale_data[i],clus_data[,3] == 2] <- 0
  periodic_waveL2E$Ana_Wave$Power[MV1_new[,data_length] == scale_data[i],clus_data[,3] == 2] <- 0
  periodic_waveL2E_2$Ana_Wave$Wave[MV2_new[,data_length] == scale_data[i],clus_data[,3] == 2] <- 0
  periodic_waveL2E_2$Ana_Wave$Power[MV2_new[,data_length] == scale_data[i],clus_data[,3] == 2] <- 0
}

#Plot Reconstructed Time-Seires

recon_org1 <- WaveletComp::reconstruct(periodic_waveL2E$Ana_Wave, plot.waves = FALSE, lwd = c(1,2),
                                       legend.coords = "topright",plot.rec = FALSE, verbose = F)

wt.image(periodic_waveL2E$Ana_Wave, periodlab = " ",timelab = "  " , main = " ",
         legend.params = list(lab = "wavelet power levels", mar = 5.1, cex = 4, n.ticks = 10), 
         color.key = "quantile", lwd = 2, plot.ridge = FALSE, color.palette = "topl_colors(n.levels)")


recon_org2 <- WaveletComp::reconstruct(periodic_waveL2E_2$Ana_Wave, plot.waves = FALSE, lwd = c(1,2),
                                       legend.coords = "topright",plot.rec = FALSE, verbose = F)

wt.image(periodic_waveL2E_2$Ana_Wave, periodlab = " ",timelab = "  " , main = " ",
         legend.params = list(lab = "wavelet power levels", mar = 5.1, cex = 4, n.ticks = 5), 
         color.key = "quantile", lwd = 2, plot.ridge = FALSE, color.palette = "topl_colors(n.levels)")

#Wavelet Squared Coherence
recon_1 <- cbind(1:(length(recon_org1$series$x.r)), recon_org1$series$x.r)
recon_2 <- cbind(1:(length(recon_org2$series$x.r)), recon_org2$series$x.r)

# or 

recon_1 <- cbind(1:(length(periodic_waveL2E$recon1$series$x.r)), periodic_waveL2E$recon1$series$x.r)
recon_2 <- cbind(1:(length(periodic_waveL2E_2$recon1$series$x.r)), periodic_waveL2E_2$recon1$series$x.r)


wtc.12=wtc(recon_1, recon_2, quiet = TRUE, nrands = 100, mother = "morlet")

#Add the dates to the axis of the squared coherence plot
wtc.12$xaxis <- date1
par(oma=c(0, 0, 0, 1), mar=c(5, 5, 5, 5) + 0.1)

#Plotting the wavelet Squared Coherence
plot(wtc.12, plot.cb=TRUE, plot.phase=TRUE, xlab = " ", ylab = " ", cex = 1.6, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, cex.lab = 1.8, fill.cols= topl_colors(n))

#Add annual lines and lines to distinguish between investment horizons
n = length(Data[, 1])
abline(v = seq(12, n, 12), h = 1:16, col = "brown", lty = 1, lwd = 1)
title("Water-Energy Nexus", cex.main = 1.8)

