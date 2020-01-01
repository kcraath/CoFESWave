############################
##
## Methods Paper Figures
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

#Latex expersions
library(latex2exp)

#Latex tables
library(knitr)

#Interactive Charts
library(highcharter)

#### Figure 1 #####

# Identify the tickers of interest
tickers <- c("CGW","XLE", "SPY")

# Download these tickers from Yahoo for the dates in the presentation
getSymbols(tickers,src="yahoo", from = "2007-06-01",to = "2019-10-31")

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

# Data <- cbind(rCGW[,1], rXLE[,1], rSPY[,1])
# 
# Data <- AllPrices[, c("CGW.Close", "XLE.Close","SPY.Close")]





###################### Try plotWCO_Yifan

## Before you run them, please read my notes in plotWCO_Yifan first.
Data <- AllPrices[, c("CGW.Close", "XLE.Close")]
plotWCO_Yifan(Data)
plotWCO_Yifan(Data, L2E = TRUE)

Data <- AllPrices[, c("CGW.Close", "XLE.Close","SPY.Close")]
plotWCO_Yifan(Data)
plotWCO_Yifan(Data, L2E = TRUE)


index(Data)


d1 <- as.matrix(index(Data),Data[,1])


d1 <- cbind(1:3127, Data[,1])
d2 <- cbind(1:3127, Data[,2])

tmp.biw <- wtc(d1, d2, nrands = 30)

plot.biwavelet(tmp.biw, plot.cb = TRUE)




# NOT RUN {
x<- 1:10
y<- 1:15
z<- outer( x,y,"+") 
image.plot(x,y,z) 

# or 
obj<- list( x=x,y=y,z=z)
image.plot(obj, legend.lab="Sverdrups")

# add some points on diagonal using standard plot function
#(with some clipping beyond 10 anticipated)
points( 5:12, 5:12, pch="X", cex=3)

# adding breaks and distinct colors for intervals of z
# with and without lab.breaks
brk<- quantile( c(z))
image.plot(x,y,z, breaks=brk, col=rainbow(4))
# annotate legend strip  at break values and add a label
image.plot(x,y,z, breaks=brk, col=rainbow(4),
           lab.breaks=names(brk))
#
# compare to 
zp <-quantile(c(z), c( .05, .1,.5, .9,.95))
image.plot(x,y,z, 
           axis.args=list( at=zp, labels=names(zp) ) )
# a log scaling for the colors
ticks<- c( 1, 2,4,8,16,32)
image.plot(x,y,log(z), axis.args=list( at=log(ticks), labels=ticks))

# see help file for designer.colors to generate a color scale that adapts to 
# quantiles of z. 
# Two add some color scales together here is an example of  5 blues to white to 5 reds
# with white being a specific size.
colorTable<- designer.colors(11, c( "blue","white", "red") )
# breaks with a gap of 10 to 17 assigned the white color
brks<- c(seq( 1, 10,,6), seq( 17, 25,,6)) 
image.plot( x,y,z,breaks=brks, col=colorTable)
#
#fat (5 characters wide) and short (50% of figure)  color bar on the bottom
image.plot( x,y,z,legend.width=5, legend.shrink=.5, horizontal=TRUE) 

# adding a label with all kinds of additional arguments.
# use side=4 for vertical legend and side= 1 for horizontal legend
# to be parallel to axes. See help(mtext).

image.plot(x,y,z, 
           legend.args=list( text="unknown units",
                             col="magenta", cex=1.5, side=4, line=2))

# and finally add some grid lines
dx <- x[2] - x[1]  
dy <- y[2] - y[1]  
xtemp<- seq(  min( x)- dx/2, max(x)+ dx/2,
              length.out = length(x) +1) 
ytemp<- seq(  min( y)- dy/2, max(y)+ dy/2,
              length.out = length(y) +1)
xline( xtemp, col="grey", lwd=2)
yline( ytemp, col="grey", lwd=2)
#

#### example using a irregular quadrilateral grid
data( RCMexample)

image.plot( RCMexample$x, RCMexample$y, RCMexample$z[,,1])
ind<- 50:75 # make a smaller image to show bordering lines
image.plot( RCMexample$x[ind,ind], RCMexample$y[ind,ind], RCMexample$z[ind,ind,1],
            border="grey50", lwd=2)


#### multiple images with a common legend

set.panel()

# Here is quick but quirky way to add a common legend to several plots. 
# The idea is leave some room in the margin and then over plot in this margin

par(oma=c( 0,0,0,4)) # margin of 4 spaces width at right hand side
set.panel( 2,2) # 2X2 matrix of plots

# now draw all your plots using usual image command
for (  k in 1:4){
  data<- matrix( rnorm(150), 10,15)
  image( data, zlim=c(-4,4), col=tim.colors())
  # and just for fun add a contour plot  
  contour( data, add=TRUE)
}

par(oma=c( 0,0,0,1))# reset margin to be much smaller.
image.plot( legend.only=TRUE, zlim=c(-4,4)) 

# image.plot tricked into  plotting in margin of old setting 

set.panel() # reset plotting device

#
# Here is a more learned strategy to add a common legend to a panel of
# plots  consult the split.screen help file for more explanations.
# For this example we draw two
# images top and bottom and add a single legend color bar on the right side 

# first divide screen into the figure region (left) and legend region (right)
split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))

# now subdivide up the figure region into two parts
split.screen(c(2,1), screen=1)-> ind
zr<- range( 2,35)
# first image
screen( ind[1])
image( x,y,z, col=tim.colors(), zlim=zr)

# second image
screen( ind[2])
image( x,y,z+10, col=tim.colors(), zlim =zr)

# move to skinny region on right and draw the legend strip 
screen( 2)
image.plot( zlim=zr,legend.only=TRUE, smallplot=c(.1,.2, .3,.7),
            col=tim.colors())

close.screen( all=TRUE)


# you can always add a legend arbitrarily to any plot;
# note that here the plot is too big for the vertical strip but the
# horizontal fits nicely.
plot( 1:10, 1:10)
image.plot( zlim=c(0,25), legend.only=TRUE)
image.plot( zlim=c(0,25), legend.only=TRUE, horizontal =TRUE)

# combining the  usual image function and adding a legend
# first change margin for some more room
# }
# NOT RUN {
par( mar=c(10,5,5,5))
image( x,y,z, col=topo.colors(64))
image.plot( zlim=c(0,25), nlevel=64,legend.only=TRUE, horizontal=TRUE,
            col=topo.colors(64))
# }
# NOT RUN {
#
# 
# sorting out the difference in formatting between matrix storage 
# and the image plot depiction 
# this really has not much to do with image.plot but I hope it is useful

A<- matrix( 1:48, ncol=6, nrow=8)
# first column of A will be 1:8 
#   ...  second  is  9:16 

image.plot(1:8, 1:6, A)
# add labels to each box 
text( c( row(A)), c( col(A)), A)
# and the indices ...
text( c( row(A)), c( col(A))-.25,  
      paste( "(", c(row(A)), ",",c(col(A)),")", sep=""), col="grey")

# "columns" of A are horizontal and rows are ordered from bottom to top!
#
# matrix in its usual tabular form where the rows are y  and columns are x
image.plot( t( A[8:1,]), axes=FALSE)

# }










