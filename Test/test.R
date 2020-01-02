## This file is used to test the counterparts of
## analysis.wavelet, reconstruct and wt.image (from R package "WaveletComp")
## in R package "CoFESWave"

## The names are temporary.
## I now use CoFESWave.whatwhat, so when the names are decided,
## it would be easy to replace them with new names.


library(WaveletComp)
library(tictoc)

#### data ####
load("H:/R/R_files/RA_wavelet/CoFESWave/Test/ExampleData.Rdata")
x = Data[,1]
Data_raw <- data.frame(date = date_Sector, x = x)

#### analysis.wavelet ####
tic()
for (i in c(1:10)) {
  my.w_original <- WaveletComp::analyze.wavelet(Data_raw, "x",
                                              loess.span = 0.0, make.pval = F, verbose = F ,
                                              n.sim = 10000, upperPeriod = length(x))
}
toc()

tic()
for (i in c(1:10)) {
  my.w_rewrite <- CoFESWave.Transform(Data_raw, "x",
                                      loess.span = 0.0, make.pval = F, verbose = F ,
                                      n.sim = 10000, upperPeriod = length(x))
}
toc()

#### reconstruct ####
recon_org0 <- WaveletComp::reconstruct(my.w_original, plot.waves = FALSE, lwd = c(1,2),
                                       legend.coords = "topright",plot.rec = FALSE, verbose = F)

recon_org1 <- CoFESWave.reconstruct(my.w_rewrite, plot.waves = FALSE, lwd = c(1,2),
                                    legend.coords = "topright",plot.rec = FALSE, verbose = F)

#### wt.image ####
library(latex2exp)
topl_colors <- colorRampPalette(c("#64d8cb","#26a69a", "#90a4ae","#5f5fc4", "#283593"))

WaveletComp::wt.image(my.w_original, periodlab = " ",timelab = "  " , main = " ",
                      legend.params = list(lab = "wavelet power levels", mar = 5.1, cex = 6, n.ticks = 10),
                      color.key = "quantile", lwd = 2, plot.ridge = FALSE, color.palette = "topl_colors(n.levels)",
                      show.date = TRUE)

CoFESWave.image(my.w_rewrite, periodlab = " ",timelab = "  " , main = " ",
                legend.params = list(lab = "wavelet power levels", mar = 5.1, cex = 6, n.ticks = 10),
                color.key = "quantile", lwd = 2, plot.ridge = FALSE, color.palette = "topl_colors(n.levels)",
                show.date = TRUE)






