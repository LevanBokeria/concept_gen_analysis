# This script will plot the power by various max N per group

# Clear the environment
rm(list=ls())

# Libraries
library(ggplot2)
library(gapminder)
library(plyr)
library(tidyverse)

# Define global variables

nIterEv <- 10000
d1      <- 0.5
d1_str  <- '05'
nLimit  <- 200
crit1   <- 6
crit2   <- 1/6

nFrom <- 80
nTo   <- 128
nBy   <- 8

altNs   <- seq(nFrom,nTo,by = nBy)

saveFig <- T

# Construct a function to get stats from simulations 
loadFile  <- paste('resultsByManyNs_d_',
                   d1_str,
                   '_crit1_', crit1, '_',
                   altNs[1], '_to_', altNs[length(altNs)],
                   '_by_', nBy,
                   '.RData',sep='')

# load data
load(loadFile)

outStats1 <- outData1$outStats
outStats2 <- outData2$outStats

plotStats <- function(dataPlotted){
  
  d_plotted <- dataPlotted$d[1]
  
  if (saveFig){
    png(file=paste('./power_by_n_plots/power_by_n_d_', d_plotted,
                   '_BF10_', crit1, '_BF01_', 1/crit2,'_',
                   altNs[1], '_to_', altNs[length(altNs)],
                   '_by_', nBy,
                   '.png',
                   sep=''))
  }
  
  matplot(dataPlotted[,3:5], type = c("b"), pch=1,col = 3:5, lwd = 2,
          xaxt="n",yaxt="n",
          main=paste('d=', d_plotted, ' BF10=', crit1, ' BF01=', 1/crit2, sep=''),
          xlab='N per group', ylab='% simulations')
  legend("topleft", legend = c('H0','H1','Und'), col=3:5, pch=1,cex = 0.7) # optional legend
  
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
       lwd = par("lwd"), equilogs = TRUE)
  
  axis(side=1,at=1:nrow(dataPlotted),labels=paste(dataPlotted$maxN))
  axis(side=2,at=seq(0,10000,by=1000),labels=paste(seq(0,100,by=10)),cex.axis=1)
  
  if (saveFig){
    dev.off()
  }
  
}

plotStats(outStats1)
plotStats(outStats2)