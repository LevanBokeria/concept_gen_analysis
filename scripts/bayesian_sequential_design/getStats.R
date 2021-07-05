# Description

# Function to calculate how many simulations supported H1, H0 or were undecided.

# Clear the environment
rm(list=ls())

# Libraries
library(ggplot2)
library(gganimate)
library(gapminder)
library(gridExtra)
library(grid)
library(plyr)
library(BayesFactor)
library(assortedRFunctions)

# Standard values
crit1        <- 10
crit2        <- 1/6

nIter        <- 10000


# Global settings
fontSize <- 15
theme_set(theme_gray(base_size = fontSize))


d1_str <- "02"
d1     <- 0.2

nLimit <- 150
nLimSimulated <- 150

# r_scale <- 'wide'
# 
# loadFile  <- paste('_rslurm_SeqDesLimit_d1_',
#                    d1_str,'_limpg_',nLimit,'_rscale_',r_scale, '_crit1_', crit1,
#                    '/simulationResults_',d1_str,'_limpg_',nLimit, '_crit1_', crit1,
#                    '.RData',sep='')
# saveName1 <- paste('_rslurm_SeqDesLimit_d1_',
#                    d1_str,'_limpg_',nLimit,'_rscale_',r_scale, '_crit1_', crit1,
#                     '/figures/limit1.png',sep='')
# saveName2 <- paste('_rslurm_SeqDesLimit_d1_',
#                    d1_str,'_limpg_',nLimit,'_rscale_',r_scale, '_crit1_', crit1,
#                    '/figures/limit2.png',sep='')

loadFile  <- paste('_rslurm_d1_',
                   d1_str,'_limpg_',nLimit, '_crit1_', crit1,
                   '/simulationResults_',d1_str,'_limpg_',nLimit, '_crit1_', crit1,
                   '.RData',sep='')

# load data
load(loadFile)

# /* 
# ----------------------------- Sequential design with limit 1---------------------------
# */
# For plotting
minN         <- 10
maxN         <- nLimSimulated

tempDF <- subset(df3, d == 0.0)

tempDF$trans_bf <- NA
tempDF$trans_bf[tempDF$bf < 1] <- -1/tempDF$bf[tempDF$bf < 1] + 1
tempDF$trans_bf[tempDF$bf > 1] <- tempDF$bf[tempDF$bf > 1] - 1

tempDF_agg <- ddply(tempDF, c('id'), summarise, n = n[length(n)], bf = bf[length(bf)])
tempDF_agg$support <- 'undecided'
tempDF_agg$support[tempDF_agg$bf > crit1] <- 'H1'
tempDF_agg$support[tempDF_agg$bf < crit2] <- 'H0'

d0_undecided <- tempDF_agg[tempDF_agg$support == 'undecided',]
d0_H1 <- tempDF_agg[tempDF_agg$support == 'H1',]
d0_H0 <- tempDF_agg[tempDF_agg$support == 'H0',]

# /*
# ----------------------------- Sequential design with limit 2---------------------------
# */
# For plotting
minN         <- 10
bar_yAxis    <- 0.22
hist_yAxis   <- bar_yAxis*nIter

tempDF <- subset(df3, d == d1)

tempDF$trans_bf <- NA
tempDF$trans_bf[tempDF$bf < 1] <- -1/tempDF$bf[tempDF$bf < 1] + 1
tempDF$trans_bf[tempDF$bf > 1] <- tempDF$bf[tempDF$bf > 1] - 1

tempDF_agg <- ddply(tempDF, c('id'), summarise, n = n[length(n)], bf = bf[length(bf)])
tempDF_agg$support <- 'undecided'
tempDF_agg$support[tempDF_agg$bf > crit1] <- 'H1'
tempDF_agg$support[tempDF_agg$bf < crit2] <- 'H0'

d1_undecided <- tempDF_agg[tempDF_agg$support == 'undecided',]
d1_H1        <- tempDF_agg[tempDF_agg$support == 'H1',]
d1_H0        <- tempDF_agg[tempDF_agg$support == 'H0',]

