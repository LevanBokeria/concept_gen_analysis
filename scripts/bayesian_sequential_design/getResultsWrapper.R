# This script is a wrapper to get bunch of simulation results and stats quickly 
# in a table.

# It will take the rslurm output, turn it into a nice dataframe (called df3), 
# and then analyze probabilities of supporting H1 and H0.

# Clear the environment
rm(list=ls())

# Libraries
library(plyr)
library(BayesFactor)
library(assortedRFunctions)
library(tidyverse)

# Define global variables
nIter   <- 10000
d1      <- 0.5
d1_str  <- '05'
nLimit  <- 200
crit1   <- 6
crit2   <- 1/6
minN    <- 24
batchSize <- 8

saveDF <- TRUE # should df3 be saved? If already done, assign FALSE
readDF <- FALSE # if DF was saved previously, just read it without recreating it

saveOutData <- TRUE # save the resulting table containing probabilities of H1 and H0?

# Function to construct dataframe from slurm results
getRslurmResults <- function (nIter,d1_str,nLimit,crit1,saveDF,readDF){
  
  loadFolder <- paste('_rslurm_d1_',d1_str,'_limpg_',
                      nLimit, '_crit1_', crit1, '_minN_', minN,
                      '_batchSize_', batchSize, sep='')
  saveName   <- paste(loadFolder, '/simulationResults_',d1_str,'_limpg_',
                      nLimit, '_crit1_', crit1, '_minN_', minN,
                      '_batchSize_', batchSize, '.RData',sep='')
  if(readDF){
    
    df3 <- load(saveName)
    
  } else {

    # Load data
    paths      <- loadFolder
    tempList   <- readRDS(paste0(paths, '/results_0.RDS'), refhook = NULL)
    tempUnList <- unlist(tempList)
    
    tempUnllist_names_org <- names(tempUnList)
    tempUnllist_names     <- tempUnllist_names_org
    tempUnllist_names[tempUnllist_names_org == 'bf'] <- 'bf0' # replace bf with bf0 otherwise it's get removed.
    tempUnllist_names <- tempUnllist_names[!(tempUnllist_names == 'crit1' | tempUnllist_names == 'crit2')] # Remove crit1 and crit2
    tempUnllist_names <- gsub("[^0-9.-]", "", tempUnllist_names) # Remove all non-numerical info
    
    tempUnllist_seq   <- as.numeric(tempUnllist_names[tempUnllist_names != ''])
    tempUnllist_seq[tempUnllist_seq == 0] <- 1 # Replace zero with 1 again so that the id will repeated once
    tempUnllist_seq   <- tempUnllist_seq - c(tempUnllist_seq[2:length(tempUnllist_seq)], 1)
    tempUnllist_seq   <- tempUnllist_seq[tempUnllist_seq != -1] + 1
    
    # Get columns for df
    d         <- tempUnList[tempUnllist_names_org == 'd']
    n_raw     <- tempUnList[tempUnllist_names_org == 'n']
    crit1     <- tempUnList[tempUnllist_names_org == 'crit1']
    crit2     <- tempUnList[tempUnllist_names_org == 'crit2']
    batchSize <- tempUnList[tempUnllist_names_org == 'batchSize']
    limit     <- tempUnList[tempUnllist_names_org == 'limit']
    
    # Repeat by how many steps were taken
    id        <- rep(1:length(tempList), times = tempUnllist_seq)
    d         <- rep(d, times = tempUnllist_seq)
    crit1     <- rep(crit1, times = tempUnllist_seq)
    crit2     <- rep(crit2, times = tempUnllist_seq)
    batchSize <- rep(batchSize, times = tempUnllist_seq)
    limit     <- rep(limit, times = tempUnllist_seq)
    
    # Create sample size sequences
    n <- c()
    for(i in 1:length(n_raw)){
      n <- c(n, seq(minN, n_raw[i], batchSize[1]))
    }
    
    # Make DF
    df3 <- data.frame(id = id,
                      d = d,
                      n = n,
                      crit1 = crit1,
                      crit2 = crit2,
                      batchSize = batchSize,
                      limit = limit, 
                      bf = tempUnList[grep("bf", tempUnllist_names_org)])
  
    # Save dfs
    if (saveDF){
      save(list = c('df3','nIter'), file = saveName)  
    }
  
  } #if loadDF
  
  return(df3)
  
}

# Function to get stats from simulations 
getStats <- function(crit1,crit2,nIter,d1_str,d1,nLimit){
  loadFile  <- paste('_rslurm_d1_',
                     d1_str,'_limpg_',nLimit, '_crit1_', crit1, '_minN_', minN,
                     '_batchSize_', batchSize,
                     '/simulationResults_',d1_str,'_limpg_',nLimit, '_crit1_', crit1, 
                     '_minN_', minN, '_batchSize_', batchSize,
                     '.RData',sep='')
  
  # load data
  load(loadFile)
  
  # Process for d = 0
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
  
 
  # Process for d = d1
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
  
  # Create a table to record all the data in one place
  outData <- c(length(d0_H0$id),length(d0_H1$id),length(d0_undecided$id),
               length(d1_H0$id),length(d1_H1$id),length(d1_undecided$id))
  
  names(outData) <- c('d0_H0','d0_H1','d0_undecided','d1_H0','d1_H1','d1_undecided')
  
  outData <- as_tibble(as.list(outData))
  
  return(outData)
}

# Call the function to construct a dataframe from slurm results
df3 <- getRslurmResults(nIter,d1_str,nLimit,crit1,saveDF,readDF)

# Call the function to get the probabilities data
outData <- getStats(crit1,crit2,nIter,d1_str,d1,nLimit)

# Save outData
loadFolder <- paste('_rslurm_d1_',d1_str,'_limpg_',
                    nLimit, '_crit1_', crit1, '_minN_', minN,
                    '_batchSize_', batchSize, sep='')
saveNameOutData   <- paste(loadFolder, '/outData.RData',sep='')

if (saveOutData){
  save(outData, file = saveNameOutData)
}


