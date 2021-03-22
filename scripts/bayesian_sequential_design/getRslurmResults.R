# This is a script reading the output from the rslurm for the talk Bayesian sequential designs Wed 29 July 2020

# clear things
rm(list=ls())

# Define parameters
d1     <- 0.5
d1_str <- '05' # because this is how slurm names the folders
crit1  <- 10
nLimit <- 135
nIter  <- 10000

# Where are the simulation results
loadFolder <- paste('_rslurm_d1_',d1_str,'_limpg_',
                    nLimit, '_crit1_', crit1, sep='')

# Where should we save the postprocessed dataframe df3?
saveName   <- paste(loadFolder, '/sanity_simulationResults_',d1_str,'_limpg_',
                    nLimit, '_crit1_', crit1, '.RData',sep='')

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
  n <- c(n, seq(10, n_raw[i], batchSize[1]))
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
save(list = c('df3','nIter'), file = saveName)


