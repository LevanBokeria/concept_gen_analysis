# A script to just quickly see probabilities of supporting H1 and H2 dependent 
# on sample size, effect size, etc.

library(BayesFactor)
library(assortedRFunctions)
nSim <- 10000
d1 <- 0.5
nMax <- 180
bf <- c()
for(i in 1:nSim){
    
    if (i %% 500 == 0){
        print(i)
    }
    
    bf[i] <- reportBF(ttestBF(rnorm(nMax), rnorm(nMax, d1), paired = FALSE))
}
bf10 <- mean(bf > 10)
bf01 <- mean(bf < 1/6)

save(bf,bf10,bf01,file=paste('./quickSimResults/d_',d1,'_nMax_',nMax,'.RData',sep=''))