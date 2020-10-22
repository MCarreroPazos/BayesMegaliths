library(rcarbon)
library(here)
source(here('src','sampleCounts.R'))
x = read.csv(here('oxcalresults','MCMC_uniform.csv'))
x = x[,-c(1,ncol(x))] # eliminate pass number and last column
cntMatrix = sampleCounts(x) #this will take a couple of minutes
png(file=here('figures','testPlot.png'),width = 480,height=480)
plotCntMat(cntMatrix)
dev.off()



