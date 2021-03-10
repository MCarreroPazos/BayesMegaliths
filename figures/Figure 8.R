#######################################################################
### FIGURE 8. OBSERVED TIME SERIES PER REGION AND PERMUTATION TEST ###
#####################################################################

# Load Required Libraries
library(readxl)
library(maptools)
library(spatstat)
library(rgdal)
library(rcarbon)
library(dplyr)
library(readr)
setwd("~/BayesMegaliths-master") # Set first the working directory at the BayesMegaliths-master folder, as this script is coded to have that folder as the root
library(here)

# Read C14 database with only those sites that contains more than 1 radiocarbon date
dates = read.csv(here('data','C14dates_Iberia1.csv'), sep=";",na='n/a')

# Read bespoke R functions (located in "src" folder), and MCMC output (located in "Oxcalresults" folder)
source(here('src','markTest.R'))
source(here('src','sampleCounts.R'))
x = read.csv(here('oxcalresults','MCMC_uniform.csv'))

# Eliminate pass number and last column
x = x[,-c(1,ncol(x))]

# Group the dates by region
sites_region = unique(select(dates,SiteID,Region))

# Reorder sites_region to match MCMC output
mcmc_sitenames = unlist(lapply(strsplit(colnames(x),'\\.'),function(x){x[[2]]}))
nsites = length(mcmc_sitenames)/2
sites_region=sites_region[match(mcmc_sitenames[seq(1,by=2,length.out = nsites)],sites_region$SiteID),]
marks = sites_region$Region

# Compute observed time series by region and run permutation test with 50 random simulations
res=markTest(x=x,marks=marks,nsim=50) # Note this takes several hours of processing with just 50 simulations. Increasing the number will raise the computational cost enormously.

# Plot Figure 8 and save it into "figures" folder
png(file=here('figures','Figure 8.png'),width = 1000,height=600)
par(mar=c(5,4,2,1),mfrow=c(3,2))
for (i in 1:length(unique(marks)))
{
  plotMarkTest(res,index=i,main=unique(marks)[i])
  legend('topleft',legend=c(paste0('n=',sum(marks==unique(marks)[i])),paste0('P-Value=',round(res$pValueList[i],4))),bty='n',cex=1)
}
dev.off()