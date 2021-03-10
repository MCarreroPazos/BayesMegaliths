################################################################################################
### FIGURE 9. OBSERVED TIME SERIES PER TYPE OF MONUMENT (ARQUITECTURE) AND PERMUTATION TEST ###
##############################################################################################

# Load Required Libraries
library(rcarbon)
library(dplyr)
setwd("~/BayesMegaliths-master") # Set first the working directory at the BayesMegaliths-master folder, as this script is coded to have that folder as the root
library(here)

# Read C14 database with only those sites which have more than 1 radiocarbon date
dates = read.csv2(here('data','C14dates_Iberia1.csv'), sep=";",na='n/a')

# Read bespoke R functions (located in "src" folder), and MCMC output (located in "Oxcalresults" folder)
source(here('src','markTest.R'))
source(here('src','sampleCounts.R'))
x = read.csv(here('oxcalresults','MCMC_uniform.csv'))

# Eliminate pass number and last column
x = x[,-c(1,ncol(x))]

# Group the dates by architecture
sites_type = unique(select(dates,SiteID,Architecture))

# Reorder sites_type to match MCMC output
mcmc_sitenames = unlist(lapply(strsplit(colnames(x),'\\.'),function(x){x[[2]]}))
nsites = length(mcmc_sitenames)/2
sites_type=sites_type[match(mcmc_sitenames[seq(1,by=2,length.out = nsites)],sites_type$SiteID),]
marks = sites_type$Architecture

# Compute observed time series by region and run permutation test with 50 random simulations
res1=markTest(x=x,marks=marks,nsim=50) # Note this takes several hours of processing with just 50 simulations. Increasing the number will raise the computational cost enormously.

# Plot Figure 9 and save it into "figures" folder
png(file=here('figures','Figure 9.png'),width = 1000,height=600)
par(mar=c(4,2,4,2),mfrow=c(2,4))
for (i in 1:length(unique(marks)))
{
  plotMarkTest(res1,index=i,main=unique(marks)[i])
  legend('topright',legend=c(paste0('n=',sum(marks==unique(marks)[i])),paste0('P-Value=',round(res1$pValueList[i],4))),bty='n',cex=1)
}
dev.off()