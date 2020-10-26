library(rcarbon)
library(here)
library(dplyr)
source(here('src','markTest.R'))

# Read raw data
dates = read.csv(here('data','C14dates_Iberia.csv'), sep=";")
sites_region = unique(select(dates,SiteID,Region))

# Read MCMC output
x = read.csv(here('oxcalresults','MCMC_uniform.csv'))
x = x[,-c(1,ncol(x))] # eliminate pass number and last column
x=x[4900:5000,] #run this only on 100 posterior samples


# Reorder sites_region to match MCMC output
mcmc_sitenames = unlist(lapply(strsplit(colnames(x),'\\.'),function(x){x[[2]]}))
nsites = length(mcmc_sitenames)/2
sites_region=sites_region[match(mcmc_sitenames[seq(1,by=2,length.out = nsites)],sites_region$SiteID),]
marks = sites_region$Region

res=markTest(x=x,marks=marks,nsim=50)


png(file=here('figures','testPlotMarkTest.png'),width = 720,height=480)
par(mar=c(5,4,2,1),mfrow=c(2,3))
for (i in 1:length(unique(marks)))
{
plotMarkTest(res,index=i,main=unique(marks)[i])
legend('topleft',legend=c(paste0('n=',sum(marks==unique(marks)[i])),paste0('P-Value=',round(res$pValueList[i],4))),bty='n',cex=1)
}


dev.off()
