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

# Reorder sites_region to match MCMC output
mcmc_sitenames = unlist(lapply(strsplit(colnames(x),'\\.'),function(x){x[[2]]}))
nsites = length(mcmc_sitenames)/2
sites_region=sites_region[match(mcmc_sitenames[seq(1,by=2,length.out = nsites)],sites_region$SiteID),]
marks = sites_region$Region

res=markTest(x=x,marks=marks,nsim=50)

png(file=here('figures','PlotMark.png'),width = 1000,height=600)
par(mar=c(5,4,2,1),mfrow=c(2,3))
for (i in 1:length(unique(marks)))
{
plotMarkTest(res,index=i,main=unique(marks)[i])
legend('topleft',legend=c(paste0('n=',sum(marks==unique(marks)[i])),paste0('P-Value=',round(res$pValueList[i],4))),bty='n',cex=1)
}
par(mai=c(0,0,0,0))
plot.new()
legend(x="center",legend=c("Observed SPD","Null hypothesis (pan-regional trend)", "Positive deviation","Negative deviation"),
       col =c("black","lightgrey","indianred","royalblue"),lty=c(1,1,1,1),
       lwd=c(1,5,5,5),cex=2,xjust = 0.5)
dev.off()
