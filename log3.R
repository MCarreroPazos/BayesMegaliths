library(rcarbon)
library(doParallel)
library(here)

# Load data
dates = read.csv(here('data','C14dates_Iberia.csv'), sep=";",na='n/a')

# Calibrate dates
dates_cal<- calibrate(x=dates$C14, 
                      errors=dates$STD, normalised=FALSE, 
                      calMatrix=TRUE)

# Create the bin-width (here h = 100)
bin_width <- binPrep(sites=as.character(dates$Site),
                      ages=dates$C14, h=100)

# Setup SPD limits, computer cores and number of simulations
start_sumprob <- 7000
end_sumprob <- 2000

cores <- detectCores()
cores_utils <- 6  

nsim=5000

#Iberia SPD (Null model = exponential). Notice that this takes a few hours of processing
Iberia_exp_h100 <- modelTest(dates_cal, errors=dates$STD, bins=bin_width, 
                             nsim=nsim, runm=100,timeRange=c(start_sumprob,end_sumprob),raw = TRUE,
                             model="exponential",datenormalised=FALSE,ncores=cores_utils)

saveRDS(Iberia_exp_h100, file = here("data","Iberia_SPD.rds"))
#To recall Iberia SPDs: Iberia_exp_h100 <-readRDS(here("Iberia_SPD.rds",refhook = NULL))

#Figure SPD
plot(Iberia_exp_h100, xlim=c(7000, 2000))
pan <-layout(matrix(c(1,1,2,2), ncol=2, byrow=TRUE), heights=c(4, 1))
par(mai=rep(0.64, 4))
layout.show(pan)
plot(Iberia_exp_h100, xlim=c(7000, 2000))
lines(Iberia_exp_h100$fit$calBP,Iberia_exp_h100$result$PrDens, type = "l", col="blue")
lines(Iberia_exp_h100$fit$calBP,Iberia_exp_h100$fit$PrDens,type="l",lty=2,col="red")
# insert plot title
title("Iberia SPD based on exponentinal null model")
mtext("null model: p.value < 0.001", side = 3, line = 0.25, cex = 0.80)
#insert legend
par(mai=c(0,0,0,0))
plot.new()
# panel.first = grid()
legend(x="center", ncol=5,legend=c("SPD","Conf. Interval", "Null model", "positive dev.","negative dev."),
       col =c("blue","lightgrey","red","indianred","royalblue"),lty=c(1,1,2,1,1),
       lwd=c(1,5,1,5,5),cex=0.72,xjust = -2)