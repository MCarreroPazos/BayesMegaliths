# Load Required Libraries
library(readxl)
library(oxcAAR)
library(rcarbon)
library(readr)
# Read 14C Dates
dates = read.csv("~data\\C14dates_Iberia.csv", sep=";")

# Read bespoke R functions
source('~/src/oxcalScriptCreator.R')

# Create OxCalScript (and save it into oxcalscripts file)
oxcalScriptGen(id=dates$LabNumber,c14age=dates$C14,errors=dates$STD,group=NULL,site=dates$SiteID,fn="~oxcalscripts\\script.oxcal",interval=100,mcnsim=5000,mcname="MCMC_uniform",model=c("uniform"))

# Run OxCalScript Locally (Notice this requires about 30-50 hours of processing)
#This will create an output (MCMC_uniform.csv) stored in a local temporary file (in Windows)
quickSetupOxcal()
oxcalscript <- read_file('~oxcalscripts\\script.oxcal')
result_file <- executeOxcalScript(oxcalscript)

#Create the MCMC samples posteriors
mcmcoutput=read.csv("~oxcalresults\\MCMC_uniform.csv")
mcmcoutput = mcmcoutput[,-c(1,ncol(mcmcoutput))] #Remove unused columns
index=which(apply(mcmcoutput,1,max)>1950) #identify instances with posterior outside calibration range
mcmcoutput = mcmcoutput[-index,]
mcmcoutput = apply(mcmcoutput,2,BCADtoBP) # Convert to BP

# Each row of mcmcoutput is a possible start and end date for each site. If you try running the script below you can see what I mean. 
rindex = sample(nrow(mcmcoutput),size=1)
plot(0,0,type='l',xlim=rev(range(mcmcoutput[rindex,])),ylim=c(0.5,3.5),axes=FALSE,xlab='cal BP',ylab='Sites')
axis(2,at=1:3,labels=unique(test$SiteIDNew))
axis(1)
lines(c(mcmcoutput[rindex,1:2]),c(1,1),lwd=2)
lines(c(mcmcoutput[rindex,3:4]),c(2,2),lwd=2)
lines(c(mcmcoutput[rindex,5:6]),c(3,3),lwd=2)

# The idea is to count the number of sites for a given chronological block (say of 100 years), for each of these possible results. In this particular example we might have a maximum of only one site per time-block. Like in this particular case:

plot(0,0,type='l',xlim=rev(range(mcmcoutput[1,])),ylim=c(0.5,3.5),axes=FALSE,xlab='cal BP',ylab='Sites')
axis(2,at=1:3,labels=unique(test$SiteIDNew))
axis(1)
lines(c(mcmcoutput[1,1:2]),c(1,1),lwd=2)
lines(c(mcmcoutput[1,3:4]),c(2,2),lwd=2)
lines(c(mcmcoutput[1,5:6]),c(3,3),lwd=2)

# or we might have an instance of coexistence like in the example below wherer S2 and S3 coexists towards the end of the 4th millennium.  

plot(0,0,type='l',xlim=rev(range(mcmcoutput[1049,])),ylim=c(0.5,3.5),axes=FALSE,xlab='cal BP',ylab='Sites')
axis(2,at=1:3,labels=unique(test$SiteIDNew))
axis(1)
lines(c(mcmcoutput[1049,1:2]),c(1,1),lwd=2)
lines(c(mcmcoutput[1049,3:4]),c(2,2),lwd=2)
lines(c(mcmcoutput[1049,5:6]),c(3,3),lwd=2)

# the next step is to run this for the entire data so that we can have a probabilistic estimate in the number of used sites for any time-interval.