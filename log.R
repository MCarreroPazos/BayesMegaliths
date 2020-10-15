# Load Required Libraries
library(readxl)
library(oxcAAR)
library(rcarbon)
library(readr)
# Read 14C Dates
# I noticed a couple of issues here:
# *SiteID does not seem to be actually a SiteID but an unique identified for each date.
# *The fields C14 and STD are read as text, and not as number. Theis seems to be because there are some hidden characters in various places. This should be cleaned, and stored in a proper format (i.e. CSV, not xlsx).
dates = as.data.frame(read_xlsx('../C14dates_Iberia.xlsx'))

# Read bespoke R functions
source('../src/oxcalScriptCreator.R')

# Fix Site ID
dates$SiteIDNew = paste0('S',as.numeric(as.factor(dates$Site)))

# Assign U

# Test Run with just three sites
test = subset(dates, SiteIDNew %in% c('S1','S2','S3'))
test$C14 = as.numeric(test$C14)
test$STD = as.numeric(test$STD)

# Create OxCalScript. The script below generates an oxcal script file (samplescript.oxcal) that can be copied on the browser for an online execution, or alternatively run locally.
# Notice that:
# The command can be run for each site (faster) or all sites together. The latter might be convenient but perhaps the former is easier for troubleshooting and identifyng instances of low agreement
# The argument group should be dedicated for instances where radiocarbon dates are from the same object, and hence when we expect the exact same date. Did not found anything relevant so left it blank.
oxcalScriptGen(id=test$LabNumber, c14age=test$C14,errors=test$STD,group=NULL,site=test$SiteIDNew,fn='../oxcalscripts/samplescript.oxcal',interval=100,mcnsim=5000,mcname='testMCMC',model=c("uniform"))

# Run OxCalScript Locally (I need to upgrade my linux install but should work elsewhere)
# quickSetupOxcal()
# oxcalscript <- read_file('../oxcalscripts/samplescript.oxcal')
# result_file <- executeOxcalScript(oxcalscript)

#Whether locally or remotely the script produces a file called testMCMC.csv which contains the sample posteriors. This is what we actually need for the nextrstage of analyses. I stored the one I just executed in the folder oxcalresults

mcmcoutput=read.csv('../oxcalresults/testMCMC.csv')
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







