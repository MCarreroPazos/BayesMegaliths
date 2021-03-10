##############################
### FIGURE 6. SITE COUNTS ###
############################

# Load Required Libraries
library(rcarbon)
setwd("~/BayesMegaliths-master") # Set first the working directory at the BayesMegaliths-master folder, as this script is coded to have that folder as the root
library(here)

## A) CREATE THE SAMPLE POSTERIORS

# Read bespoke R functions (located in "src" folder)
source(here('src','oxcalScriptCreator.R'))

# Read C14 database with only those sites that contains more than 1 radiocarbon date
dates = read.csv2(here('data','C14dates_Iberia1.csv'), sep=";",na='n/a')

# Create OxCalScript (and save it into oxcalscripts folder)
oxcalScriptGen(id=dates$LabNumber,c14age=dates$C14,errors=dates$STD,group=NULL,site=dates$SiteID,fn=here('oxcalscripts','script.oxcal'),interval=100,mcnsim=5000,mcname="MCMC_uniform",model=c("uniform"))

# Run OxCalScript Locally (Notice this may take several hours).
quickSetupOxcal()
oxcalscript <- read_file(here("oxcalscripts","script.oxcal"))
result_file <- executeOxcalScript(oxcalscript) # This will create an output (MCMC_uniform.csv) stored in a local temporary file (in Windows OS). Please, manually search for it and save it on "oxcalresults" folder

# Create the MCMC samples posteriors
mcmcoutput=read.csv(here("oxcalresults", "MCMC_uniform.csv"))
mcmcoutput = mcmcoutput[,-c(1,ncol(mcmcoutput))] # Remove unused columns
index=which(apply(mcmcoutput,1,max)>1950) # identify instances with posterior outside calibration range
mcmcoutput = mcmcoutput[-index,]
mcmcoutput = apply(mcmcoutput,2,BCADtoBP) # Convert to BP

# Each row of mcmcoutput is a possible start and end date for each site
# You can run this as an example:
rindex = sample(nrow(mcmcoutput),size=1)
plot(0,0,type='l',xlim=rev(range(mcmcoutput[rindex,])),ylim=c(0.5,3.5),axes=FALSE,xlab='cal BP',ylab='Sites')
axis(2,at=1:3,labels=unique(test$SiteIDNew))
axis(1)
lines(c(mcmcoutput[rindex,1:2]),c(1,1),lwd=2)
lines(c(mcmcoutput[rindex,3:4]),c(2,2),lwd=2)
lines(c(mcmcoutput[rindex,5:6]),c(3,3),lwd=2)
# On this simple plot you can see how we accounted for site duration in the case of three sites

## B) SITE COUNTS

# Read bespoke R functions (located in "src" folder)
source(here('src','sampleCounts.R')) 

# Read MCMC output
x = read.csv(here('oxcalresults','MCMC_uniform.csv'))

# Eliminate pass number and last column
x = x[,-c(1,ncol(x))]

# Count the number of sites for a given chronological block, for each of the possible results
cntMatrix = sampleCounts(x) # This may take a couple of minutes

# Plot Figure 6 and save it to "figures" folder
png(file=here('figures','Figure 6.png'),width = 480,height=480)
plotCntMat(cntMatrix)
dev.off()
