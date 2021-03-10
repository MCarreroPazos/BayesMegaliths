# Load Required Libraries
library(readxl)
library(maptools)
library(spatstat)
library(rgdal)
library(oxcAAR)
library(rcarbon)
library(dplyr)
library(animation)
library(doParallel)
library(readr)
library(here)

### Initial comments ###
  ## See "data" folder
    # C14dates_Iberia.csv refers to the whole database (962 C14 dates belonging to 319 archaeological sites).
    # C14dates_Iberia1.csv refers to only those cases where we have more than 1 radiocarbon date (823 C14 dates, belonging to 178 archaeological sites)

##########################################################
### PART 1. COMPOSITE KERNEL DENSITY ESTIMATES (CKDE) ### 
###                AND STACKED SPDs                  ###
#######################################################

## COMPOSITE KERNEL DENSITY ESTIMATES (CKDE)
# Read raw data
dates = read.csv2(here('data','C14dates_Iberia.csv'), sep=";",na='n/a')

### Subset the data by region and calibrate with corresponding curve ###

## EAST ##
# Subset data
East <- subset(dates,Region=="East")
dates_terrestrial_East <- subset(East, Curve=="intcal20")
dates_marine_East <- subset(East, Curve=="marine20")
# Calibrate dates with the corresponding curve
East_Terr_cal <- calibrate(x=dates_terrestrial_East$C14,errors=dates_terrestrial_East$STD,normalised=FALSE, 
                           calMatrix=FALSE, calCurves = dates_terrestrial_East$Curve)

East_Marine_cal <- calibrate(x=dates_marine_East$C14, errors=dates_marine_East$STD,
                             normalised=FALSE, calMatrix=FALSE, calCurves = dates_marine_East$Curve, 
                             resOffsets = dates_marine_East$DeltaR, resErrors = dates_marine_East$DeltaRErr)

dates_cal_East = rcarbon::combine(East_Terr_cal,East_Marine_cal, fixIDs = TRUE)

## NORTH ##
# Subset data
North <- subset(dates,Region=="North")
# Calibrate dates with the corresponding curve
dates_cal_North <- calibrate(x=North$C14,errors=North$STD,normalised=FALSE, 
                             calMatrix=FALSE, calCurves = North$Curve)

## SOUTH ##
# Subset data
South <- subset(dates,Region=="South")
South_terrestrial_East <- subset(South, Curve=="intcal20")
South_marine_East <- subset(South, Curve=="marine20")
South_mixed_East <- subset(South, Curve=="mixed")
# Calibrate dates with the corresponding curve
South_Terr_cal <- calibrate(x=South_terrestrial_East$C14,errors=South_terrestrial_East$STD,normalised=FALSE, 
                            calMatrix=FALSE, calCurves = South_terrestrial_East$Curve)

South_Marine_cal <- calibrate(x=South_marine_East$C14, errors=South_marine_East$STD,
                              normalised=FALSE, calMatrix=FALSE, calCurves = South_marine_East$Curve, 
                              resOffsets = South_marine_East$DeltaR, resErrors = South_marine_East$DeltaRErr)

mixedCurve <- mixCurves('intcal20','marine20',p=0.8) # Set Mixed curve

South_Mixed_cal <- calibrate(x=South_mixed_East$C14,
                             errors=South_mixed_East$STD,
                             normalised=FALSE, calMatrix=FALSE, calCurves = mixedCurve, 
                             resOffsets = South_mixed_East$DeltaR, 
                             resErrors = South_mixed_East$DeltaRErr)
# Combine in one object
dates_cal_South = rcarbon::combine(South_Terr_cal, South_Marine_cal, South_Mixed_cal, fixIDs = TRUE)

## EBRO BASIN ##
# Subset data
Ebro_basin <- subset(dates,Region=="Ebro basin")
# Calibrate dates with the corresponding curve
dates_cal_Ebro <- calibrate(x=Ebro_basin$C14,errors=Ebro_basin$STD,normalised=FALSE, 
                            calMatrix=FALSE, calCurves = Ebro_basin$Curve)

## TAGUS BASIN ##
# Subset data
Tagus_basins <- subset(dates,Region=="Tagus and western basins")
# Calibrate dates with the corresponding curve
dates_cal_Tagus <- calibrate(x=Tagus_basins$C14,errors=Tagus_basins$STD,normalised=FALSE, 
                             calMatrix=FALSE, calCurves = Tagus_basins$Curve)

## DUERO BASIN ##
# Subset data
Duero_basin <- subset(dates,Region=="Douro/Duero and Mondego basins")
# Calibrate dates with the corresponding curve
dates_cal_Duero <- calibrate(x=Duero_basin$C14,errors=Duero_basin$STD,normalised=FALSE, 
                             calMatrix=FALSE, calCurves = Duero_basin$Curve)

### Composite Kernel Density Estimates (CKDE)
## Create bins
bins_East <- binPrep(sites=East$SiteID, ages=East$C14, h=100)
bins_North <- binPrep(sites=North$SiteID, ages=North$C14, h=100)
bins_South <- binPrep(sites=South$SiteID, ages=South$C14, h=100)
bins_Ebro <- binPrep(sites=Ebro_basin$SiteID, ages=Ebro_basin$C14, h=100)
bins_Tagus <- binPrep(sites=Tagus_basins$SiteID, ages=Tagus_basins$C14, h=100)
bins_Duero <- binPrep(sites=Duero_basin$SiteID, ages=Duero_basin$C14, h=100)

## Prepare sample Dates
# East
dates.randates_East <- sampleDates(dates_cal_East,bins= bins_East, nsim=1000, verbose=TRUE)
dates.ckde_East <- ckde(dates.randates_East, timeRange=c(7000,2000), bw=200)
# North
dates.randates_North <- sampleDates(dates_cal_North,bins= bins_North, nsim=1000, verbose=TRUE)
dates.ckde_North <- ckde(dates.randates_North, timeRange=c(7000,2000), bw=200)
# South
dates.randates_South <- sampleDates(dates_cal_South,bins= bins_South, nsim=1000, verbose=TRUE)
dates.ckde_South <- ckde(dates.randates_South, timeRange=c(7000,2000), bw=200)
# Ebro basin
dates.randates_Ebro <- sampleDates(dates_cal_Ebro,bins= bins_Ebro, nsim=1000, verbose=TRUE)
dates.ckde_Ebro <- ckde(dates.randates_Ebro, timeRange=c(7000,2000), bw=200)
# Tagus basins
dates.randates_Tagus <- sampleDates(dates_cal_Tagus,bins= bins_Tagus, nsim=1000, verbose=TRUE)
dates.ckde_Tagus <- ckde(dates.randates_Tagus, timeRange=c(7000,2000), bw=200)
# Duero basins
dates.randates_Duero <- sampleDates(dates_cal_Duero,bins= bins_Duero, nsim=1000, verbose=TRUE)
dates.ckde_Duero <- ckde(dates.randates_Duero, timeRange=c(7000,2000), bw=200)

# Plot Composite Kernel Density Estimates (CKDE) (Figure 4)
par(mar=c(4, 4, 4, 4), mfrow=c(2,3))
plot(dates.ckde_North, type="multiline", main="North")
plot(dates.ckde_Duero, type="multiline", main="Douro/Duero and Mondego basins")
plot(dates.ckde_Tagus, type="multiline", main="Tagus and western basins")
plot(dates.ckde_South, type="multiline", main="South")
plot(dates.ckde_East, type="multiline", main="East")
plot(dates.ckde_Ebro, type="multiline", main="Ebro basin")
par(mfrow=c(1,1))

## STACKED SPDs
# Read raw data
dates = read.csv2(here('data','C14dates_Iberia.csv'), sep=";",na='n/a')

# Subset dates by curve
dates_terrestrial <- subset(dates, Curve=="intcal20")
dates_marine <- subset(dates, Curve=="marine20")
dates_Mixed <- subset(dates, Curve=="mixed")

# Set Mixed curve
mixedCurve <- mixCurves('intcal20','marine20',p=0.8) 

# Calibrate dates using the corresponding curve
dates_Terr_cal <- calibrate(x=dates_terrestrial$C14,errors=dates_terrestrial$STD,normalised=TRUE, 
                            calMatrix=FALSE, calCurves = dates_terrestrial$Curve)

dates_Marine_cal <- calibrate(x=dates_marine$C14, errors=dates_marine$STD,
                              normalised=TRUE, calMatrix=FALSE, calCurves = dates_marine$Curve, 
                              resOffsets = dates_marine$DeltaR, resErrors = dates_marine$DeltaRErr)

dates_Mixed_cal <- calibrate(x=dates_Mixed$C14,
                             errors=dates_Mixed$STD,
                             normalised=TRUE, calMatrix=FALSE, calCurves = mixedCurve, 
                             resOffsets = dates_Mixed$DeltaR, 
                             resErrors = dates_Mixed$DeltaRErr)

# Combine calibrated dates
dates_cal = rcarbon::combine(dates_Terr_cal,dates_Marine_cal,dates_Mixed_cal, fixIDs = TRUE)

# Create the bin-width (h=100)
bins.meg <- binPrep(sites=dates$SiteID,
                    ages=dates$C14, h=100)

# SPD stacked visual comparison by regions
meg.spd=stackspd(x=dates_cal,
                 group=dates$Region,timeRange=c(7000,2000),
                 bins=bins.meg,runm=100)

# Plot Stacked SPDs (Figure 5)
par(mar=c(5.1, 4.1, 2, 5))
plot(meg.spd,type='multipanel',legend=F, gapFactor=0.3)
legend(x="bottomright", box.lty = 0, # Insert custom legend, to match the order in graph
       legend=c("East","South", "Ebro basin", "North","Tagus and western basins", "Douro/Duero and \nMondego basins"),
       cex=0.9, fill= c("#E6AB02", "#66A61E", "#E7298A", "#7570B3", "#D95F02","#1B9E77"))

##################################################
### PART 2. BAYESIAN MODELING AND SITE COUNTS ###
################################################

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
# On this simple plot you can see how we accounted for site duration

## B) SITE COUNTS

# Read bespoke R functions (located in "src" folder)
source(here('src','sampleCounts.R')) 

# Read MCMC output
x = read.csv(here('oxcalresults','MCMC_uniform.csv'))

# Eliminate pass number and last column
x = x[,-c(1,ncol(x))]

# Count the number of sites for a given chronological block, for each of the possible results
cntMatrix = sampleCounts(x) # This may take a couple of minutes

# Plot Sites count (Figure 6)
plotCntMat(cntMatrix)

## C) COMPARE EMPIRICAL SPDs AGAINS EACH OTHER

# Read raw data and study area
dates = read.csv(here('data','C14dates_Iberia.csv'), sep=";",na='n/a')
studyarea <- readOGR(dsn=here("shp","area.shp"), layer="area")
base=as.owin(studyarea)

# Subset dates by curve
dates_terrestrial <- subset(dates, Curve=="intcal20")
dates_marine <- subset(dates, Curve=="marine20")
dates_Mixed <- subset(dates, Curve=="mixed")

# Set Mixed curve
mixedCurve <- mixCurves('intcal20','marine20', p=0.8) 

# Calibrate dates using corresponding curve
dates_Terr_cal <- calibrate(x=dates_terrestrial$C14,errors=dates_terrestrial$STD,normalised=TRUE, 
                            calMatrix=FALSE, calCurves = dates_terrestrial$Curve)

dates_Marine_cal <- calibrate(x=dates_marine$C14, errors=dates_marine$STD,
                              normalised=TRUE, calMatrix=FALSE, calCurves = dates_marine$Curve, 
                              resOffsets = dates_marine$DeltaR, resErrors = dates_marine$DeltaRErr)

dates_Mixed_cal <- calibrate(x=dates_Mixed$C14,
                             errors=dates_Mixed$STD,
                             normalised=TRUE, calMatrix=FALSE, calCurves = mixedCurve, 
                             resOffsets = dates_Mixed$DeltaR, 
                             resErrors = dates_Mixed$DeltaRErr)

# Combine calibrated dates in one object
dates_cal = rcarbon::combine(dates_Terr_cal,dates_Marine_cal,dates_Mixed_cal, fixIDs = TRUE)

# Define number of simulations
nsim=5000

# Run permutation test
perm.meg=permTest(x=dates_cal,marks=dates$Region,timeRange=c(7000,2000),bins=bins1,
                  nsim=nsim,runm=100)

# Plot results (Figure 7)
par(mar=c(4, 4, 4, 4), mfrow=c(2,3))
plot(perm.meg,focalm = 1,main="East", xlim=c(7000, 2000), ylim=c(0, 0.10))
plot(perm.meg,focalm = 2,main="South", xlim=c(7000, 2000), ylim=c(0, 0.10))
plot(perm.meg,focalm = 3,main="North", xlim=c(7000, 2000), ylim=c(0, 0.10))
plot(perm.meg,focalm = 4,main="Ebro basin", xlim=c(7000, 2000), ylim=c(0, 0.10))
plot(perm.meg,focalm = 5,main="Tagus and western basins", xlim=c(7000, 2000), ylim=c(0, 0.10))
plot(perm.meg,focalm = 6,main="Douro/Duero and Mondego basins", xlim=c(7000, 2000), ylim=c(0, 0.10))
par(mfrow=c(1,1))

#####################################################################
### PART 3. OBSERVED TIME SERIES PER REGION AND PERMUTATION TEST ###
###################################################################

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

# Compute observed time series by region and run permutation test
res=markTest(x=x,marks=marks,nsim=50) # Note this takes several hours of processing with just 50 simulations. Increasing the number will raise the computational cost enormously.

# Plot results (Figure 8)
par(mar=c(5,4,2,1),mfrow=c(3,2))
for (i in 1:length(unique(marks)))
{
  plotMarkTest(res,index=i,main=unique(marks)[i])
  legend('topleft',legend=c(paste0('n=',sum(marks==unique(marks)[i])),paste0('P-Value=',round(res$pValueList[i],4))),bty='n',cex=1)
}

##############################################################################################
### PART 4. OBSERVED TIME SERIES PER TYPE OF MONUMENT (ARQUITECTURE) AND PERMUTATION TEST ###
############################################################################################

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

# Compute observed time series by region and run permutation test
res1=markTest(x=x,marks=marks,nsim=50) # Note this takes several hours of processing with just 50 simulations. Increasing the number will raise the computational cost enormously.

# Plot results (Figure 9)
par(mar=c(4,2,4,2),mfrow=c(2,4))
for (i in 1:length(unique(marks)))
{
  plotMarkTest(res1,index=i,main=unique(marks)[i])
  legend('topright',legend=c(paste0('n=',sum(marks==unique(marks)[i])),paste0('P-Value=',round(res1$pValueList[i],4))),bty='n',cex=1)
}


###########################################################
### PART 5. SPATIAL PERMUTATION TEST - FOCAL INTENSITY ###
#########################################################

# Read raw data and study area
dates = read.csv(here('data','C14dates_Iberia.csv'), sep=";",na='n/a')
studyarea <- readOGR(dsn=here("shp","area.shp"), layer="area")
base=as.owin(studyarea)

# Subset dates by curve
dates_terrestrial <- subset(dates, Curve=="intcal20")
dates_marine <- subset(dates, Curve=="marine20")
dates_Mixed <- subset(dates, Curve=="mixed")

# Set Mixed curve
mixedCurve <- mixCurves('intcal20','marine20', p=0.8) 

# Calibrate dates using corresponding curve
dates_Terr_cal <- calibrate(x=dates_terrestrial$C14,errors=dates_terrestrial$STD,normalised=TRUE, 
                            calMatrix=FALSE, calCurves = dates_terrestrial$Curve)

dates_Marine_cal <- calibrate(x=dates_marine$C14, errors=dates_marine$STD,
                              normalised=TRUE, calMatrix=FALSE, calCurves = dates_marine$Curve, 
                              resOffsets = dates_marine$DeltaR, resErrors = dates_marine$DeltaRErr)

dates_Mixed_cal <- calibrate(x=dates_Mixed$C14,
                             errors=dates_Mixed$STD,
                             normalised=TRUE, calMatrix=FALSE, calCurves = mixedCurve, 
                             resOffsets = dates_Mixed$DeltaR, 
                             resErrors = dates_Mixed$DeltaRErr)

# Combine calibrated dates in one object
dates_cal = rcarbon::combine(dates_Terr_cal,dates_Marine_cal,dates_Mixed_cal, fixIDs = TRUE)

# Create bin-width
bins1 <- binPrep(sites=dates$SiteID, ages=dates$C14, h=100)

# Map the spatio-temporal intensity of the radiocarbon dates between 7000-2000 (focal year = 500y)
stkde1 <- stkde(x=dates_cal, coords=dates[,c("Eastings_X","Northings_Y")], 
                win=base, sbw=40000, cellres=5000, 
                focalyears=seq(7000, 2000, -500), tbw=100, 
                bins=bins1, backsight=400, outdir = (here("SpatialPerm_Results")), 
                amount=1, verbose=TRUE)

# Plot results (Figure 10)
par(mar=c(0.2, 0.2, 2, 1),mfrow=c(3,4))
plot(stkde1, 7000, main="7000 calBP", cex.main=2)
plot(stkde1, 6500, main="6500 calBP", cex.main=2)
plot(stkde1, 6000, main="6000 calBP", cex.main=2)
plot(stkde1, 5500, main="5500 calBP", cex.main=2)
plot(stkde1, 5000, main="5000 calBP", cex.main=2)
plot(stkde1, 4500, main="4500 calBP", cex.main=2)
plot(stkde1, 4000, main="4000 calBP", cex.main=2)
plot(stkde1, 3500, main="3500 calBP", cex.main=2)
plot(stkde1, 3000, main="3000 calBP", cex.main=2)
plot(stkde1, 2500, main="2500 calBP", cex.main=2)
plot(stkde1, 2000, main="2000 calBP", cex.main=2)
par(mfrow=c(1,1))