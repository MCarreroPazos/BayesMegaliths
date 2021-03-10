##############################################
### FIGURE 7. PERMUTATION TEST BY REGIONS ###
############################################

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

# Comparing empirical SPDs against pan-regional trend
nsim=10000
perm.meg=permTest(x=dates_cal,marks=dates$Region,timeRange=c(7000,2000),bins=bins1,
                  nsim=nsim,runm=100)
summary(perm.meg) # to get specific positive and negative deviations

png(file=here('figures','Figure 7.png'),width = 1000,height=600)
par(mar=c(4, 4, 4, 4), mfrow=c(2,3))
plot(perm.meg,focalm = 1,main="Douro/Duero and Mondego basins", xlim=c(7000, 2000), ylim=c(0, 0.10))
plot(perm.meg,focalm = 2,main="Tagus and western basins", xlim=c(7000, 2000), ylim=c(0, 0.10))
plot(perm.meg,focalm = 3,main="North", xlim=c(7000, 2000), ylim=c(0, 0.10))
plot(perm.meg,focalm = 4,main="Ebro basin", xlim=c(7000, 2000), ylim=c(0, 0.10))
plot(perm.meg,focalm = 5,main="South", xlim=c(7000, 2000), ylim=c(0, 0.10))
plot(perm.meg,focalm = 6,main="East", xlim=c(7000, 2000), ylim=c(0, 0.10))
par(mfrow=c(1,1))
dev.off()