###########################################################
### Figure 4. COMPOSITE KERNEL DENSITY ESTIMATES (CKDE) ##
#########################################################

# Load libraries
library(maptools)
library(spatstat)
library(rgdal)
library(rcarbon)
setwd(~BayesMegaliths-master) # Set first the working directory at the BayesMegaliths-master folder, as this script is coded to have that folder as the root
library(here)

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

## Plot figure 4 and save it to figures folder 
png(file=here('figures','Figure 4.png'),width = 1000,height=600)
par(mar=c(4, 4, 4, 4), mfrow=c(2,3))
plot(dates.ckde_North, type="multiline", main="North")
plot(dates.ckde_Duero, type="multiline", main="Douro/Duero and Mondego basins")
plot(dates.ckde_Tagus, type="multiline", main="Tagus and western basins")
plot(dates.ckde_South, type="multiline", main="South")
plot(dates.ckde_East, type="multiline", main="East")
plot(dates.ckde_Ebro, type="multiline", main="Ebro basin")
par(mfrow=c(1,1))
dev.off()