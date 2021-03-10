##############################################################
### Figure 10. SPATIAL PERMUTATION TEST - FOCAL INTENSITY ###
############################################################

# Load Required Libraries
library(readxl)
library(maptools)
library(spatstat)
library(rgdal)
library(rcarbon)
library(dplyr)
library(animation)
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

# Map the Spatio-temporal intensity of the 14C dates between 7000-2000 calBP (focal year = 500y) and save rfiles in the "SpatialPerm_Results" folder
stkde1 <- stkde(x=dates_cal, coords=dates[,c("Eastings_X","Northings_Y")], 
                win=base, sbw=40000, cellres=5000, 
                focalyears=seq(7000, 2000, -500), tbw=100, 
                bins=bins1, backsight=400, outdir = (here("SpatialPerm_Results")), 
                amount=1, verbose=TRUE)

# Plot figure 10 and save it to figures folder
png(file=here('figures','Figure 10.png'),width = 800,height=400)
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
dev.off()

# Bonus: Save results for all the focal years in an animated file (.gif). The file will be saved on current working directory "~/BayesMegaliths-master".
saveGIF(for (i in 11)
{
  par(mar=c(0.5, 0.5, 2.5, 2))
  plot(stkde1,7000, type="all")
  plot(stkde1,6500, type="all")
  plot(stkde1,6000, type="all")
  plot(stkde1,5500, type="all")
  plot(stkde1,5000, type="all")
  plot(stkde1,4500, type="all")
  plot(stkde1,4000, type="all")
  plot(stkde1,3500, type="all")
  plot(stkde1,3000, type="all")
  plot(stkde1,2500, type="all")
  plot(stkde1,2000, type="all")
},
movie.name = "Focal years of Kernel Density Maps.gif",
ani.width  = 800,
ani.height = 308,
ani.res = 100,
interval = 2)