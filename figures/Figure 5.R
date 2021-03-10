##########################################
### Figure 5. COMPARING EMPIRICAL SPDs ##
########################################

library(maptools)
library(spatstat)
library(rgdal)
library(rcarbon)
library(RColorBrewer)
setwd(~BayesMegaliths-master) # Set first the working directory at the BayesMegaliths-master folder, as this script is coded to have that folder as the root
library(here)

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

# Make the plot and save it in "Figures" folder as Figure 5
png(file=here('figures','Figure 5.png'),width = 800,height=500)
par(mar=c(5.1, 4.1, 2, 5))
plot(meg.spd,type='multipanel',legend=F, gapFactor=0.3)
legend(x="bottomright", box.lty = 0, # Insert custom legend, to match the order in graph
       legend=c("East","South", "Ebro basin", "North","Tagus and western basins", "Douro/Duero and \nMondego basins"),
       cex=0.9, fill= c("#E6AB02", "#66A61E", "#E7298A", "#7570B3", "#D95F02","#1B9E77"))
dev.off()
