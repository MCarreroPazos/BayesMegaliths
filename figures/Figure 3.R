########################################################################
### Figure 3. DISTRIBUTION OF 14C DATES AND KERNEL DENSITY ANALYSIS ###
######################################################################

# Load Required Libraries
library("ggplot2")
library("ggpubr")
library("ggspatial")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
setwd(~BayesMegaliths-master) # Set first the working directory at the BayesMegaliths-master folder, as this script is coded to have that folder as the root
library("here")

# Read and prepare data
dates = read.csv2(here('data','C14dates_Iberia.csv'), sep=";",na='n/a')

# Set black and white theme and load data
theme_set(theme_bw())
sites <- as.data.frame(dates[,c("SiteID", "Eastings_X", "Northings_Y", "Region")])
polys <- st_read(here("shp", "Regions.shp")) # Defined regions
area <- st_read(here("shp", "area.shp")) # Study area (Iberia and Balearic Islands)

# Plot A (Study area with sites)
plot1 <- ggplot(data = area) +
        xlab("Longitude") + ylab("Latitude") +
        annotation_scale(location = "br", width_hint = 0.3) + #Set scale bar
        annotation_north_arrow(location = "br", height = unit(0.3, "in"), #Set north arrow
                         width = unit(0.4, "in"),
                         which_north = "true", 
                         pad_x = unit(0.12, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_minimal)+
        geom_point(data = sites, #Plot points (from C14dates_Iberia.csv)
            aes(x = Eastings_X, y = Northings_Y, 
            colour = Region)) + ggtitle("A") + 
        geom_sf(data = polys, colour = "gray47", fill = NA, lwd=0.5) 

# Plot B (2D kernel density estimation)
plot2 <- ggplot(data = area) +
        geom_sf() +
        xlab("Longitude") + ylab("Latitude") +
        annotation_scale(location = "br", width_hint = 0.3) + #Set scale bar
        annotation_north_arrow(location = "br", height = unit(0.3, "in"), #Set north arrow
                         width = unit(0.4, "in"),
                         which_north = "true", 
                         pad_x = unit(0.12, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_minimal)+
        geom_density_2d_filled(data = sites,
              aes(x = Eastings_X, y = Northings_Y), alpha=0.5, contour_var = "ndensity")+
        geom_point(data = sites, #Plot points (from C14dates_Iberia.csv)
             aes(x = Eastings_X, y = Northings_Y)) + ggtitle("B") + theme(legend.position='none')

# Arrange both plots in 1
ggarrange(plot1, plot2, nrow=1, ncol=2, widths = c(3,2))

# Save Figure 3 to a png file under "figures" folder
ggsave(here("figures","Figure 1.png"), width = 12, height = 5)
dev.off()