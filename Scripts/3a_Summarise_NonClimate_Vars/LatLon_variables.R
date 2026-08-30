
###############################################
# Load libraries -----------------------------------------------------------
################################################

library(terra)
library(here)

# Load domains ------------------------------------------------------------

coast_domain <- rast(here("Data/coast_domain.tif"))

ice_free_domain <- rast(here("Data/ice_free_domain.tif"))

# Get coords for every grid cell
coords <- xyFromCell(coast_domain, 1:ncell(coast_domain)) 
# coords_geo <- project(coords, from = crs(coast_domain), to = "EPSG:4326")

###################################################
# Latitude ----------------------------------------------------------------
###################################################

latitude <- setValues(coast_domain, coords[,2]) 
names(latitude) <- "Latitude" 

# Mask to non-NA grid cells
latitude_coast <- mask(latitude, coast_domain, maskvalue = NA)
latitude_ice_free <- mask(latitude, ice_free_domain, maskvalue = NA)

writeRaster(latitude_coast, here("Outputs/NonClimate_Vars/latitude_COASTLINE.tif"), overwrite = T)
writeRaster(latitude_ice_free, here("Outputs/NonClimate_Vars/latitude_ICEFREE.tif"), overwrite = T)


###################################################
# Longitude ----------------------------------------------------------------
###################################################

longitude <- setValues(coast_domain, coords[,1]) 
names(longitude) <- "Longitude" 

# Mask to non-NA grid cells
longitude_coast <- mask(longitude, coast_domain, maskvalue = NA)
longitude_ice_free <- mask(longitude, ice_free_domain, maskvalue = NA)

writeRaster(longitude_coast, here("Outputs/NonClimate_Vars/longitude_COASTLINE.tif"), overwrite = T)
writeRaster(longitude_ice_free, here("Outputs/NonClimate_Vars/longitude_ICEFREE.tif"), overwrite = T)

