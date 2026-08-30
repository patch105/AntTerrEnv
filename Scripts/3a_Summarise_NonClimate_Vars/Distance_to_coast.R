
###############################################
# Load libraries -----------------------------------------------------------
################################################

library(terra)
library(here)
library(sf)


# Load domains ------------------------------------------------------------

coast_domain <- rast(here("Data/coast_domain.tif"))
ice_free_domain <- rast(here("Data/ice_free_domain.tif"))

###############################################
# Distance to Antarctic coastline -------------------------------------------
################################################

# Using the Antarctic Digital Database high resolution seamask polygon v. 7.12 (Ireland, 2026)

# Load the seamask layer
seamask <- vect(here("Data/add_seamask_high_res_v7_12.shp"))
crs(seamask) <- "EPSG:3031"

# Dissolve into a single (multi-part) geometry to speed up distance calculation
seamask_agg <- terra::aggregate(seamask)

dist_coastline_coast <- terra::distance(coast_domain, seamask_agg)
dist_coastline_coast <- mask(dist_coastline_coast, coast_domain, maskvalue = NA)

dist_coastline_ice_free <- terra::distance(ice_free_domain, seamask_agg)
dist_coastline_ice_free <- mask(dist_coastline_ice_free, ice_free_domain, maskvalue = NA)

writeRaster(dist_coastline_coast, here("Outputs/NonClimate_Vars/dist_to_coast_seamask_v7_10_COASTLINE.tif"), overwrite = T)
writeRaster(dist_coastline_ice_free, here("Outputs/NonClimate_Vars/dist_to_coast_seamask_v7_10_ICEFREE.tif"), overwrite = T)


###############################################
# References --------------------------------------------------------------
################################################

# Ireland, L. (2026). High resolution vector polygon seamask for areas south of 60S (Version 7.12) [Data set]. NERC EDS UK Polar Data Centre. https://doi.org/10.5285/9460bb2e-4f1e-4b77-9b6a-d0da4c991aae [Accessed 02 July 2026]

