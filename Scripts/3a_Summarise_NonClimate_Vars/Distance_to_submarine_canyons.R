
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
# Distance to submarine canyons -------------------------------------------
################################################

# Using a submarine canyon boundary layer from Santora et al. (2020)
# Downloaded from https://doi.org/10.7291/D1NT0S
canyons <- vect(here("Data/doi_10_7291_D1NT0S__v20210211/MarineGeology.shp/Canyons.shp"))
canyons <- project(canyons, "EPSG:3031")

dist_canyons_coast <- terra::distance(coast_domain, canyons)
dist_canyons_coast <- mask(dist_canyons_coast, coast_domain, maskvalue = NA)

dist_canyons_ice_free <- terra::distance(ice_free_domain, canyons)
dist_canyons_ice_free <- mask(dist_canyons_ice_free, ice_free_domain, maskvalue = NA)

writeRaster(dist_canyons_coast, here("Outputs/NonClimate_Vars/dist_to_submarine_canyons_COASTLINE.tif"), overwrite = T)
writeRaster(dist_canyons_ice_free, here("Outputs/NonClimate_Vars/dist_to_submarine_canyons_ICEFREE.tif"), overwrite = T)


###############################################
# References --------------------------------------------------------------
################################################

# Santora, J.A., LaRue, M.A., Ainley, D.G., 2020. Geographic structuring of Antarctic penguin populations. Global Ecology and Biogeography 29, 1716–1728. https://doi.org/10.1111/geb.13144
