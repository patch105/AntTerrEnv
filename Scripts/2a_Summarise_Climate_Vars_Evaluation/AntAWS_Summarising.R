
library(terra)

AntAWS_locs <- vect("C:/Users/n11222026/OneDrive - Queensland University of Technology/Code/AntTerrEnv/Data/antaws-dataset-x70w9q1u/AWS_location_shapefiles/Shp/267AWS.shp")

AntAWS_locs
AntAWS_locs <- project(AntAWS_locs, "EPSG:3031")
plot(AntAWS_locs)

coast_domain <- rast(here("C:/Users/n11222026/OneDrive - Queensland University of Technology/Code/AntTerrEnv/Data/coast_domain.tif"))
ice_free_domain <- rast(here("C:/Users/n11222026/OneDrive - Queensland University of Technology/Code/AntTerrEnv/Data/ice_free_domain.tif"))

# Find out how far AntAWS_locs are from ice-free areas (at 100m resolution)
dist_raster <- rast(here("C:/Users/n11222026/OneDrive - Queensland University of Technology/Code/AntTerrEnv/Data/Dist_to_Icefree_100m.tif"))

dist_raster <- distance(ice_free_domain)
plot(dist_raster)

AntAWS_locs$dist_to_ice_free <- terra::extract(dist_raster, AntAWS_locs)[, 2]

library(dplyr)

AntAWS_locs <- st_as_AntAWS_locs

test <- AntAWS_locs %>% filter(dist_to_ice_free == 0)

test <- AntAWS_locs$dist_to_ice_free
test0 <- test>0

sum(test>0)


