
###############################################
# Load libraries -----------------------------------------------------------
################################################

library(terra)
library(here)

# Load domains ------------------------------------------------------------

coast_domain <- rast(here("Data/coast_domain.tif"))

ice_free_domain <- rast(here("Data/ice_free_domain.tif"))

ice_free_domain_SPVE <- vect("C:/Users/n11222026/OneDrive - Queensland University of Technology/Code/AntTerrEnv/Data/rocks_Union_Land.shp")

# Data from the Reference Elevation Model for Antarctica (Howat, 2019)
# 100m filled mosaic
# URL: https://data.pgc.umn.edu/elev/dem/setsm/REMA/mosaic/latest/100m/rema_mosaic_100m_v2.0_filled_cop30.tar.gz

rema_100m <- rast(here("Data/rema_100m.tif"))
twi_100m <- rast(here("Data/twi_100m.tif"))

# twi_100m <- rast("C:/Users/n11222026/OneDrive - Queensland University of Technology/Data/raw/Gabi_100m_topography/twi_100m.tif")
# 
# bunger <- vect("C:/Users/n11222026/OneDrive - Queensland University of Technology/Code/Objective_3/AntarcticFutureHabitat/Data/Environmental_predictors/bunger_boundary.shp")
# 
# twi_100m <- crop(twi_100m, ext(bunger))
# twi_1km <- aggregate(twi_100m, fact = 10, fun = "mean", na.rm = T)
# 
# ice_free_domain_SPVE <- crop(ice_free_domain_SPVE, ext(bunger))

##########################################
# Make some 100m vars with terra -------------------------------------------
##########################################

fileBlocksize(rast("rema.tif")) # To find out block size for the below calculation

# Slope
terrain(rast(here("Data/rema_100m.tif")), v = "slope", neighbors = 8, filename = here("Data/slope_100m.tif"),
        gdal = c("TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512"))

# Aspect
terrain(rast(here("Data/rema_100m.tif")), v = "aspect", neighbors = 8, filename = here("Data/aspect_100m.tif"),
        gdal = c("TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512"))

# TRI
terrain(rast(here("Data/rema_100m.tif")), v = "TRI", filename = here("Data/TRI_100m.tif"),
        gdal = c("TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512"))

# TPI
terrain(rast(here("Data/rema_100m.tif")), v = "TPI", filename = here("Data/TPI_100m.tif"),
        gdal = c("TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512"))

# Roughness
terrain(rast(here("Data/rema_100m.tif")), v = "roughness", filename = here("Data/roughness_100m.tif"),
        gdal = c("TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512"))


# Now some derived variables from Amatulli et al. (2018) --------------------

slope_100m <- rast(here("Data/slope_100m.tif"))
aspect_100m <- rast(here("Data/aspect_100m.tif"))
  
northness_100m <- sin(slope_100m) * cos(aspect_100m)
eastness_100m <- sin(slope_100m) * sin(aspect_100m)

# Upscaling from 100m to 1km ----------------------------------------------

# Function to replace values in the original (not upscaled) 1km raster cells that are present in the 1km upscaled (from 100m) raster
overlayRaster <- function(raster1, raster2, epsg){
  if (ext(raster1) != ext(raster2)) {
    raster1 <- resample(raster1, raster2)
  }
  
  condition <- !is.na(raster1) # where the non NA ice-free area values are
  mask <- terra::mask(raster2, condition, maskvalues = TRUE) # remove the cells where the new ice-free layer will be
  
  raster_merged <- rast(nrows = nrow(raster2), ncols = ncol(raster2), crs=paste0("epsg:", epsg)) # create empty raster
  ext(raster_merged) <- ext(raster2)
  
  stopifnot(ext(raster_merged) == ext(raster2))
  
  raster_merged[] <- ifelse(is.na(mask[]), raster1[], mask[])
  return(raster_merged)
}


# Slope
slope_icefree_100m <- mask(slope_100m, ice_free_domain_SPVE)
slope_icefree_1km <- aggregate(slope_icefree_100m, fact = 10, fun = "mean", na.rm = T)
slope_1km <- aggregate(slope_100m, fact = 10, fun = "mean", na.rm = T)
slope_icefree_1km_merged <- overlayRaster(slope_icefree_1km, slope_1km, 3031)

# Aspect
# TRI
# TPI
# Roughness
# TWI




##########################################
# Format some vars pre-processed in QGIS ----------------------------------
##########################################

# Slope 
# Aspect
# TWI




test <- rast("C:/Users/n11222026/OneDrive - Queensland University of Technology/Data/raw/Gabi_100m_topography/twi_100m.tif")
oasis_boundary <- vect(here("C:/Users/n11222026/OneDrive - Queensland University of Technology/Code/Objective_3/AntarcticFutureHabitat/Data/Environmental_predictors/bunger_boundary.shp"))
tet <- crop(test)


###############################################
# References --------------------------------------------------------------
################################################

# Howat, I.M., Porter, C., Smith, B.E., Noh, M.-J., Morin, P., 2019. The Reference Elevation 
# Model of Antarctica. The Cryosphere 13, 665–674. https://doi.org/10.5194/tc-13-665
# 2019 