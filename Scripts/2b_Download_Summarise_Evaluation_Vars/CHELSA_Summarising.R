
library(dplyr)
library(purrr)
library(terra)
library(here)


# Set inpath -------------------------------------------------------------

# inpath <- file.path("Z://AntarcticFutureHabitat/Data/CHELSA")
# inpath <- here("Data/CHELSA")
inpath <- "/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/CHELSA"

# Set outpath -------------------------------------------------------------

outpath <- "/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/CHELSA/Validation"

# Load ice-free domain ----------------------------------------------------

domain <- rast(here("Data/Environmental_predictors/ice_free_upsamp_1km.tif"))

# Set domain values
domain <- ifel(not.na(domain), 1, NA)


###################################
# HISTORICAL (1981 - 2010) ------------------------------------------------
###################################

hist_files <- list.files(inpath, pattern = "CHELSA_tas_\\d{2}_1981-2010_V\\.2\\.1\\.tif$", 
                         full.names = TRUE) %>% 
  sort() # ensures months 01-12 order

chelsa_hist <- rast(hist_files)

# Name layers by month
names(chelsa_hist) <- sprintf("tas_%02d", 1:12)

# Crop to Antarctic extent
bbox_4326 <- ext(-180, 180, -90, -60)
chelsa_hist <- crop(chelsa_hist, bbox_4326)

# Reproject to Antarctic
chelsa_hist <- project(chelsa_hist, "epsg:3031") 


# Annual temp -------------------------------------------------------------

mean_annual_temp_K <- app(chelsa_hist, mean, na.rm = TRUE)
# Convert from Kelvin to Celsius
mean_annual_temp <- mean_annual_temp_K - 273.15

writeRaster(mean_annual_temp, file.path(outpath, "Mean_Annual_Temperature_HISTORICAL_1981_2010.tif"), overwrite = T)


# Crop to ice-free areas --------------------------------------------------

mean_annual_temp_IF <- terra::project(mean_annual_temp, domain, method = "bilinear")
mean_annual_temp_IF <- mask(mean_annual_temp_IF, domain, maskvalue = NA)
writeRaster(mean_annual_temp_IF, file.path(outpath, "Mean_Annual_Temperature_HISTORICAL_1981_2010_ICEFREE.tif"), overwrite = T)

##########################
# Summer temp (DJF) --------------------------------------------------------
##########################

chelsa_hist_DJF <- chelsa_hist[[c(1, 2, 12)]]

mean_summer_temp_K <- app(chelsa_hist_DJF, mean, na.rm = TRUE)
# Convert from Kelvin to Celsius
mean_summer_temp <- mean_summer_temp_K - 273.15

writeRaster(mean_summer_temp, file.path(outpath, "Mean_Summer_Temperature_HISTORICAL_1981_2010.tif"), overwrite = T)


# Crop to ice-free areas --------------------------------------------------

mean_summer_temp_IF <- terra::project(mean_summer_temp, domain, method = "bilinear")
mean_summer_temp_IF <- mask(mean_summer_temp_IF, domain, maskvalue = NA)
writeRaster(mean_summer_temp_IF, file.path(outpath, "Mean_Summer_Temperature_HISTORICAL_1981_2010_ICEFREE.tif"), overwrite = T)


######################
# Winter temp (DJF) --------------------------------------------------------
######################

chelsa_hist_JJA <- chelsa_hist[[c(6, 7, 8)]]

mean_winter_temp_K <- app(chelsa_hist_JJA, mean, na.rm = TRUE)
# Convert from Kelvin to Celsius
mean_winter_temp <- mean_winter_temp_K - 273.15

writeRaster(mean_winter_temp, file.path(outpath, "Mean_Winter_Temperature_HISTORICAL_1981_2010.tif"), overwrite = T)


# Crop to ice-free areas --------------------------------------------------

mean_winter_temp_IF <- terra::project(mean_winter_temp, domain, method = "bilinear")
mean_winter_temp_IF <- mask(mean_winter_temp_IF, domain, maskvalue = NA)
writeRaster(mean_winter_temp_IF, file.path(outpath, "Mean_Winter_Temperature_HISTORICAL_1981_2010_ICEFREE.tif"), overwrite = T)



###################################
# FUTURE (2071 - 2100) ------------------------------------------------
###################################

# Extract arguments from command line
args <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])

future_models  <- c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0",  "ukesm1-0-ll")

model <- future_models[job_index]

future_files <- list.files(inpath, pattern = paste0("CHELSA_", model, "_r1i1p1f1_w5e5_ssp370_tas_\\d{2}_2071-2100_V\\.2\\.1\\.tif$"),
                           full.names = TRUE)

chelsa_future <- rast(future_files)

# Name layers by month
names(chelsa_future) <- sprintf("tas_%02d", 1:12)

# Crop to Antarctic extent
bbox_4326 <- ext(-180, 180, -90, -60)
chelsa_future <- crop(chelsa_future, bbox_4326)

# Reproject to Antarctic
chelsa_future <- project(chelsa_future, "epsg:3031") 


# Annual temp -------------------------------------------------------------
mean_annual_temp_K <- app(chelsa_future, mean, na.rm = TRUE)
# Convert from Kelvin to Celsius
mean_annual_temp <- mean_annual_temp_K - 273.15
writeRaster(mean_annual_temp, file.path(outpath, paste0("Mean_Annual_Temperature_FUTURE_", model, "_2071_2100.tif")), overwrite = T)
# Crop to ice-free areas --------------------------------------------------
mean_annual_temp_IF <- terra::project(mean_annual_temp, domain, method = "bilinear")
mean_annual_temp_IF <- mask(mean_annual_temp_IF, domain, maskvalue = NA)
writeRaster(mean_annual_temp_IF, file.path(outpath, paste0("Mean_Annual_Temperature_FUTURE_", model, "_2071_2100_ICEFREE.tif")), overwrite = T)

##########################
# Summer temp (DJF) --------------------------------------------------------
##########################
chelsa_future_DJF <- chelsa_future[[c(1, 2, 12)]]
mean_summer_temp_K <- app(chelsa_future_DJF, mean, na.rm = TRUE)
# Convert from Kelvin to Celsius
mean_summer_temp <- mean_summer_temp_K - 273.15
writeRaster(mean_summer_temp, file.path(outpath, paste0("Mean_Summer_Temperature_FUTURE_", model, "_2071_2100.tif")), overwrite = T)
# Crop to ice-free areas --------------------------------------------------
mean_summer_temp_IF <- terra::project(mean_summer_temp, domain, method = "bilinear")
mean_summer_temp_IF <- mask(mean_summer_temp_IF, domain, maskvalue = NA)
writeRaster(mean_summer_temp_IF, file.path(outpath, paste0("Mean_Summer_Temperature_FUTURE_", model, "_2071_2100_ICEFREE.tif")), overwrite = T)

######################
# Winter temp (JJA) --------------------------------------------------------
######################
chelsa_future_JJA <- chelsa_future[[c(6, 7, 8)]]
mean_winter_temp_K <- app(chelsa_future_JJA, mean, na.rm = TRUE)
# Convert from Kelvin to Celsius
mean_winter_temp <- mean_winter_temp_K - 273.15
writeRaster(mean_winter_temp, file.path(outpath, paste0("Mean_Winter_Temperature_FUTURE_", model, "_2071_2100.tif")), overwrite = T)
# Crop to ice-free areas --------------------------------------------------
mean_winter_temp_IF <- terra::project(mean_winter_temp, domain, method = "bilinear")
mean_winter_temp_IF <- mask(mean_winter_temp_IF, domain, maskvalue = NA)
writeRaster(mean_winter_temp_IF, file.path(outpath, paste0("Mean_Winter_Temperature_FUTURE_", model, "_2071_2100_ICEFREE.tif")), overwrite = T)
