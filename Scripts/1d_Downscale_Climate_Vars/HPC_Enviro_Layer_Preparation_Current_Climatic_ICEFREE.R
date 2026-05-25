
library(dplyr)
library(purrr)
library(terra)
library(here)

######################################################
############ MAR-MPI-ESM Layer final preparation #########
######################################################

annual_temp_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Temperature_HISTORICAL_1995_2014.tif"))
annual_temp_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Temperature_FUTURE_2081_2100.tif"))


summer_temp_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Summer_Temperature_HISTORICAL_1995_2014.tif"))
summer_temp_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Summer_Temperature_FUTURE_2081_2100.tif"))

annual_precip_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Total_Annual_Precipitation_HISTORICAL_1995_2014.tif"))
annual_precip_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Total_Annual_Precipitation_FUTURE_2081_2100.tif"))

summer_precip_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Summer_Total_Precipitation_HISTORICAL_1995_2014.tif"))
summer_precip_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Summer_Total_Precipitation_FUTURE_2081_2100.tif"))


nov_sea_ice_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_November_Sea_Ice_Concentration_HISTORICAL_1995_2014.tif"))
nov_sea_ice_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_November_Sea_Ice_Concentration_FUTURE_2081_2100.tif"))

nov_sea_ice_100km_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_November_Sea_Ice_Concentration_100km_HISTORICAL_1995_2014.tif"))
nov_sea_ice_100km_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_November_Sea_Ice_Concentration_100km_FUTURE_2081_2100.tif"))

# snow_area_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Snow_Area_Percentage_HISTORICAL_1995_2014.tif"))
# snow_area_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Snow_Area_Percentage_FUTURE_2081_2100.tif"))
# 
# summer_snow_area_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Summer_Snow_Area_Percentage_HISTORICAL_1995_2014.tif"))
# summer_snow_area_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Summer_Snow_Area_Percentage_FUTURE_2081_2100.tif"))

degree_days_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Total_Degree_Days-5_HISTORICAL_1995_2014.tif"))
degree_days_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Total_Degree_Days-5_FUTURE_2081_2100.tif"))


wind_speed_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Wind_Speed_HISTORICAL_1995_2014.tif"))
wind_speed_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Wind_Speed_FUTURE_2081_2100.tif"))


solar_rad_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Solar_Radiation_HISTORICAL_1995_2014.tif"))
solar_rad_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Solar_Radiation_FUTURE_2081_2100.tif"))


winter_temp_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Winter_Temperature_HISTORICAL_1995_2014.tif"))
winter_temp_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Winter_Temperature_FUTURE_2081_2100.tif"))


octfeb_sea_ice_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Oct_Feb_Sea_Ice_Concentration_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Oct_Feb_Sea_Ice_Concentration_FUTURE_2081_2100.tif"))

octfeb_sea_ice_32km_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Oct_Feb_Sea_Ice_Concentration_32km_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_32km_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Oct_Feb_Sea_Ice_Concentration_32km_FUTURE_2081_2100.tif"))

octfeb_sea_ice_36km_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Oct_Feb_Sea_Ice_Concentration_36km_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_36km_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Oct_Feb_Sea_Ice_Concentration_36km_FUTURE_2081_2100.tif"))

octfeb_sea_ice_42km_hist <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Oct_Feb_Sea_Ice_Concentration_42km_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_42km_future <- rast(here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Oct_Feb_Sea_Ice_Concentration_42km_FUTURE_2081_2100.tif"))


# Spatial interpolation and cropping to ice-free --------------------------

# Load target domain
domain <- rast(here("Data/Environmental_predictors/ice_free_upsamp_1km.tif"))

# Set domain values
domain <- ifel(not.na(domain), 1, NA)

annual_temp_hist <- terra::project(annual_temp_hist, domain, method = "near")
annual_temp_hist <- mask(annual_temp_hist, domain, maskvalue = NA)
writeRaster(annual_temp_hist, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

annual_temp_future <- terra::project(annual_temp_future, domain, method = "near")
annual_temp_future <- mask(annual_temp_future, domain, maskvalue = NA)
writeRaster(annual_temp_future, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


summer_temp_hist <- terra::project(summer_temp_hist, domain, method = "near")
summer_temp_hist <- mask(summer_temp_hist, domain, maskvalue = NA)
writeRaster(summer_temp_hist, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Summer_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

summer_temp_future <- terra::project(summer_temp_future, domain, method = "near")
summer_temp_future <- mask(summer_temp_future, domain, maskvalue = NA)
writeRaster(summer_temp_future, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Summer_Temperature_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


annual_precip_hist <- terra::project(annual_precip_hist, domain, method = "near")
annual_precip_hist <- mask(annual_precip_hist, domain, maskvalue = NA)
writeRaster(annual_precip_hist, here("Data/Environmental_predictors/MAR_MPI_ESM/Total_Annual_Precipitation_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

annual_precip_future <- terra::project(annual_precip_future, domain, method = "near")
annual_precip_future <- mask(annual_precip_future, domain, maskvalue = NA)
writeRaster(annual_precip_future, here("Data/Environmental_predictors/MAR_MPI_ESM/Total_Annual_Precipitation_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


# snow_area_hist <- terra::project(snow_area_hist, domain, method = "near")
# snow_area_hist <- mask(snow_area_hist, domain, maskvalue = NA)
# writeRaster(snow_area_hist, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Snow_Area_Percentage_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)
# 
# snow_area_future <- terra::project(snow_area_future, domain, method = "near")
# snow_area_future <- mask(snow_area_future, domain, maskvalue = NA)
# writeRaster(snow_area_future, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Snow_Area_Percentage_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


degree_days_hist <- terra::project(degree_days_hist, domain, method = "near")
degree_days_hist <- mask(degree_days_hist, domain, maskvalue = NA)
writeRaster(degree_days_hist, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Total_Degree_Days-5_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

degree_days_future <- terra::project(degree_days_future, domain, method = "near")
degree_days_future <- mask(degree_days_future, domain, maskvalue = NA)
writeRaster(degree_days_future, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Total_Degree_Days-5_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)



wind_speed_hist <- terra::project(wind_speed_hist, domain, method = "near")
wind_speed_hist <- mask(wind_speed_hist, domain, maskvalue = NA)
writeRaster(wind_speed_hist, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Wind_Speed_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

wind_speed_future <- terra::project(wind_speed_future, domain, method = "near")
wind_speed_future <- mask(wind_speed_future, domain, maskvalue = NA)
writeRaster(wind_speed_future, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Wind_Speed_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


solar_rad_hist <- terra::project(solar_rad_hist, domain, method = "near")
solar_rad_hist <- mask(solar_rad_hist, domain, maskvalue = NA)
writeRaster(solar_rad_hist, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Solar_Radiation_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

solar_rad_future <- terra::project(solar_rad_future, domain, method = "near")
solar_rad_future <- mask(solar_rad_future, domain, maskvalue = NA)
writeRaster(solar_rad_future, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Annual_Solar_Radiation_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


# nov_sea_ice_100km_hist <- terra::project(nov_sea_ice_100km_hist, domain, method = "near")
# nov_sea_ice_100km_hist <- mask(nov_sea_ice_100km_hist, domain, maskvalue = NA)
writeRaster(nov_sea_ice_100km_hist, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_November_Sea_Ice_Concentration_100km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# nov_sea_ice_100km_future <- terra::project(nov_sea_ice_100km_future, domain, method = "near")
# nov_sea_ice_100km_future <- mask(nov_sea_ice_100km_future, domain, maskvalue = NA)
writeRaster(nov_sea_ice_100km_future, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_November_Sea_Ice_Concentration_100km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


winter_temp_hist <- terra::project(winter_temp_hist, domain, method = "near")
winter_temp_hist <- mask(winter_temp_hist, domain, maskvalue = NA)
writeRaster(winter_temp_hist, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Winter_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

winter_temp_future <- terra::project(winter_temp_future, domain, method = "near")
winter_temp_future <- mask(winter_temp_future, domain, maskvalue = NA)
writeRaster(winter_temp_future, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Winter_Temperature_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


summer_precip_hist <- terra::project(summer_precip_hist, domain, method = "near")
summer_precip_hist <- mask(summer_precip_hist, domain, maskvalue = NA)
writeRaster(summer_precip_hist, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Summer_Total_Precipitation_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

summer_precip_future <- terra::project(summer_precip_future, domain, method = "near")
summer_precip_future <- mask(summer_precip_future, domain, maskvalue = NA)
writeRaster(summer_precip_future, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Summer_Total_Precipitation_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)



# summer_snow_area_hist <- terra::project(summer_snow_area_hist, domain, method = "near")
# summer_snow_area_hist <- mask(summer_snow_area_hist, domain, maskvalue = NA)
# writeRaster(summer_snow_area_hist, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Summer_Snow_Area_Percentage_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)
# 
# summer_snow_area_future <- terra::project(summer_snow_area_future, domain, method = "near")
# summer_snow_area_future <- mask(summer_snow_area_future, domain, maskvalue = NA)
# writeRaster(summer_snow_area_future, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Summer_Snow_Area_Percentage_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)




# octfeb_sea_ice_32km_hist <- terra::project(octfeb_sea_ice_32km_hist, domain, method = "near")
# octfeb_sea_ice_32km_hist <- mask(octfeb_sea_ice_32km_hist, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_32km_hist, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Oct_Feb_Sea_Ice_Concentration_32km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# octfeb_sea_ice_32km_future <- terra::project(octfeb_sea_ice_32km_future, domain, method = "near")
# octfeb_sea_ice_32km_future <- mask(octfeb_sea_ice_32km_future, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_32km_future, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Oct_Feb_Sea_Ice_Concentration_32km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


# octfeb_sea_ice_36km_hist <- terra::project(octfeb_sea_ice_36km_hist, domain, method = "near")
# octfeb_sea_ice_36km_hist <- mask(octfeb_sea_ice_36km_hist, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_36km_hist, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Oct_Feb_Sea_Ice_Concentration_36km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# octfeb_sea_ice_36km_future <- terra::project(octfeb_sea_ice_36km_future, domain, method = "near")
# octfeb_sea_ice_36km_future <- mask(octfeb_sea_ice_36km_future, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_36km_future, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Oct_Feb_Sea_Ice_Concentration_36km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


# octfeb_sea_ice_42km_hist <- terra::project(octfeb_sea_ice_42km_hist, domain, method = "near")
# octfeb_sea_ice_42km_hist <- mask(octfeb_sea_ice_42km_hist, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_42km_hist, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Oct_Feb_Sea_Ice_Concentration_42km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# octfeb_sea_ice_42km_future <- terra::project(octfeb_sea_ice_42km_future, domain, method = "near")
# octfeb_sea_ice_42km_future <- mask(octfeb_sea_ice_42km_future, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_42km_future, here("Data/Environmental_predictors/MAR_MPI_ESM/Mean_Oct_Feb_Sea_Ice_Concentration_42km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)



######################################################
############ MAR-CESM2 Layer final preparation #########
######################################################

annual_temp_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Temperature_HISTORICAL_1995_2014.tif"))
annual_temp_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Temperature_FUTURE_2081_2100.tif"))


summer_temp_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Summer_Temperature_HISTORICAL_1995_2014.tif"))
summer_temp_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Summer_Temperature_FUTURE_2081_2100.tif"))

annual_precip_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Total_Annual_Precipitation_HISTORICAL_1995_2014.tif"))
annual_precip_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Total_Annual_Precipitation_FUTURE_2081_2100.tif"))

summer_precip_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Summer_Total_Precipitation_HISTORICAL_1995_2014.tif"))
summer_precip_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Summer_Total_Precipitation_FUTURE_2081_2100.tif"))


nov_sea_ice_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_November_Sea_Ice_Concentration_HISTORICAL_1995_2014.tif"))
nov_sea_ice_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_November_Sea_Ice_Concentration_FUTURE_2081_2100.tif"))

nov_sea_ice_100km_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_November_Sea_Ice_Concentration_100km_HISTORICAL_1995_2014.tif"))
nov_sea_ice_100km_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_November_Sea_Ice_Concentration_100km_FUTURE_2081_2100.tif"))

# snow_area_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Snow_Area_Percentage_HISTORICAL_1995_2014.tif"))
# snow_area_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Snow_Area_Percentage_FUTURE_2081_2100.tif"))

# summer_snow_area_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Summer_Snow_Area_Percentage_HISTORICAL_1995_2014.tif"))
# summer_snow_area_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Summer_Snow_Area_Percentage_FUTURE_2081_2100.tif"))

degree_days_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Total_Degree_Days-5_HISTORICAL_1995_2014.tif"))
degree_days_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Total_Degree_Days-5_FUTURE_2081_2100.tif"))


wind_speed_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Wind_Speed_HISTORICAL_1995_2014.tif"))
wind_speed_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Wind_Speed_FUTURE_2081_2100.tif"))


solar_rad_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Solar_Radiation_HISTORICAL_1995_2014.tif"))
solar_rad_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Solar_Radiation_FUTURE_2081_2100.tif"))


winter_temp_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Winter_Temperature_HISTORICAL_1995_2014.tif"))
winter_temp_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Winter_Temperature_FUTURE_2081_2100.tif"))



octfeb_sea_ice_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_FUTURE_2081_2100.tif"))

octfeb_sea_ice_32km_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_32km_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_32km_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_32km_FUTURE_2081_2100.tif"))

octfeb_sea_ice_36km_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_36km_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_36km_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_36km_FUTURE_2081_2100.tif"))

octfeb_sea_ice_42km_hist <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_42km_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_42km_future <- rast(here("Data/Environmental_predictors/MAR_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_42km_FUTURE_2081_2100.tif"))


# Spatial interpolation and cropping to ice-free --------------------------

# Load target domain
domain <- rast(here("Data/Environmental_predictors/ice_free_upsamp_1km.tif"))

# Set domain values
domain <- ifel(not.na(domain), 1, NA)

annual_temp_hist <- terra::project(annual_temp_hist, domain, method = "near")
annual_temp_hist <- mask(annual_temp_hist, domain, maskvalue = NA)
writeRaster(annual_temp_hist, here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

annual_temp_future <- terra::project(annual_temp_future, domain, method = "near")
annual_temp_future <- mask(annual_temp_future, domain, maskvalue = NA)
writeRaster(annual_temp_future, here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


summer_temp_hist <- terra::project(summer_temp_hist, domain, method = "near")
summer_temp_hist <- mask(summer_temp_hist, domain, maskvalue = NA)
writeRaster(summer_temp_hist, here("Data/Environmental_predictors/MAR_CESM2/Mean_Summer_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

summer_temp_future <- terra::project(summer_temp_future, domain, method = "near")
summer_temp_future <- mask(summer_temp_future, domain, maskvalue = NA)
writeRaster(summer_temp_future, here("Data/Environmental_predictors/MAR_CESM2/Mean_Summer_Temperature_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


annual_precip_hist <- terra::project(annual_precip_hist, domain, method = "near")
annual_precip_hist <- mask(annual_precip_hist, domain, maskvalue = NA)
writeRaster(annual_precip_hist, here("Data/Environmental_predictors/MAR_CESM2/Total_Annual_Precipitation_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

annual_precip_future <- terra::project(annual_precip_future, domain, method = "near")
annual_precip_future <- mask(annual_precip_future, domain, maskvalue = NA)
writeRaster(annual_precip_future, here("Data/Environmental_predictors/MAR_CESM2/Total_Annual_Precipitation_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


# snow_area_hist <- terra::project(snow_area_hist, domain, method = "near")
# snow_area_hist <- mask(snow_area_hist, domain, maskvalue = NA)
# writeRaster(snow_area_hist, here("Data/Environmental_predictors/MAR_CESM2/Mean_Snow_Area_Percentage_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)
# 
# snow_area_future <- terra::project(snow_area_future, domain, method = "near")
# snow_area_future <- mask(snow_area_future, domain, maskvalue = NA)
# writeRaster(snow_area_future, here("Data/Environmental_predictors/MAR_CESM2/Mean_Snow_Area_Percentage_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


degree_days_hist <- terra::project(degree_days_hist, domain, method = "near")
degree_days_hist <- mask(degree_days_hist, domain, maskvalue = NA)
writeRaster(degree_days_hist, here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Total_Degree_Days-5_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

degree_days_future <- terra::project(degree_days_future, domain, method = "near")
degree_days_future <- mask(degree_days_future, domain, maskvalue = NA)
writeRaster(degree_days_future, here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Total_Degree_Days-5_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)



wind_speed_hist <- terra::project(wind_speed_hist, domain, method = "near")
wind_speed_hist <- mask(wind_speed_hist, domain, maskvalue = NA)
writeRaster(wind_speed_hist, here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Wind_Speed_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

wind_speed_future <- terra::project(wind_speed_future, domain, method = "near")
wind_speed_future <- mask(wind_speed_future, domain, maskvalue = NA)
writeRaster(wind_speed_future, here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Wind_Speed_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


solar_rad_hist <- terra::project(solar_rad_hist, domain, method = "near")
solar_rad_hist <- mask(solar_rad_hist, domain, maskvalue = NA)
writeRaster(solar_rad_hist, here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Solar_Radiation_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

solar_rad_future <- terra::project(solar_rad_future, domain, method = "near")
solar_rad_future <- mask(solar_rad_future, domain, maskvalue = NA)
writeRaster(solar_rad_future, here("Data/Environmental_predictors/MAR_CESM2/Mean_Annual_Solar_Radiation_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


nov_sea_ice_100km_hist <- terra::project(nov_sea_ice_100km_hist, domain, method = "near")
nov_sea_ice_100km_hist <- mask(nov_sea_ice_100km_hist, domain, maskvalue = NA)
writeRaster(nov_sea_ice_100km_hist, here("Data/Environmental_predictors/MAR_CESM2/Mean_November_Sea_Ice_Concentration_100km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

nov_sea_ice_100km_future <- terra::project(nov_sea_ice_100km_future, domain, method = "near")
nov_sea_ice_100km_future <- mask(nov_sea_ice_100km_future, domain, maskvalue = NA)
writeRaster(nov_sea_ice_100km_future, here("Data/Environmental_predictors/MAR_CESM2/Mean_November_Sea_Ice_Concentration_100km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


winter_temp_hist <- terra::project(winter_temp_hist, domain, method = "near")
winter_temp_hist <- mask(winter_temp_hist, domain, maskvalue = NA)
writeRaster(winter_temp_hist, here("Data/Environmental_predictors/MAR_CESM2/Mean_Winter_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

winter_temp_future <- terra::project(winter_temp_future, domain, method = "near")
winter_temp_future <- mask(winter_temp_future, domain, maskvalue = NA)
writeRaster(winter_temp_future, here("Data/Environmental_predictors/MAR_CESM2/Mean_Winter_Temperature_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


summer_precip_hist <- terra::project(summer_precip_hist, domain, method = "near")
summer_precip_hist <- mask(summer_precip_hist, domain, maskvalue = NA)
writeRaster(summer_precip_hist, here("Data/Environmental_predictors/MAR_CESM2/Mean_Summer_Total_Precipitation_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

summer_precip_future <- terra::project(summer_precip_future, domain, method = "near")
summer_precip_future <- mask(summer_precip_future, domain, maskvalue = NA)
writeRaster(summer_precip_future, here("Data/Environmental_predictors/MAR_CESM2/Mean_Summer_Total_Precipitation_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)



# summer_snow_area_hist <- terra::project(summer_snow_area_hist, domain, method = "near")
# summer_snow_area_hist <- mask(summer_snow_area_hist, domain, maskvalue = NA)
# writeRaster(summer_snow_area_hist, here("Data/Environmental_predictors/MAR_CESM2/Mean_Summer_Snow_Area_Percentage_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)
# 
# summer_snow_area_future <- terra::project(summer_snow_area_future, domain, method = "near")
# summer_snow_area_future <- mask(summer_snow_area_future, domain, maskvalue = NA)
# writeRaster(summer_snow_area_future, here("Data/Environmental_predictors/MAR_CESM2/Mean_Summer_Snow_Area_Percentage_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)



# octfeb_sea_ice_32km_hist <- terra::project(octfeb_sea_ice_32km_hist, domain, method = "near")
# octfeb_sea_ice_32km_hist <- mask(octfeb_sea_ice_32km_hist, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_32km_hist, here("Data/Environmental_predictors/MAR_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_32km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# octfeb_sea_ice_32km_future <- terra::project(octfeb_sea_ice_32km_future, domain, method = "near")
# octfeb_sea_ice_32km_future <- mask(octfeb_sea_ice_32km_future, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_32km_future, here("Data/Environmental_predictors/MAR_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_32km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


# octfeb_sea_ice_36km_hist <- terra::project(octfeb_sea_ice_36km_hist, domain, method = "near")
# octfeb_sea_ice_36km_hist <- mask(octfeb_sea_ice_36km_hist, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_36km_hist, here("Data/Environmental_predictors/MAR_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_36km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# octfeb_sea_ice_36km_future <- terra::project(octfeb_sea_ice_36km_future, domain, method = "near")
# octfeb_sea_ice_36km_future <- mask(octfeb_sea_ice_36km_future, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_36km_future, here("Data/Environmental_predictors/MAR_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_36km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


# octfeb_sea_ice_42km_hist <- terra::project(octfeb_sea_ice_42km_hist, domain, method = "near")
# octfeb_sea_ice_42km_hist <- mask(octfeb_sea_ice_42km_hist, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_42km_hist, here("Data/Environmental_predictors/MAR_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_42km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# octfeb_sea_ice_42km_future <- terra::project(octfeb_sea_ice_42km_future, domain, method = "near")
# octfeb_sea_ice_42km_future <- mask(octfeb_sea_ice_42km_future, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_42km_future, here("Data/Environmental_predictors/MAR_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_42km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)



######################################################
############ HCLIM-MPI-ESM Layer final preparation #########
######################################################

annual_temp_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Temperature_HISTORICAL_1995_2014.tif"))
annual_temp_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Temperature_FUTURE_2081_2100.tif"))


summer_temp_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Summer_Temperature_HISTORICAL_1995_2014.tif"))
summer_temp_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Summer_Temperature_FUTURE_2081_2100.tif"))

annual_precip_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Total_Annual_Precipitation_HISTORICAL_1995_2014.tif"))
annual_precip_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Total_Annual_Precipitation_FUTURE_2081_2100.tif"))

summer_precip_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Summer_Total_Precipitation_HISTORICAL_1995_2014.tif"))
summer_precip_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Summer_Total_Precipitation_FUTURE_2081_2100.tif"))


nov_sea_ice_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_November_Sea_Ice_Concentration_HISTORICAL_1995_2014.tif"))
nov_sea_ice_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_November_Sea_Ice_Concentration_FUTURE_2081_2100.tif"))

nov_sea_ice_100km_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_November_Sea_Ice_Concentration_100km_HISTORICAL_1995_2014.tif"))
nov_sea_ice_100km_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_November_Sea_Ice_Concentration_100km_FUTURE_2081_2100.tif"))

snow_area_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Snow_Cover_HISTORICAL_1995_2014.tif"))
snow_area_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Snow_Cover_FUTURE_2081_2100.tif"))

summer_snow_area_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Summer_Snow_Cover_HISTORICAL_1995_2014.tif"))
summer_snow_area_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Summer_Snow_Cover_FUTURE_2081_2100.tif"))

degree_days_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Total_Degree_Days-5_HISTORICAL_1995_2014.tif"))
degree_days_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Total_Degree_Days-5_FUTURE_2081_2100.tif"))


wind_speed_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Wind_Speed_HISTORICAL_1995_2014.tif"))
wind_speed_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Wind_Speed_FUTURE_2081_2100.tif"))


solar_rad_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Solar_Radiation_HISTORICAL_1995_2014.tif"))
solar_rad_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Solar_Radiation_FUTURE_2081_2100.tif"))


winter_temp_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Winter_Temperature_HISTORICAL_1995_2014.tif"))
winter_temp_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Winter_Temperature_FUTURE_2081_2100.tif"))


octfeb_sea_ice_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Oct_Feb_Sea_Ice_Concentration_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Oct_Feb_Sea_Ice_Concentration_FUTURE_2081_2100.tif"))

octfeb_sea_ice_32km_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Oct_Feb_Sea_Ice_Concentration_32km_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_32km_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Oct_Feb_Sea_Ice_Concentration_32km_FUTURE_2081_2100.tif"))

octfeb_sea_ice_36km_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Oct_Feb_Sea_Ice_Concentration_36km_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_36km_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Oct_Feb_Sea_Ice_Concentration_36km_FUTURE_2081_2100.tif"))

octfeb_sea_ice_42km_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Oct_Feb_Sea_Ice_Concentration_42km_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_42km_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Oct_Feb_Sea_Ice_Concentration_42km_FUTURE_2081_2100.tif"))

mean_melt_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Melt_HISTORICAL_1995_2014.tif"))
mean_melt_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Melt_FUTURE_2081_2100.tif"))

total_melt_hist <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Total_Annual_Melt_HISTORICAL_1995_2014.tif"))
total_melt_future <- rast(here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Total_Annual_Melt_FUTURE_2081_2100.tif"))




# Spatial interpolation and cropping to ice-free --------------------------

# Load target domain
domain <- rast(here("Data/Environmental_predictors/ice_free_upsamp_1km.tif"))

# Set domain values
domain <- ifel(not.na(domain), 1, NA)

annual_temp_hist <- terra::project(annual_temp_hist, domain, method = "near")
annual_temp_hist <- mask(annual_temp_hist, domain, maskvalue = NA)
writeRaster(annual_temp_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

annual_temp_future <- terra::project(annual_temp_future, domain, method = "near")
annual_temp_future <- mask(annual_temp_future, domain, maskvalue = NA)
writeRaster(annual_temp_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)




summer_temp_hist <- terra::project(summer_temp_hist, domain, method = "near")
summer_temp_hist <- mask(summer_temp_hist, domain, maskvalue = NA)
writeRaster(summer_temp_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Summer_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

summer_temp_future <- terra::project(summer_temp_future, domain, method = "near")
summer_temp_future <- mask(summer_temp_future, domain, maskvalue = NA)
writeRaster(summer_temp_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Summer_Temperature_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


annual_precip_hist <- terra::project(annual_precip_hist, domain, method = "near")
annual_precip_hist <- mask(annual_precip_hist, domain, maskvalue = NA)
writeRaster(annual_precip_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Total_Annual_Precipitation_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

annual_precip_future <- terra::project(annual_precip_future, domain, method = "near")
annual_precip_future <- mask(annual_precip_future, domain, maskvalue = NA)
writeRaster(annual_precip_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Total_Annual_Precipitation_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


snow_area_hist <- terra::project(snow_area_hist, domain, method = "near")
snow_area_hist <- mask(snow_area_hist, domain, maskvalue = NA)
writeRaster(snow_area_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Snow_Area_Percentage_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

snow_area_future <- terra::project(snow_area_future, domain, method = "near")
snow_area_future <- mask(snow_area_future, domain, maskvalue = NA)
writeRaster(snow_area_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Snow_Area_Percentage_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


degree_days_hist <- terra::project(degree_days_hist, domain, method = "near")
degree_days_hist <- mask(degree_days_hist, domain, maskvalue = NA)
writeRaster(degree_days_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Total_Degree_Days-5_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

degree_days_future <- terra::project(degree_days_future, domain, method = "near")
degree_days_future <- mask(degree_days_future, domain, maskvalue = NA)
writeRaster(degree_days_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Total_Degree_Days-5_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)



wind_speed_hist <- terra::project(wind_speed_hist, domain, method = "near")
wind_speed_hist <- mask(wind_speed_hist, domain, maskvalue = NA)
writeRaster(wind_speed_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Wind_Speed_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

wind_speed_future <- terra::project(wind_speed_future, domain, method = "near")
wind_speed_future <- mask(wind_speed_future, domain, maskvalue = NA)
writeRaster(wind_speed_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Wind_Speed_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


solar_rad_hist <- terra::project(solar_rad_hist, domain, method = "near")
solar_rad_hist <- mask(solar_rad_hist, domain, maskvalue = NA)
writeRaster(solar_rad_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Solar_Radiation_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

solar_rad_future <- terra::project(solar_rad_future, domain, method = "near")
solar_rad_future <- mask(solar_rad_future, domain, maskvalue = NA)
writeRaster(solar_rad_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Solar_Radiation_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


# nov_sea_ice_100km_hist <- terra::project(nov_sea_ice_100km_hist, domain, method = "near")
# nov_sea_ice_100km_hist <- mask(nov_sea_ice_100km_hist, domain, maskvalue = NA)
writeRaster(nov_sea_ice_100km_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_November_Sea_Ice_Concentration_100km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# nov_sea_ice_100km_future <- terra::project(nov_sea_ice_100km_future, domain, method = "near")
# nov_sea_ice_100km_future <- mask(nov_sea_ice_100km_future, domain, maskvalue = NA)
writeRaster(nov_sea_ice_100km_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_November_Sea_Ice_Concentration_100km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


winter_temp_hist <- terra::project(winter_temp_hist, domain, method = "near")
winter_temp_hist <- mask(winter_temp_hist, domain, maskvalue = NA)
writeRaster(winter_temp_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Winter_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

winter_temp_future <- terra::project(winter_temp_future, domain, method = "near")
winter_temp_future <- mask(winter_temp_future, domain, maskvalue = NA)
writeRaster(winter_temp_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Winter_Temperature_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


summer_precip_hist <- terra::project(summer_precip_hist, domain, method = "near")
summer_precip_hist <- mask(summer_precip_hist, domain, maskvalue = NA)
writeRaster(summer_precip_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Summer_Total_Precipitation_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

summer_precip_future <- terra::project(summer_precip_future, domain, method = "near")
summer_precip_future <- mask(summer_precip_future, domain, maskvalue = NA)
writeRaster(summer_precip_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Summer_Total_Precipitation_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)



summer_snow_area_hist <- terra::project(summer_snow_area_hist, domain, method = "near")
summer_snow_area_hist <- mask(summer_snow_area_hist, domain, maskvalue = NA)
writeRaster(summer_snow_area_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Summer_Snow_Area_Percentage_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

summer_snow_area_future <- terra::project(summer_snow_area_future, domain, method = "near")
summer_snow_area_future <- mask(summer_snow_area_future, domain, maskvalue = NA)
writeRaster(summer_snow_area_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Summer_Snow_Area_Percentage_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)




# octfeb_sea_ice_32km_hist <- terra::project(octfeb_sea_ice_32km_hist, domain, method = "near")
# octfeb_sea_ice_32km_hist <- mask(octfeb_sea_ice_32km_hist, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_32km_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Oct_Feb_Sea_Ice_Concentration_32km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# octfeb_sea_ice_32km_future <- terra::project(octfeb_sea_ice_32km_future, domain, method = "near")
# octfeb_sea_ice_32km_future <- mask(octfeb_sea_ice_32km_future, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_32km_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Oct_Feb_Sea_Ice_Concentration_32km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


# octfeb_sea_ice_36km_hist <- terra::project(octfeb_sea_ice_36km_hist, domain, method = "near")
# octfeb_sea_ice_36km_hist <- mask(octfeb_sea_ice_36km_hist, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_36km_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Oct_Feb_Sea_Ice_Concentration_36km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# octfeb_sea_ice_36km_future <- terra::project(octfeb_sea_ice_36km_future, domain, method = "near")
# octfeb_sea_ice_36km_future <- mask(octfeb_sea_ice_36km_future, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_36km_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Oct_Feb_Sea_Ice_Concentration_36km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


# octfeb_sea_ice_42km_hist <- terra::project(octfeb_sea_ice_42km_hist, domain, method = "near")
# octfeb_sea_ice_42km_hist <- mask(octfeb_sea_ice_42km_hist, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_42km_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Oct_Feb_Sea_Ice_Concentration_42km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# octfeb_sea_ice_42km_future <- terra::project(octfeb_sea_ice_42km_future, domain, method = "near")
# octfeb_sea_ice_42km_future <- mask(octfeb_sea_ice_42km_future, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_42km_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Oct_Feb_Sea_Ice_Concentration_42km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


mean_melt_hist <- terra::project(mean_melt_hist, domain, method = "near")
mean_melt_hist <- mask(mean_melt_hist, domain, maskvalue = NA)
writeRaster(mean_melt_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Melt_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

mean_melt_future <- terra::project(mean_melt_future, domain, method = "near")
mean_melt_future <- mask(mean_melt_future, domain, maskvalue = NA)
writeRaster(mean_melt_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Mean_Annual_Melt_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)

total_melt_hist <- terra::project(total_melt_hist, domain, method = "near")
total_melt_hist <- mask(total_melt_hist, domain, maskvalue = NA)
writeRaster(total_melt_hist, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Total_Annual_Melt_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

total_melt_future <- terra::project(total_melt_future, domain, method = "near")
total_melt_future <- mask(total_melt_future, domain, maskvalue = NA)
writeRaster(total_melt_future, here("Data/Environmental_predictors/HCLIM_MPI_ESM1/Total_Annual_Melt_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


######################################################
############ HCLIM-CESM2 Layer final preparation #########
######################################################

annual_temp_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Temperature_HISTORICAL_1995_2014.tif"))
annual_temp_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Temperature_FUTURE_2081_2100.tif"))


summer_temp_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Summer_Temperature_HISTORICAL_1995_2014.tif"))
summer_temp_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Summer_Temperature_FUTURE_2081_2100.tif"))

annual_precip_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Total_Annual_Precipitation_HISTORICAL_1995_2014.tif"))
annual_precip_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Total_Annual_Precipitation_FUTURE_2081_2100.tif"))

summer_precip_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Summer_Total_Precipitation_HISTORICAL_1995_2014.tif"))
summer_precip_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Summer_Total_Precipitation_FUTURE_2081_2100.tif"))


nov_sea_ice_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_November_Sea_Ice_Concentration_HISTORICAL_1995_2014.tif"))
nov_sea_ice_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_November_Sea_Ice_Concentration_FUTURE_2081_2100.tif"))

nov_sea_ice_100km_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_November_Sea_Ice_Concentration_100km_HISTORICAL_1995_2014.tif"))
nov_sea_ice_100km_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_November_Sea_Ice_Concentration_100km_FUTURE_2081_2100.tif"))

snow_area_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Snow_Cover_HISTORICAL_1995_2014.tif"))
snow_area_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Snow_Cover_FUTURE_2081_2100.tif"))

summer_snow_area_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Summer_Snow_Cover_HISTORICAL_1995_2014.tif"))
summer_snow_area_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Summer_Snow_Cover_FUTURE_2081_2100.tif"))

degree_days_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Total_Degree_Days-5_HISTORICAL_1995_2014.tif"))
degree_days_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Total_Degree_Days-5_FUTURE_2081_2100.tif"))


wind_speed_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Wind_Speed_HISTORICAL_1995_2014.tif"))
wind_speed_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Wind_Speed_FUTURE_2081_2100.tif"))


solar_rad_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Solar_Radiation_HISTORICAL_1995_2014.tif"))
solar_rad_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Solar_Radiation_FUTURE_2081_2100.tif"))


winter_temp_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Winter_Temperature_HISTORICAL_1995_2014.tif"))
winter_temp_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Winter_Temperature_FUTURE_2081_2100.tif"))


octfeb_sea_ice_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_FUTURE_2081_2100.tif"))

octfeb_sea_ice_32km_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_32km_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_32km_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_32km_FUTURE_2081_2100.tif"))

octfeb_sea_ice_36km_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_36km_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_36km_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_36km_FUTURE_2081_2100.tif"))

octfeb_sea_ice_42km_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_42km_HISTORICAL_1995_2014.tif"))
octfeb_sea_ice_42km_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_42km_FUTURE_2081_2100.tif"))

mean_melt_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Melt_HISTORICAL_1995_2014.tif"))
mean_melt_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Melt_FUTURE_2081_2100.tif"))

total_melt_hist <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Total_Annual_Melt_HISTORICAL_1995_2014.tif"))
total_melt_future <- rast(here("Data/Environmental_predictors/HCLIM_CESM2/Total_Annual_Melt_FUTURE_2081_2100.tif"))




# Spatial interpolation and cropping to ice-free --------------------------

# Load target domain
domain <- rast(here("Data/Environmental_predictors/ice_free_upsamp_1km.tif"))

# Set domain values
domain <- ifel(not.na(domain), 1, NA)

annual_temp_hist <- terra::project(annual_temp_hist, domain, method = "near")
annual_temp_hist <- mask(annual_temp_hist, domain, maskvalue = NA)
writeRaster(annual_temp_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

annual_temp_future <- terra::project(annual_temp_future, domain, method = "near")
annual_temp_future <- mask(annual_temp_future, domain, maskvalue = NA)
writeRaster(annual_temp_future, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


summer_temp_hist <- terra::project(summer_temp_hist, domain, method = "near")
summer_temp_hist <- mask(summer_temp_hist, domain, maskvalue = NA)
writeRaster(summer_temp_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Summer_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

summer_temp_future <- terra::project(summer_temp_future, domain, method = "near")
summer_temp_future <- mask(summer_temp_future, domain, maskvalue = NA)
writeRaster(summer_temp_future, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Summer_Temperature_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


annual_precip_hist <- terra::project(annual_precip_hist, domain, method = "near")
annual_precip_hist <- mask(annual_precip_hist, domain, maskvalue = NA)
writeRaster(annual_precip_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Total_Annual_Precipitation_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

annual_precip_future <- terra::project(annual_precip_future, domain, method = "near")
annual_precip_future <- mask(annual_precip_future, domain, maskvalue = NA)
writeRaster(annual_precip_future, here("Data/Environmental_predictors/HCLIM_CESM2/Total_Annual_Precipitation_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


snow_area_hist <- terra::project(snow_area_hist, domain, method = "near")
snow_area_hist <- mask(snow_area_hist, domain, maskvalue = NA)
writeRaster(snow_area_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Snow_Area_Percentage_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

snow_area_future <- terra::project(snow_area_future, domain, method = "near")
snow_area_future <- mask(snow_area_future, domain, maskvalue = NA)
writeRaster(snow_area_future, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Snow_Area_Percentage_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


degree_days_hist <- terra::project(degree_days_hist, domain, method = "near")
degree_days_hist <- mask(degree_days_hist, domain, maskvalue = NA)
writeRaster(degree_days_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Total_Degree_Days-5_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

degree_days_future <- terra::project(degree_days_future, domain, method = "near")
degree_days_future <- mask(degree_days_future, domain, maskvalue = NA)
writeRaster(degree_days_future, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Total_Degree_Days-5_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)



wind_speed_hist <- terra::project(wind_speed_hist, domain, method = "near")
wind_speed_hist <- mask(wind_speed_hist, domain, maskvalue = NA)
writeRaster(wind_speed_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Wind_Speed_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

wind_speed_future <- terra::project(wind_speed_future, domain, method = "near")
wind_speed_future <- mask(wind_speed_future, domain, maskvalue = NA)
writeRaster(wind_speed_future, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Wind_Speed_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


solar_rad_hist <- terra::project(solar_rad_hist, domain, method = "near")
solar_rad_hist <- mask(solar_rad_hist, domain, maskvalue = NA)
writeRaster(solar_rad_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Solar_Radiation_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

solar_rad_future <- terra::project(solar_rad_future, domain, method = "near")
solar_rad_future <- mask(solar_rad_future, domain, maskvalue = NA)
writeRaster(solar_rad_future, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Solar_Radiation_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


# nov_sea_ice_100km_hist <- terra::project(nov_sea_ice_100km_hist, domain, method = "near")
# nov_sea_ice_100km_hist <- mask(nov_sea_ice_100km_hist, domain, maskvalue = NA)
writeRaster(nov_sea_ice_100km_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_November_Sea_Ice_Concentration_100km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# nov_sea_ice_100km_future <- terra::project(nov_sea_ice_100km_future, domain, method = "near")
# nov_sea_ice_100km_future <- mask(nov_sea_ice_100km_future, domain, maskvalue = NA)
writeRaster(nov_sea_ice_100km_future, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_November_Sea_Ice_Concentration_100km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


winter_temp_hist <- terra::project(winter_temp_hist, domain, method = "near")
winter_temp_hist <- mask(winter_temp_hist, domain, maskvalue = NA)
writeRaster(winter_temp_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Winter_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

winter_temp_future <- terra::project(winter_temp_future, domain, method = "near")
winter_temp_future <- mask(winter_temp_future, domain, maskvalue = NA)
writeRaster(winter_temp_future, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Winter_Temperature_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


summer_precip_hist <- terra::project(summer_precip_hist, domain, method = "near")
summer_precip_hist <- mask(summer_precip_hist, domain, maskvalue = NA)
writeRaster(summer_precip_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Summer_Total_Precipitation_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

summer_precip_future <- terra::project(summer_precip_future, domain, method = "near")
summer_precip_future <- mask(summer_precip_future, domain, maskvalue = NA)
writeRaster(summer_precip_future, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Summer_Total_Precipitation_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)



summer_snow_area_hist <- terra::project(summer_snow_area_hist, domain, method = "near")
summer_snow_area_hist <- mask(summer_snow_area_hist, domain, maskvalue = NA)
writeRaster(summer_snow_area_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Summer_Snow_Area_Percentage_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

summer_snow_area_future <- terra::project(summer_snow_area_future, domain, method = "near")
summer_snow_area_future <- mask(summer_snow_area_future, domain, maskvalue = NA)
writeRaster(summer_snow_area_future, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Summer_Snow_Area_Percentage_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)




# octfeb_sea_ice_32km_hist <- terra::project(octfeb_sea_ice_32km_hist, domain, method = "near")
# octfeb_sea_ice_32km_hist <- mask(octfeb_sea_ice_32km_hist, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_32km_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_32km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# octfeb_sea_ice_32km_future <- terra::project(octfeb_sea_ice_32km_future, domain, method = "near")
# octfeb_sea_ice_32km_future <- mask(octfeb_sea_ice_32km_future, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_32km_future, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_32km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


# octfeb_sea_ice_36km_hist <- terra::project(octfeb_sea_ice_36km_hist, domain, method = "near")
# octfeb_sea_ice_36km_hist <- mask(octfeb_sea_ice_36km_hist, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_36km_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_36km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# octfeb_sea_ice_36km_future <- terra::project(octfeb_sea_ice_36km_future, domain, method = "near")
# octfeb_sea_ice_36km_future <- mask(octfeb_sea_ice_36km_future, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_36km_future, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_36km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


# octfeb_sea_ice_42km_hist <- terra::project(octfeb_sea_ice_42km_hist, domain, method = "near")
# octfeb_sea_ice_42km_hist <- mask(octfeb_sea_ice_42km_hist, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_42km_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_42km_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

# octfeb_sea_ice_42km_future <- terra::project(octfeb_sea_ice_42km_future, domain, method = "near")
# octfeb_sea_ice_42km_future <- mask(octfeb_sea_ice_42km_future, domain, maskvalue = NA)
writeRaster(octfeb_sea_ice_42km_future, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Oct_Feb_Sea_Ice_Concentration_42km_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)


mean_melt_hist <- terra::project(mean_melt_hist, domain, method = "near")
mean_melt_hist <- mask(mean_melt_hist, domain, maskvalue = NA)
writeRaster(mean_melt_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Melt_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

mean_melt_future <- terra::project(mean_melt_future, domain, method = "near")
mean_melt_future <- mask(mean_melt_future, domain, maskvalue = NA)
writeRaster(mean_melt_future, here("Data/Environmental_predictors/HCLIM_CESM2/Mean_Annual_Melt_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)

total_melt_hist <- terra::project(total_melt_hist, domain, method = "near")
total_melt_hist <- mask(total_melt_hist, domain, maskvalue = NA)
writeRaster(total_melt_hist, here("Data/Environmental_predictors/HCLIM_CESM2/Total_Annual_Melt_HISTORICAL_1995_2014_ICEFREE.tif"), overwrite = T)

total_melt_future <- terra::project(total_melt_future, domain, method = "near")
total_melt_future <- mask(total_melt_future, domain, maskvalue = NA)
writeRaster(total_melt_future, here("Data/Environmental_predictors/HCLIM_CESM2/Total_Annual_Melt_FUTURE_2081_2100_ICEFREE.tif"), overwrite = T)







# 
# #################################################
# ############ Polar Res Interpolating 1km #########
# ################################################
# 
# # Temperature is in K
# annual_temp <- rast(here("Data/Environmental_predictors/Mean_Annual_Temperature_ALL_YEARS.tif"))
# annual_temp <- annual_temp - 273.15 # Convert to Celsius
# 
# summer_temp <- rast(here("Data/Environmental_predictors/Mean_Summer_Temperature_2001_2018.tif"))
# summer_temp <- summer_temp - 273.15 # Convert to Celsius
# 
# # Degree days above -5
# deg_day <- rast(here("Data/Environmental_predictors/Mean_of_Total_Annual_Degree_Days-5_ALL_YEARS.tif"))
# 
# # Precip is in kg m-2 s-1
# # The units "kg m⁻² s⁻¹" represent precipitation as a mass flux per unit area per unit time, meaning the amount of water (in kilograms) that falls over a given surface area (in square meters) within a specified time period (in seconds). Essentially, it describes how much water is falling on a particular area within a given time frame
# 
# annual_precip_snow <- rast(here("Data/Environmental_predictors/Mean_of_Total_Annual_Precipitation_SNOW_ALL_YEARS.tif"))
# # annual_precip_snow <- annual_precip_snow * 60 * 60 * 24 * 365.25 # Convert to mm/year
# 
# annual_precip_rain <- rast(here("Data/Environmental_predictors/Mean_of_Total_Annual_Precipitation_RAIN_ALL_YEARS.tif"))
# # annual_precip_rain <- annual_precip_rain * 60 * 60 * 24 * 365.25 # Convert to mm/year
# 
# annual_precip <- annual_precip_snow + annual_precip_rain
# 
# # Wind speed is in m/s
# annual_wind <- rast(here("Data/Environmental_predictors/Mean_Annual_Wind_Speed_ALL_YEARS.tif"))
# 
# 
# #### Spatial interpolation ####
# 
# # Load target domain
# domain <- rast(here("Data/Environmental_predictors/ice_free_upsamp_1km.tif"))
# 
# # Set domain values
# domain <- ifel(not.na(domain), 1, NA)
# 
# annual_temp <- terra::project(annual_temp, domain, method = "near")
# annual_temp <- mask(annual_temp, domain, maskvalue = NA)
# writeRaster(annual_temp, here("Data/Environmental_predictors/Mean_Annual_Temperature_ALL_YEARS_1km.tif"), overwrite = T)
# 
# summer_temp <- terra::project(summer_temp, domain, method = "near")
# summer_temp <- mask(summer_temp, domain, maskvalue = NA)
# writeRaster(summer_temp, here("Data/Environmental_predictors/Mean_Summer_Temperature_ALL_YEARS_1km.tif"), overwrite = T)
# 
# deg_day <- terra::project(deg_day, domain, method = "near")
# deg_day <- mask(deg_day, domain, maskvalue = NA)
# writeRaster(deg_day, here("Data/Environmental_predictors/Mean_of_Total_Annual_Degree_Days-5_ALL_YEARS_1km.tif"), overwrite = T)
# 
# annual_wind <- terra::project(annual_wind, domain, method = "near")
# annual_wind <- mask(annual_wind, domain, maskvalue = NA)
# writeRaster(annual_wind, here("Data/Environmental_predictors/Mean_Annual_Wind_Speed_ALL_YEARS_1km.tif"), overwrite = T)
# 
# annual_precip_snow <- terra::project(annual_precip_snow, domain, method = "near")
# annual_precip_snow <- mask(annual_precip_snow, domain, maskvalue = NA)
# writeRaster(annual_precip_snow, here("Data/Environmental_predictors/Mean_of_Total_Annual_Precipitation_SNOW_ALL_YEARS_1km.tif"), overwrite = T)
# 
# annual_precip_rain <- terra::project(annual_precip_rain, domain, method = "near")
# annual_precip_rain <- mask(annual_precip_rain, domain, maskvalue = NA)
# writeRaster(annual_precip_rain, here("Data/Environmental_predictors/Mean_of_Total_Annual_Precipitation_RAIN_ALL_YEARS_1km.tif"), overwrite = T)
# 
# annual_precip <- terra::project(annual_precip, domain, method = "near")
# annual_precip <- mask(annual_precip, domain, maskvalue = NA)
# writeRaster(annual_precip, here("Data/Environmental_predictors/Mean_of_Total_Annual_Precipitation_ALL_YEARS_1km.tif"), overwrite = T)
# 
# ######################
# 
# # terra::density(annual_temp, maxcells=ncell(annual_temp), plot = T)
# # 
# # 
# # test <- project(annual_temp, domain, method = "cubic")
# # test <- mask(test, domain, maskvalue = NA)
# # 
# # terra::density(test,plot = T)
# # 
# # writeRaster(test, "test3.tif", overwrite = T)
# # writeRaster(annual_temp, "test2.tif", overwrite = T)
# # 
# # plot(annual_temp)
# # plot(test)
# # 
# # 
# # 
# # # Load elevation
# # elevation <- rast(here("Data/Environmental_predictors/elevation_ICEFREE.tif"))
# # 
# # devtools::install_github("https://github.com/ErikKusch/KrigR")
# # library(KrigR)
# # 
# 
# 
# 
# # ###############################################
# # ######## Sea Ice Extent and Concentration #########
# # ################################################
# 
# ## Data accessed 19th March 2025
# 
# # Sea Ice Index, Version 3 https://nsidc.org/data/g02135/versions/3#anchor-data-access-tools
# # Downloaded from: https://noaadata.apps.nsidc.org/NOAA/G02135/
# 
# ## Script to download in R from : https://nsidc.org/data/user-resources/help-center/how-access-and-download-noaansidc-data
# 
# library(rvest)
# 
# setwd("C:/Users/n11222026/OneDrive - Queensland University of Technology/Data/raw/Sea_Ice")
# 
# url <- list("https://noaadata.apps.nsidc.org/NOAA/G02135/south/monthly/geotiff/12_Dec/",
#             "https://noaadata.apps.nsidc.org/NOAA/G02135/south/monthly/geotiff/11_Nov/",
#             "https://noaadata.apps.nsidc.org/NOAA/G02135/south/monthly/geotiff/10_Oct/",
#             "https://noaadata.apps.nsidc.org/NOAA/G02135/south/monthly/geotiff/09_Sep/",
#             "https://noaadata.apps.nsidc.org/NOAA/G02135/south/monthly/geotiff/08_Aug/",
#             "https://noaadata.apps.nsidc.org/NOAA/G02135/south/monthly/geotiff/07_Jul/",
#             "https://noaadata.apps.nsidc.org/NOAA/G02135/south/monthly/geotiff/06_Jun/",
#             "https://noaadata.apps.nsidc.org/NOAA/G02135/south/monthly/geotiff/05_May/",
#             "https://noaadata.apps.nsidc.org/NOAA/G02135/south/monthly/geotiff/04_Apr/",
#             "https://noaadata.apps.nsidc.org/NOAA/G02135/south/monthly/geotiff/03_Mar/",
#             "https://noaadata.apps.nsidc.org/NOAA/G02135/south/monthly/geotiff/02_Feb/",
#             "https://noaadata.apps.nsidc.org/NOAA/G02135/south/monthly/geotiff/01_Jan/"
# )
# 
# map(url, function(x){
#   
#   page <- read_html(x)
#   
#   files <- page %>% 
#     html_nodes("a") %>% 
#     html_attr("href")
#   
#   for(i in 2:length(files)){
#     
#     u <- paste(x, files[i], sep="/")
#     
#     download.file(u, files[i], mode = "wb")
#     
#   }
#   
# })
# 
# 
# ## Loading the saved files:
# 
# years <- as.character(1979:2021)
# 
# pattern <- paste(years, collapse = "|")
# 
# concentration <- list.files("C:/Users/n11222026/OneDrive - Queensland University of Technology/Data/raw/Sea_Ice", pattern = "concentration", full.names = TRUE)
# 
# # Subset the concentration list to only include file paths containing a year from the list
# concentration <- concentration[grepl(pattern, concentration)]
# 
# sea_ice_conc_month <- concentration %>% 
#   map(~rast(.x)) %>% 
#   rast() 
#   
# # Get layer names for November rasters (had to add _ to avoid 202112 year)
# nov <- grep("11_", names(sea_ice_conc_month), value = T)
# 
# # Subset the raster object to keep only the selected layers
# sea_ice_conc_nov <- sea_ice_conc_month[[nov]]
# 
# # Divide by 10 to get percentage
# sea_ice_conc_nov <- sea_ice_conc_nov/ 10
# 
# # Calculate mean sea ice concentration for November over all years
# 
# sea_ice_conc_nov_mean <- app(sea_ice_conc_nov, mean, na.rm=TRUE) 
# 
# sea_ice_conc_nov_mean <- project(sea_ice_conc_nov_mean, "EPSG:3031")
# 
# # Set values of > 100 to NA 
# sea_ice_conc_nov_mean <- ifel(sea_ice_conc_nov_mean > 100, NA, sea_ice_conc_nov_mean)
# 
# 
# 
# # Save the mean sea ice concentration for November
# writeRaster(sea_ice_conc_nov_mean, here("Data/Environmental_predictors/sea_ice_conc_nov_mean.tif"), overwrite = TRUE)
# 
# 
# # Calculate average sea-ice concentration within 100 km of each pixel 
# 
# domain <- rast(here("Data/Environmental_predictors/ice_free_upsamp_1km.tif"))
# 
# domain <- ifel(not.na(domain), 1, NA)
# 
# 
# 
# ########## TO DECIDE:
# # What month to use per species and what percentage ice concentration to use
# 
# ##############################################################################
# 



