
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
# Distance to research facility -------------------------------------------
################################################

# Data from COMNAP (2024)
# NOTE: Downloaded the zipped file and unzipped it manually

facilities_df <- read.csv(here("Data/comnap-antarctic-facilities-3.5.0/dist/csv/COMNAP_Antarctic_Facilities_Master.csv"))
facilities <- st_as_sf(facilities_df,
                     coords = c("Longitude..DD.", "Latitude..DD."),
                     crs = 4326) # Set as WGS 1984

# * NOTE - if you wanted you could optionally filter out for certain facility types (e.g., only stations)
# facilities <- facilities %>% filter(Type == "Station")

facilities <- st_transform(facilities, 3031) # Project to WGS_1984_Stereographic_South_Pole
facilities_SPVE <- vect(facilities)


dist_facility_coast <- terra::distance(coast_domain, facilities_SPVE)
dist_facility_coast <- mask(dist_facility_coast, coast_domain, maskvalue = NA)

dist_facility_ice_free <- terra::distance(ice_free_domain, facilities_SPVE)
dist_facility_ice_free <- mask(dist_facility_ice_free, ice_free_domain, maskvalue = NA)

writeRaster(dist_facility_coast, here("Outputs/NonClimate_Vars/distance_to_facility_COASTLINE.tif"), overwrite = T)
writeRaster(dist_facility_ice_free, here("Outputs/NonClimate_Vars/distance_to_facility_ICEFREE.tif"), overwrite = T)


# # Archived code looking at peak population:
# library(dplyr)
# library(ggplot2)
# 
# # Convert Peak.population to numeric, handling non-numeric characters
# # (commas, "unknown", blanks, etc. will become NA automatically with a warning)
# facilities <- facilities %>%
#   mutate(
#     Peak.population.num = as.numeric(gsub(",", "", Peak.Population)),
#     has_data = !is.na(Peak.population.num)
#   )
# 
# ggplot() +
#   geom_sf(data = facilities, aes(size = Peak.population.num, color = has_data)) +
#   scale_size_continuous(name = "Peak Population", range = c(1, 8)) +
#   scale_color_manual(
#     values = c("TRUE" = "firebrick", "FALSE" = "grey70"),
#     labels = c("TRUE" = "Data available", "FALSE" = "No data"),
#     name = NULL
#   ) +
#   theme_minimal() +
#   labs(title = "Antarctic Facilities Sized by Peak Population")


###############################################
# Distance to Antarctic Specially Protected Areas -----------------------------
################################################

# Data from Terauds et al. (2024)
# NOTE: Downloaded the zipped file and unzipped it manually

ASPAs_SPVE <- vect(here("Data/AAS_4296_Updated_ASPAs_2024/AAS_4296_Updated_ASPAs_2024/ASPAs_polygons_v5_2024/ASPAs_polygons_v5_2024.shp")) #EPSG:3031

# # Dissolve into a single (multi-part) geometry to speed up distance calculation
# ASPAs_agg <- terra::aggregate(ASPAs_SPVE)

dist_ASPA_coast <- terra::distance(coast_domain, ASPAs_SPVE)
dist_ASPA_coast <- mask(dist_ASPA_coast, coast_domain, maskvalue = NA)

dist_ASPA_ice_free <- terra::distance(ice_free_domain, ASPAs_SPVE)
dist_ASPA_ice_free <- mask(dist_ASPA_ice_free, ice_free_domain, maskvalue = NA)

writeRaster(dist_ASPA_coast, here("Outputs/NonClimate_Vars/distance_to_ASPA_COASTLINE.tif"), overwrite = T)
writeRaster(dist_ASPA_ice_free, here("Outputs/NonClimate_Vars/distance_to_ASPA_ICEFREE.tif"), overwrite = T)

# DRAFT STUFF!!!

# Convert non-NA cells to points
ice_free_pts <- as.points(ice_free_domain)

# Now both sides are SpatVector -> use_nodes is valid here
ice_free_pts$dist_ASPA <- terra::distance(ice_free_pts, ASPAs_SPVE, use_nodes = TRUE)

# Rasterize back onto the original grid
dist_ASPA_ice_free <- rasterize(ice_free_pts, ice_free_domain, field = "dist_ASPA")
dist_ASPA_ice_free <- mask(dist_ASPA_ice_free, ice_free_domain, maskvalue = NA)

writeVector(as.points(ASPAs_SPVE), "ASPA_points.shp")



###############################################
# References --------------------------------------------------------------
################################################

# COMNAP. (2024). COMNAP List of Antarctic Facilities v3.5.0. https://github.com/PolarGeospatialCenter/comnap-antarctic-facilities/releases/tag/v3.5.0 [Accessed 1 July 2026]

# Terauds, A., Wauchope, H.S., Wen, W. and Lee, J.R. (2024) Antarctic Specially Protected Areas (Points and Polygons) 2024 Update, Ver. 1, Australian Antarctic Data Centre – https://doi:10.26179/4qk0-cz71 [Accessed 12 September 2025]

