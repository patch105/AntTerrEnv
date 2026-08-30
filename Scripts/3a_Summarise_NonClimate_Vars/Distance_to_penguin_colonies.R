

###############################################
# Load libraries -----------------------------------------------------------
################################################

library(terra)
library(here)
library(sf)
library(readr)
library(dplyr)
library(reshape)


# Load domains ------------------------------------------------------------

coast_domain <- rast(here("Data/coast_domain.tif"))
ice_free_domain <- rast(here("Data/ice_free_domain.tif"))


#########################################################
###### Distance to penguin colony - Pygoscelis Penguins ########
#########################################################

# Penguin colony data come from the Antarctic Penguin Biogeographic Project (Che-Castaldo et al., 2023)
# Downloaded at: https://doi.org/10.48361/zftxkr
# Downloads as a zipped DwC-A file with separate 'occurrence' and 'event' text files

penguin_event <- read_tsv(here("Data/dwca-mapppd_count_data-v2.3/event.txt"), quote = "", col_types = cols(.default = "c")) # read everything as character first, then convert as needed

penguin_occurrence <- read_tsv(here("Data/dwca-mapppd_count_data-v2.3/occurrence.txt"), quote = "", col_types = cols(.default = "c"))

penguin_colony <- merge(penguin_event, penguin_occurrence, by = "id", all.x = TRUE)

# Remove emperor & king penguin records, plus macaroni from this dataset
penguin_colony <- penguin_colony %>%
  filter(vernacularName != "emperor penguin") %>%
  filter(vernacularName != "king penguin") %>% 
  filter(vernacularName != "macaroni penguin") # Doing macaroni penguins from updated dataset


# Make dataframe an sf object, then a spatVector
penguin_colony_sf <- penguin_colony %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform(3031)

penguin_colony_sf <- count(penguin_colony_sf, locality, vernacularName)

penguin_colony_SPVE <- vect(penguin_colony_sf)


##################################################
# Load new macaroni penguin breeding records ------------------------------
##################################################

# From Hallet & Lynch (2024)
# Downloaded excel spreadsheet (Supplementary File 1) and saved first page as csv file

macaroni_penguin_colony <- read.csv(here("Data/hallet_lynch_macaroni_breeding_sites.csv"))

# Make dataframe an sf object, then a spatVector
macaroni_penguin_colony_sf <- macaroni_penguin_colony %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326) %>%
  st_transform(3031)

macaroni_penguin_colony_sf <- count(macaroni_penguin_colony_sf, SITE.NAME)

macaroni_penguin_colony_SPVE <- vect(macaroni_penguin_colony_sf)



##################################################
# Combine all penguin colony locs -----------------------------------------
##################################################

penguin_colonies_ALL <- rbind(penguin_colony_SPVE, macaroni_penguin_colony_SPVE)


# Calculate distance to nearest penguin colony.
dist_penguin_colony_coast <- terra::distance(coast_domain, penguin_colonies_ALL)
dist_penguin_colony_coast <- mask(dist_penguin_colony_coast, coast_domain, maskvalue = NA)

dist_penguin_colony_ice_free <- terra::distance(ice_free_domain, penguin_colonies_ALL)
dist_penguin_colony_ice_free <- mask(dist_penguin_colony_ice_free, ice_free_domain, maskvalue = NA)

writeRaster(dist_penguin_colony_coast, here("Outputs/NonClimate_Vars/dist_to_penguin_colony_COASTLINE.tif"), overwrite = T)
writeRaster(dist_penguin_colony_ice_free, here("Outputs/NonClimate_Vars/dist_to_penguin_colony_ICEFREE.tif"), overwrite = T)



###############################################
# Load libraries -----------------------------------------------------------
################################################

library(terra)
library(here)
library(sf)
library(readr)
library(dplyr)
library(reshape)


# Load domains ------------------------------------------------------------

coast_domain <- rast(here("Data/coast_domain.tif"))
ice_free_domain <- rast(here("Data/ice_free_domain.tif"))


#########################################################
###### Distance to penguin colony - Pygoscelis Penguins ########
#########################################################

# Penguin colony data come from the Antarctic Penguin Biogeographic Project (Che-Castaldo et al., 2023)
# Downloaded at: https://doi.org/10.48361/zftxkr
# Downloads as a zipped DwC-A file with separate 'occurrence' and 'event' text files

penguin_event <- read_tsv(here("Data/dwca-mapppd_count_data-v2.3/event.txt"), quote = "", col_types = cols(.default = "c")) # read everything as character first, then convert as needed

penguin_occurrence <- read_tsv(here("Data/dwca-mapppd_count_data-v2.3/occurrence.txt"), quote = "", col_types = cols(.default = "c"))

penguin_colony <- merge(penguin_event, penguin_occurrence, by = "id", all.x = TRUE)

# Update some columns to the right format:
penguin_colony$year <- as.numeric(penguin_colony$year)
penguin_colony$organismQuantity <- as.numeric(penguin_colony$organismQuantity)

# Remove emperor & king penguin records, plus macaroni from this dataset
penguin_colony <- penguin_colony %>%
  filter(vernacularName != "emperor penguin") %>%
  filter(vernacularName != "king penguin") %>% 
  filter(vernacularName != "macaroni penguin") # Doing macaroni penguins from updated dataset


# Make dataframe an sf object, then a spatVector
penguin_colony_sf <- penguin_colony %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform(3031)

penguin_colony_SPVE <- vect(penguin_colony_sf)


# Step 1. Estimate colony radius from breeding pairs -------------------------

# Tidy dataframe
penguin_colony_sf <- penguin_colony_sf %>%
  dplyr::select("locality", "year", "organismQuantity", "organismQuantityType", "lifeStage", "vernacularName") %>%
  mutate(type = ifelse(organismQuantityType == "nests", "nests", lifeStage)) %>%
  dplyr::select(-organismQuantityType, -lifeStage) 

# Load functions to assist with penguin rookery size calculations
# From Aniko Toth: https://github.com/anikobtoth/Antarctica/tree/master/scripts

rm_outliers <- function(data){
  outliers <- boxplot(data, plot = FALSE)$out
  if(length(outliers) == 0){
    return(data)
  } else {
    data_no_outlier <- data[-which(data %in% outliers)]
    return(data_no_outlier)
  }
}

adjust_count <- function(data){
  data %>% log() %>% rm_outliers() %>% exp()
}

BP_translate <- function(nests, adults, chicks, ratios){
  nests[which(nests == 0)] <- NA
  chicks[which(chicks == 0)] <- NA
  adults[which(adults == 0)] <- NA
  
  BP <- nests
  BP <- ifelse(is.na(BP), chicks/ratios["cn"], BP)
  BP <- ifelse(is.na(BP), adults/ratios["an"], BP)
  
  return(BP)
}

namerows <- function(table){
  rownames(table) <- table[,1]
  table <- table[,2:ncol(table)]
  #   return(table)
}


# Load penguin colony size and breeding pair data from LaRue et al. (2014). ("A method for estimating colony sizes of Adelie penguins using remote sensing imagery). 
# * NOTE - LaRue_Penguin_Colonies.csv was made by manually transcribing data from Table 1 of LaRue et al., (2014).

# We fit a linear model to obtain the relationship between breeding pairs and area of the colony. The coefficient and intercept are then input later to estimate the area of the colony based on the number of breeding pairs.

larue <- read.csv(here("Data/LaRue_Penguin_Colonies.csv"))

mod <- lm(area ~ BP,  data = larue)

#summary(mod)

intercept <- mod$coefficients[[1]]

coef <- mod$coefficients[[2]]

# Estimate number of breeding pairs based on nests, or ratio of chick/adults to nests. Use this to estimate the area of the colony in ha.
#
# Code adapted from Aniko Toth: <https://github.com/anikobtoth/Antarctica/blob/master/scripts/overlays.R>

# Get mean of nests, adults, and chicks for each colony
rookeries <- penguin_colony_sf %>%
  group_by(locality, vernacularName) %>%
  summarise(start = min(year), end = max(year),
            nests = mean(organismQuantity[type == "nests"], na.rm = T),
            adults = mean(organismQuantity[type == "adult"], na.rm = T),
            chicks = mean(organismQuantity[type == "chick"], na.rm = T))


rookeries$nests[which(rookeries$nests == 0)] <- NA
rookeries$chicks[which(rookeries$chicks == 0)] <- NA
rookeries$adults[which(rookeries$adults == 0)] <- NA

# compare numbers of adults, nests, and chicks
ratios <- rookeries %>%
  st_drop_geometry() %>%
  mutate(cn = chicks/nests, an = adults/nests, ac = adults/chicks) %>%
  group_by(vernacularName) %>%
  summarise(cn = mean(adjust_count(cn), na.rm = T), # Calling the adjust_count here removes outlier ratios (outside 95% CI)
            an = mean(adjust_count(an), na.rm = T),
            ac = mean(adjust_count(ac), na.rm = T)) %>%
  data.frame() %>%  namerows() %>%  t()

ratios <- list(ratios[,1], ratios[,2], ratios[,3])


rookeriesBP <- rookeries %>%
  split(.$vernacularName) %>%
  map2(ratios, ~.x %>% mutate(BP = BP_translate(nests, adults, chicks, ratios = .y))) %>%  # Translate to Base Pairs
  bind_rows() %>% filter(!is.na(BP)) %>%
  mutate(est_HA =  ((coef*BP) + intercept)/10000,  # This equation is based on the linear regression of breeding pairs vs. colony size
         r   =  2*sqrt(est_HA*10000/pi)) # Calculate the radius of the colony in m, * by 10,000 to convert ha to m^2


# Step 2. Remove spatial outlier colonies ---------------------------------

# We will remove colonies that are \> 3km from the ice-free land layer. This means that they are \> 3km from the 1km grid cell. Those colonies that are not overlapping with ice-free land but are within 4km are snapped to the nearest ice-free land.


# Calculate how far these rookeries are from ice-free land (domain)

# Load the ice-free union layer but at a 100m raster resolution
domain100m <- rast(here("Data/Environmental_predictors/ice_free_union_reproj_100m.tif"))

# Set domain values
# domain100m <- ifel(not.na(domain100m), 1, NA)

# print("domain100m")

# dist_raster <- distance(domain100m)

# print("dist_raster")

dist_raster <- rast(here("Data/Environmental_predictors/Dist_to_Icefree_100m.tif"))




###############################################
# References --------------------------------------------------------------
################################################

# Che-Castaldo, C., Humphries, G., Lynch, H., 2023. Antarctic Penguin Biogeography Project: Database of abundance and distribution for the Adélie, chinstrap, gentoo, emperor, macaroni and king penguin south of 60 S. Biodivers Data J 11, e101476. https://doi.org/10.3897/BDJ.11.e101476

# Hallet, M. and Lynch, H.J., 2024. Update on the abundance and distribution of Macaroni Penguins (Eudyptes chrysolophus) in the Antarctic Peninsula region. Polar Biology, 47(6), pp.607-615.

# LaRue, M.A., Lynch, H.J., Lyver, P.O.B., Barton, K., Ainley, D.G., Pollard, A., Fraser, W.R., Ballard, G., 2014. A method for estimating colony sizes of Adélie penguins using 
# remote sensing imagery. Polar Biol 37, 507–517. https://doi.org/10.1007/s00300-014
# 1451-8 

