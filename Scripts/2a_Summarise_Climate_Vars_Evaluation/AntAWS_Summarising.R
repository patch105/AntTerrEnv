library(terra)
library(here)
library(dplyr)
library(stringr)
library(tidyr)

## =========================================================================
## Step 1: Station locations + distance to ice-free area
## =========================================================================

AntAWS_locs <- vect(here("Data/AntAWS/antaws-dataset-x70w9q1u/AWS_location_shapefiles/Shp/267AWS.shp"))
AntAWS_locs <- project(AntAWS_locs, "EPSG:3031")

# Precomputed distance-to-icefree raster (100m), as in your original script:
# domain100m <- rast(here("Data/ice_free_domain_100m.tif"))
# domain100m <- ifel(!is.na(domain100m), 1, NA)
# dist_raster <- distance(domain100m)
# writeRaster(dist_raster, here("Data/Dist_to_Icefree_100m.tif"))

dist_raster <- rast(here("Data/Dist_to_Icefree_100m.tif"))
AntAWS_locs$dist_to_ice_free <- terra::extract(dist_raster, AntAWS_locs)[, 2]

# Pull attributes into a plain data frame, using the projected (EPSG:3031)
# coordinates from the geometry itself rather than the original lat/lon attribute columns
coords_3031 <- crds(AntAWS_locs)

stations_df <- as.data.frame(AntAWS_locs) %>%
  select(zhandian, elevation, dist_to_ice_free) %>%
  mutate(
    zhandian = as.character(zhandian),
    x_3031 = coords_3031[, 1],
    y_3031 = coords_3031[, 2]
  )

## =========================================================================
## Step 2: Per-station daily data -> temporal coverage summary
## =========================================================================

# Standardise the mangled column names (encoding issue from check.names = FALSE)
standardize_names <- function(df) {
  nms <- names(df)
  nms[str_detect(nms, "^Temperature")] <- "Temperature_C"
  nms[str_detect(nms, "^Wind Speed")]  <- "WindSpeed_ms"
  names(df) <- nms
  df
}

# Ensure all 4 seasonal columns exist per variable even if a season is entirely
# absent from a given station's file (pivot_wider only creates columns it sees data for)
ensure_season_cols <- function(df) {
  needed <- c("temp_ndays_DJF", "temp_ndays_MAM", "temp_ndays_JJA", "temp_ndays_SON",
              "wind_ndays_DJF", "wind_ndays_MAM", "wind_ndays_JJA", "wind_ndays_SON")
  for (col in setdiff(needed, names(df))) df[[col]] <- 0
  df
}

# For one station's csv: coverage per year (1995-2014) for temp and wind,
# plus a seasonal (DJF/MAM/JJA/SON) day-count breakdown for each variable
summarise_coverage <- function(file_path, station_name) {
  df <- tryCatch(
    read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df <- standardize_names(df)
  
  needed <- c("Year", "Month", "Day", "Temperature_C", "WindSpeed_ms")
  for (mc in setdiff(needed, names(df))) df[[mc]] <- NA
  
  df <- df %>%
    mutate(Year = suppressWarnings(as.integer(Year))) %>%
    filter(!is.na(Year), Year >= 1995, Year <= 2014)
  
  if (nrow(df) == 0) return(NULL)
  
  # Austral seasons. NOTE: this groups December with Jan/Feb of the SAME
  # calendar year (e.g. "DJF 1998" = Dec 1998 + Jan/Feb 1998), not the
  # conventional austral-summer definition of Dec(year) + Jan/Feb(year+1).
  # Adjust here if you need the cross-year-boundary version instead.
  df <- df %>%
    mutate(Season = case_when(
      Month %in% c(12, 1, 2)  ~ "DJF",
      Month %in% c(3, 4, 5)   ~ "MAM",
      Month %in% c(6, 7, 8)   ~ "JJA",
      Month %in% c(9, 10, 11) ~ "SON",
      TRUE ~ NA_character_
    ))
  
  annual <- df %>%
    group_by(Year) %>%
    summarise(
      temp_months = paste(sort(unique(Month[!is.na(Temperature_C)])), collapse = ","),
      temp_ndays  = sum(!is.na(Temperature_C)),
      wind_months = paste(sort(unique(Month[!is.na(WindSpeed_ms)])), collapse = ","),
      wind_ndays  = sum(!is.na(WindSpeed_ms)),
      .groups = "drop"
    )
  
  seasonal <- df %>%
    filter(!is.na(Season)) %>%
    group_by(Year, Season) %>%
    summarise(
      temp_ndays = sum(!is.na(Temperature_C)),
      wind_ndays = sum(!is.na(WindSpeed_ms)),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from  = Season,
      values_from = c(temp_ndays, wind_ndays),
      values_fill = 0
    ) %>%
    ensure_season_cols()
  
  annual %>%
    left_join(seasonal, by = "Year") %>%
    mutate(across(matches("^(temp|wind)_ndays_(DJF|MAM|JJA|SON)$"), ~ replace(., is.na(.), 0))) %>%
    filter(temp_ndays > 0 | wind_ndays > 0) %>%
    mutate(zhandian = station_name) %>%
    select(zhandian, Year, temp_months, temp_ndays, wind_months, wind_ndays,
           starts_with("temp_ndays_"), starts_with("wind_ndays_"))
}

# Build the full diagnostics table for one folder (25% or 75%)
build_diagnostics <- function(folder_path, stations_df) {
  files <- list.files(folder_path, pattern = "_day\\.csv$", full.names = TRUE)
  file_station_names <- str_remove(basename(files), "_day\\.csv$")
  
  results <- vector("list", length(files))
  names(results) <- file_station_names
  for (i in seq_along(files)) {
    results[[i]] <- summarise_coverage(files[i], file_station_names[i])
  }
  coverage_df <- bind_rows(results)
  
  # Flag mismatches between shapefile stations and csv files
  no_csv    <- setdiff(stations_df$zhandian, file_station_names)
  no_shape  <- setdiff(file_station_names, stations_df$zhandian)
  if (length(no_csv) > 0) {
    message("Stations in shapefile with no matching csv (", length(no_csv), "): ",
            paste(no_csv, collapse = ", "))
  }
  if (length(no_shape) > 0) {
    message("CSV files with no matching station in shapefile (", length(no_shape), "): ",
            paste(no_shape, collapse = ", "))
  }
  
  stations_df %>%
    inner_join(coverage_df, by = "zhandian") %>%
    arrange(zhandian, Year)
}

folder_25 <- here("Data/AntAWS/antaws-dataset-x70w9q1u/AntAWSvers2/The AntAWS dataset/Daily_25%/Daily_25%/")
folder_75 <- here("Data/AntAWS/antaws-dataset-x70w9q1u/AntAWSvers2/The AntAWS dataset/Daily_75%/Daily_75%/")

diag_25 <- build_diagnostics(folder_25, stations_df)
diag_75 <- build_diagnostics(folder_75, stations_df)

write.csv(diag_25, here("Data/AntAWS/AntAWS_Diagnostics_25.csv"), row.names = FALSE)
write.csv(diag_75, here("Data/AntAWS/AntAWS_Diagnostics_75.csv"), row.names = FALSE)

## =========================================================================
## Step 3: Stations meeting distance + coverage criteria
## =========================================================================

max_dist_km   <- 10      # <-- update this threshold (km from ice-free area) as needed
min_days      <- 100    # minimum days of data within a year to count that year
full_years    <- 1995:2014
min_years     <- 10

# For a given diagnostics table (diag_25 or diag_75) and variable ("temp" or "wind"):
# flags whether each station has full 1995-2014 coverage (every year with >= min_days),
# and whether it has at least ~60% of the study period's years (anywhere in the
# range, not necessarily contiguous) with >= min_days
summarise_qualifying <- function(diagnostics_df, var_prefix) {
  ndays_col  <- paste0(var_prefix, "_ndays")
  max_dist_m <- max_dist_km * 1000   # dist_to_ice_free is in metres
  
  diagnostics_df %>%
    filter(dist_to_ice_free <= max_dist_m,
           .data[[ndays_col]] >= min_days) %>%
    group_by(zhandian) %>%
    summarise(
      n_qualifying_years = n(),
      qualifying_years   = paste(sort(unique(Year)), collapse = ","),
      full_coverage_1995_2014 = all(full_years %in% Year),
      min_years_met_60pct     = n_qualifying_years >= min_years,
      .groups = "drop"
    ) %>%
    arrange(desc(full_coverage_1995_2014), desc(n_qualifying_years))
}

temp_qual_25 <- summarise_qualifying(diag_25, "temp")
wind_qual_25 <- summarise_qualifying(diag_25, "wind")
temp_qual_75 <- summarise_qualifying(diag_75, "temp")
wind_qual_75 <- summarise_qualifying(diag_75, "wind")


write.csv(temp_qual_25, here("Data/AntAWS/AntAWS_TempQualifying_25.csv"), row.names = FALSE)
write.csv(wind_qual_25, here("Data/AntAWS/AntAWS_WindQualifying_25.csv"), row.names = FALSE)
write.csv(temp_qual_75, here("Data/AntAWS/AntAWS_TempQualifying_75.csv"), row.names = FALSE)
write.csv(wind_qual_75, here("Data/AntAWS/AntAWS_WindQualifying_75.csv"), row.names = FALSE)

# Names of stations meeting each criterion, e.g.:
temp_stations <- temp_qual_25 %>% filter(min_years_met_60pct) %>% pull(zhandian)
wind_stations <- wind_qual_25 %>% filter(min_years_met_60pct) %>% pull(zhandian)

AntAWS_locs <- st_read(here("Data/AntAWS/antaws-dataset-x70w9q1u/AWS_location_shapefiles/Shp/267AWS.shp"))
AntAWS_locs <- st_transform(AntAWS_locs, crs = "EPSG:3031")

AntAWS_locs_temp_subset <- AntAWS_locs %>% 
  filter(zhandian %in% temp_stations)

st_write(AntAWS_locs_temp_subset, "Data/AntAWS/AntAWS_locs_temp_subset.shp", append = F)

AntAWS_locs_wind_subset <- AntAWS_locs %>% 
  filter(zhandian %in% wind_stations)

st_write(AntAWS_locs_wind_subset, "Data/AntAWS/AntAWS_locs_wind_subset.shp", append = F)


# TO DO IF I GO AHEAD WITH IT: --------------------------------------------

# DO THE ELEVATION CORRECTION BASED ON ADABATIC LAPSE RATE AND ELEVATION FROM HCLIM, RACMO
# FIX THE WIND SOMEHOW AS WELL
# VORONOI EVALUATION FOR FINAL METRIC 




