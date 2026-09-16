library(terra)
library(here)
library(dplyr)
library(stringr)
library(tidyr)
library(sf)

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
  ) %>% 
  mutate(zhandian = ifelse(zhandian == "Mt.Fleming", "Mt. Fleming", zhandian),
         zhandian = ifelse(zhandian == "Mt.Erebus", "Mt. Erebus", zhandian))

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
  
  # Austral seasons, conventional definition: December is grouped with
  # Jan/Feb of the FOLLOWING calendar year (e.g. "DJF 1999" = Dec 1998 +
  # Jan/Feb 1999). SeasonYear carries that shifted year for grouping;
  # calendar Year is left untouched for the annual summary below.
  df <- df %>%
    mutate(
      Season = case_when(
        Month %in% c(12, 1, 2)  ~ "DJF",
        Month %in% c(3, 4, 5)   ~ "MAM",
        Month %in% c(6, 7, 8)   ~ "JJA",
        Month %in% c(9, 10, 11) ~ "SON",
        TRUE ~ NA_character_
      ),
      SeasonYear = if_else(Month == 12, Year + 1L, Year)
    )
  
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
    group_by(SeasonYear, Season) %>%
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
    ensure_season_cols() %>%
    rename(Year = SeasonYear)
  
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

max_dist_km          <- 10      # km from ice-free area
min_days_year        <- 183     # min total days/year to consider the year at all
min_days_season      <- 30      # min days in EACH qualifying season for the year to count as seasonally balanced
core_years           <- 1995:2014
extended_years       <- 1985:2014
min_years_core       <- 10      # usable, seasonally-balanced years required within 1995-2014
min_years_extended   <- 14      # usable years required if drawing on the full 1985-2014 window

# For a station-year to count as "usable", require:
#   (a) total ndays >= min_days_year, AND
#   (b) at least 3 of the 4 seasons have >= min_days_season
# (b) is what prevents a summer-heavy year from passing on total days alone,
# while tolerating one weak season (e.g. deep-winter instrument issues).

summarise_qualifying <- function(diagnostics_df, var_prefix) {
  ndays_col    <- paste0(var_prefix, "_ndays")
  season_cols  <- paste0(var_prefix, "_ndays_", c("DJF", "MAM", "JJA", "SON"))
  max_dist_m   <- max_dist_km * 1000
  
  year_level <- diagnostics_df %>%
    filter(dist_to_ice_free <= max_dist_m,
           Year %in% extended_years) %>%
    mutate(
      seasons_met          = rowSums(across(all_of(season_cols), ~ . >= min_days_season)),
      seasonally_balanced  = seasons_met >= 3,   # was: all 4 (seasons_met == 4)
      usable               = .data[[ndays_col]] >= min_days_year & seasonally_balanced,
      period               = if_else(Year %in% core_years, "core", "extension")
    )
  
  station_level <- year_level %>%
    filter(usable) %>%
    group_by(zhandian) %>%
    summarise(
      usable_years          = paste(sort(unique(Year)), collapse = ","),
      n_core_years          = sum(period == "core"),
      n_extended_years      = n(),  # total usable years across full 1985-2014 window
      full_coverage_core    = all(core_years %in% Year[period == "core"]),
      .groups = "drop"
    ) %>%
    mutate(
      qualifies_core     = n_core_years     >= min_years_core,
      qualifies_extended = n_extended_years >= min_years_extended,
      qualifies           = qualifies_core | qualifies_extended,
      tier = case_when(
        full_coverage_core  ~ "Tier 1: full core coverage",
        qualifies_core      ~ "Tier 2: sufficient core years",
        qualifies_extended  ~ "Tier 3: qualifies via extension to 1985",
        TRUE                ~ "Excluded"
      )
    ) %>%
    arrange(match(tier, c("Tier 1: full core coverage",
                          "Tier 2: sufficient core years",
                          "Tier 3: qualifies via extension to 1985",
                          "Excluded")),
            desc(n_core_years))
  
  station_level
}

temp_qual_25 <- summarise_qualifying(diag_25, "temp")
wind_qual_25 <- summarise_qualifying(diag_25, "wind")
temp_qual_75 <- summarise_qualifying(diag_75, "temp")
wind_qual_75 <- summarise_qualifying(diag_75, "wind")

write.csv(temp_qual_25, here("Data/AntAWS/AntAWS_TempQualifying_25.csv"), row.names = FALSE)
write.csv(wind_qual_25, here("Data/AntAWS/AntAWS_WindQualifying_25.csv"), row.names = FALSE)
write.csv(temp_qual_75, here("Data/AntAWS/AntAWS_TempQualifying_75.csv"), row.names = FALSE)
write.csv(wind_qual_75, here("Data/AntAWS/AntAWS_WindQualifying_75.csv"), row.names = FALSE)




AntAWS_locs <- st_read(here("Data/AntAWS/antaws-dataset-x70w9q1u/AWS_location_shapefiles/Shp/267AWS.shp"))
AntAWS_locs <- st_transform(AntAWS_locs, crs = "EPSG:3031")

AntAWS_locs_temp_subset <- AntAWS_locs %>% 
  filter(zhandian %in% temp_stations)

st_write(AntAWS_locs_temp_subset, "Data/AntAWS/AntAWS_locs_temp_subset.shp", append = F)

AntAWS_locs_wind_subset <- AntAWS_locs %>% 
  filter(zhandian %in% wind_stations)

st_write(AntAWS_locs_wind_subset, "Data/AntAWS/AntAWS_locs_wind_subset.shp", append = F)


## =========================================================================
## Step 4: Pull out and save just the usable AntAWS station data
## =========================================================================

temp_stations <- temp_qual_25 %>% filter(qualifies) %>% pull(zhandian)
wind_stations <- wind_qual_25 %>% filter(qualifies) %>% pull(zhandian)

# Re-derive the year-level usable flags (same logic as inside summarise_qualifying)
# so we can filter daily data down to just the qualifying years. Restricted to
# core_years (1995-2014) per the target study period - NOT the 1985 extension,
# even for stations that only qualified as a station via Tier 3.
get_usable_years <- function(diagnostics_df, var_prefix) {
  ndays_col   <- paste0(var_prefix, "_ndays")
  season_cols <- paste0(var_prefix, "_ndays_", c("DJF", "MAM", "JJA", "SON"))
  max_dist_m  <- max_dist_km * 1000
  
  diagnostics_df %>%
    filter(dist_to_ice_free <= max_dist_m,
           Year %in% core_years) %>%
    mutate(
      seasons_met         = rowSums(across(all_of(season_cols), ~ . >= min_days_season)),
      seasonally_balanced = seasons_met >= 3,
      usable              = .data[[ndays_col]] >= min_days_year & seasonally_balanced
    ) %>%
    filter(usable) %>%
    select(zhandian, Year)
}

temp_usable_years <- get_usable_years(diag_25, "temp")
wind_usable_years <- get_usable_years(diag_25, "wind")

# Read and combine daily records for a chosen subset of stations from a folder,
# restricted to that station's qualifying years (via inner_join to usable_years_df),
# keeping only the column relevant to `value_col`, and attaching EPSG:3031
# location + elevation + dist_to_ice_free from stations_df.
extract_daily_subset <- function(folder_path, station_names, value_col, stations_df, usable_years_df) {
  files <- list.files(folder_path, pattern = "_day\\.csv$", full.names = TRUE)
  file_station_names <- str_remove(basename(files), "_day\\.csv$")
  
  keep <- file_station_names %in% station_names
  files <- files[keep]
  file_station_names <- file_station_names[keep]
  
  missing <- setdiff(station_names, file_station_names)
  if (length(missing) > 0) {
    message("Qualifying stations with no matching csv file (", length(missing), "): ",
            paste(missing, collapse = ", "))
  }
  
  results <- vector("list", length(files))
  for (i in seq_along(files)) {
    df <- tryCatch(
      read.csv(files[i], check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(df) || nrow(df) == 0) next
    
    df <- standardize_names(df)
    
    needed <- c("Year", "Month", "Day", value_col)
    for (mc in setdiff(needed, names(df))) df[[mc]] <- NA
    
    results[[i]] <- df %>%
      mutate(Year = suppressWarnings(as.integer(Year))) %>%
      filter(!is.na(Year), !is.na(.data[[value_col]])) %>%
      mutate(zhandian = file_station_names[i]) %>%
      select(zhandian, Year, Month, Day, all_of(value_col))
  }
  
  daily_df <- bind_rows(results)
  
  # Keep only station-years that passed the qualifying-year test
  daily_df <- daily_df %>%
    inner_join(usable_years_df, by = c("zhandian", "Year"))
  
  # Attach station location/attributes (EPSG:3031 coords, elevation, dist_to_ice_free)
  daily_df %>%
    left_join(
      stations_df %>% select(zhandian, x_3031, y_3031, elevation, dist_to_ice_free),
      by = "zhandian"
    )
}

temp_eval_df <- extract_daily_subset(folder_25, temp_stations, "Temperature_C", stations_df, temp_usable_years)
wind_eval_df <- extract_daily_subset(folder_25, wind_stations, "WindSpeed_ms", stations_df, wind_usable_years)

write.csv(temp_eval_df, here("Data/AntAWS/AntAWS_TempEval_25.csv"), row.names = FALSE)
write.csv(wind_eval_df, here("Data/AntAWS/AntAWS_WindEval_25.csv"), row.names = FALSE)


# To account for the uneven spatial distribution of AWSs, daily statistics are area-weighted according to Voronoi polygons around each station. Each station's contribution is proportional to the surface area it represents, preventing that the statistics are biased towards regions where the AWS network is denser, such as Victoria Land.

