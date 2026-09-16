# ==============================================================================
# PolarRes26 -- STEP 4: EXTRACT NEAREST-CELL DAILY VALUES AT AntAWS STATIONS
# ==============================================================================
# Takes the daily comparison rasters written by Script 2 (one multi-layer
# .tif per model x variable, on the model's native grid, unmasked, full
# 1995-2014 historical period) and extracts the single nearest grid cell's
# daily time series at each qualifying AntAWS station location from Script 1
# (AntAWS_locs_temp_subset.shp / AntAWS_locs_wind_subset.shp).
#
# Deliberately does NOT reproject or resample the raster itself -- these
# files are large (thousands of daily layers, several GB each) and we only
# need values at a few dozen point locations. Instead, the station points
# are projected INTO the model's native CRS, and terra::extract() does a
# point-in-cell lookup with a windowed read, touching only the relevant
# pixels across all layers rather than resampling the whole stack. This is
# the cheap side of the round-trip; reprojecting the raster is what Script 3
# does for the climatology products, and is not needed here.
#
# Nearest cell only (method = "simple", terra's default for point extract):
# no fallback search if a station lands on an NA cell. Since these DAILY
# files are explicitly unmasked, an NA here means the station point falls
# outside the model's actual data domain -- a genuine exclusion, not a
# masking artifact to search around.
#
# One job = one (model, variable) combination, same pattern as Scripts 2/3.
# Output: one long-format CSV per job (zhandian x date x model value), plus
# a one-row-per-station diagnostic of distance from the station's true
# EPSG:3031 location to the centre of the nearest cell actually used.
#
# CALENDAR: RACMO_CESM2 specifically uses a 365-day (no-leap) calendar;
# every other job (including the other RACMO driving-GCM combinations,
# HCLIM, and MetUM) uses the standard Gregorian calendar. Rather than trust
# whatever time metadata survives Script 2's GeoTIFF write/read round-trip
# (custom calendars are exactly the kind of thing that can silently
# misdecode via GDAL's generic TIFF time tag), dates are always
# reconstructed explicitly per model, and the layer count is checked
# against the calendar-specific expected day count before trusting the
# reconstruction.
#
# ELEVATION / LAPSE-RATE CORRECTION (tas only): modelled near-surface
# temperature at the nearest grid cell is corrected for the difference
# between the model's own (coarse-resolution) orography at that cell and
# the AWS station's true elevation, using the dry-adiabatic lapse rate
# under a standard-atmosphere assumption:
#   T_corrected = T_model - Gamma_d * (z_station - z_model)
# Requires a static orography ("orog") file per RCM family (HCLIM, RACMO,
# MetUM -- shared across driving GCMs, since orography doesn't depend on
# the driving model), extracted at the same nearest-cell locations. NOT
# needed for wind, so wind jobs run independently of orography being
# available.
#
# NOT done here: joining against the AntAWS observation eval CSVs
# (AntAWS_TempEval_25.csv / AntAWS_WindEval_25.csv) from Script 1. That join
# (by zhandian + Year + Month + Day) is left to a separate combining step
# once every model job has finished, so this script stays parallel/one-job-
# per-file like Scripts 2 and 3.
# ==============================================================================

# ---- 0. Setup ------------------------------------------------------------------

library(terra)
library(here)
library(dplyr)
library(stringr)
library(lubridate)

# ---- 1. Configuration -----------------------------------------------------------

models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5",
            "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5",
            "MetUM_ERA5")

years_hist <- seq(1995, 2014, by = 1)
hist_range <- paste(min(years_hist), max(years_hist), sep = "_")

# variable -> everything that differs between tas and wind jobs:
#   label       matches Script 2's save_daily_raster() file label
#   var_name    kept for reference / sanity messages (matches Script 2)
#   stations_shp the qualifying-station subset written at the end of Script 1
var_config <- list(
  tas  = list(label = "Temperature", var_name = "tas",
              stations_shp = here("Data/AntAWS/AntAWS_locs_temp_subset.shp")),
  wind = list(label = "Wind_Speed",  var_name = "sfcWind",
              stations_shp = here("Data/AntAWS/AntAWS_locs_wind_subset.shp"))
)

input_base  <- here("Data/Environmental_predictors/PolarRes26")
output_base <- here("Data/AntAWS/PolarRes26_Extracted")
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

# Static orography ("orog") files -- one per RCM family, shared across all
# driving GCMs for that RCM (this is why the source URLs all say "ERA5
# evaluation": orography is a fixed field, downloaded once per RCM, not
# once per model job). Not yet downloaded -- expected local paths mirror
# the source THREDDS/JASMIN directory layout. Only read for tas jobs.
orog_paths <- list(
  HCLIM = here("Data/PolarRes26/orog",
               "orog_ANT-12_ERA5_evaluation_r1i1p1f1_HCLIMcom-DMI_HCLIM43-ALADIN_v1-r1_fx.nc"),
  RACMO = here("Data/PolarRes26/orog",
               "orog_ANT-12_ERA5_evaluation_r1i1p1f1_UU-IMAU_RACMO24P-NN_v1-r1_fx.nc"),
  MetUM = here("Data/PolarRes26/orog",
               "orog_ANT-12_ERA5_evaluation_r1i1p1f1_BAS_MetUM_v1-r1_fx.nc")
)

# Dry-adiabatic lapse rate (K/m), standard-atmosphere assumption. Swap for
# 0.0065 if an environmental/standard-atmosphere lapse rate is wanted
# instead of the dry-adiabatic one.
lapse_rate_K_per_m <- 0.0098

# Same HCLIM CRS-recovery template as Scripts 2/3 -- used only as a source
# of CRS metadata, never as data, and only after an extent check.
hclim_crs_template_path <- here(
  "Data/PolarRes26/HCLIM_CESM2/historical/r11i1p1f1/HCLIM43-ALADIN/v1-r1/day/hurs/v20251130/",
  "hurs_ANT-12_CESM2_historical_r11i1p1f1_HCLIMcom-DMI_HCLIM43-ALADIN_v1-r1_day_19860101-19901231.nc"
)

# Same rotated-pole CRS/extent fix as Scripts 2/3, for every non-HCLIM model.
exfix <- ext(c(144, 210, -28.1, 25))
crsfix <- "GEOGCRS[\"Rotated_pole\",
    BASEGEOGCRS[\"unknown\",
        DATUM[\"unnamed\",
            ELLIPSOID[\"Sphere\",6371229,0,
                LENGTHUNIT[\"metre\",1,
                    ID[\"EPSG\",9001]]]],
        PRIMEM[\"Greenwich\",0,
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]]],
    DERIVINGCONVERSION[\"Pole rotation (netCDF CF convention)\",
        METHOD[\"Pole rotation (netCDF CF convention)\"],
        PARAMETER[\"Grid north pole latitude (netCDF CF convention)\",5,
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]],
        PARAMETER[\"Grid north pole longitude (netCDF CF convention)\",20,
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]],
        PARAMETER[\"North pole grid longitude (netCDF CF convention)\",0,
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]]],
    CS[ellipsoidal,2],
        AXIS[\"latitude\",north,
            ORDER[1],
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]],
        AXIS[\"longitude\",east,
            ORDER[2],
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]]]"

# ---- 2. job_index selects ONE (model, variable) combination --------------------

args <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])

job_grid <- expand.grid(model = models, variable = names(var_config),
                        stringsAsFactors = FALSE)
model    <- job_grid$model[job_index]
variable <- job_grid$variable[job_index]
cfg      <- var_config[[variable]]

message("Job ", job_index, "/", nrow(job_grid), " -> Model: ", model,
        " | Variable: ", variable, " | Period: ", hist_range)

input_path <- file.path(input_base, model, "comparison",
                        sprintf("Daily_%s_%s_%s_historical_DAILY.tif",
                                cfg$label, model, hist_range))

if (!file.exists(input_path)) {
  stop("Expected Script 2 output not found: ", input_path)
}

out_csv  <- file.path(output_base, sprintf("%s_%s_station_extract.csv", model, variable))
diag_csv <- file.path(output_base, sprintf("%s_%s_station_diagnostics.csv", model, variable))

# ---- 3. Load the daily raster and confirm/fix its CRS ---------------------------
# Custom rotated-pole WKT (RACMO/MetUM) or a manually-stamped CRS (HCLIM
# quirk) can fail to round-trip cleanly through a GeoTIFF write/read, so
# this repeats the same defensive check Script 3 does, rather than
# assuming Script 2's fix survived.

r <- rast(input_path)
needs_gridfix <- !startsWith(model, "HCLIM")

if (is.na(crs(r)) || crs(r) == "") {
  
  if (needs_gridfix) {
    message("  raster has no CRS -- reapplying rotated-pole CRS/extent fix")
    set.crs(r, crsfix)
    set.ext(r, exfix)
    
  } else {
    hclim_template <- rast(hclim_crs_template_path)
    extents_match <- isTRUE(all.equal(
      as.vector(ext(r)), as.vector(ext(hclim_template)),
      tolerance = 1e-6
    ))
    if (!extents_match) {
      stop("Input raster has no CRS and its extent does not match the HCLIM ",
           "CRS template -- refusing to guess a CRS for a grid that doesn't ",
           "line up. Input: ", input_path,
           "\n  input extent:    ", paste(round(as.vector(ext(r)), 4), collapse = ", "),
           "\n  template extent: ", paste(round(as.vector(ext(hclim_template)), 4), collapse = ", "))
    }
    crs(r) <- crs(hclim_template)
    message("  raster had no CRS -- extent matched the HCLIM template, so CRS ",
            "was copied from it.")
  }
}

# ---- 4. Recover per-layer dates, calendar-aware ----------------------------------
# Script 2 saves one layer per day for the full historical window with no
# gaps/duplicates (same guarantee load_variable_series() relies on). Dates
# are always reconstructed explicitly from the model's known calendar
# (rather than trusting whatever time metadata survived the GeoTIFF
# write/read round-trip), because a no-leap calendar decoded as if it were
# standard Gregorian would silently drift by a day at every leap year.

build_calendar_dates <- function(years, calendar) {
  full <- seq(as.Date(paste0(min(years), "-01-01")),
              as.Date(paste0(max(years), "-12-31")), by = "day")
  if (calendar == "365_day") {
    full <- full[format(full, "%m-%d") != "02-29"]
  }
  full
}

# Only RACMO_CESM2 is currently known to use a no-leap calendar (the other
# RACMO driving-GCM combinations use standard Gregorian) -- everything else
# is assumed standard unless/until another specific model turns out to
# need its own entry here.
calendar <- if (model == "RACMO_CESM2") "365_day" else "standard"
dates <- build_calendar_dates(years_hist, calendar)

if (length(dates) != nlyr(r)) {
  stop("Layer count (", nlyr(r), ") does not match the expected ", calendar,
       "-calendar day count (", length(dates), ") for ", hist_range,
       ". Input: ", input_path,
       " -- check whether this model's calendar assumption is still correct.")
}

message("  calendar: ", calendar, " (", length(dates), " days)")

# ---- 5. Load qualifying station points and reproject them into the model's CRS --
# Reprojecting ~a few dozen points is essentially free; this is what lets us
# avoid reprojecting the multi-GB raster.

stations <- vect(cfg$stations_shp)
stations_native <- project(stations, crs(r))

stopifnot(nrow(stations) == nrow(stations_native))

# ---- 6. Extract the nearest cell's full daily series per station ----------------

ext_vals <- terra::extract(r, stations_native, cells = TRUE, method = "simple")
# ext_vals: ID, cell, <one column per layer>

value_cols <- setdiff(names(ext_vals), c("ID", "cell"))
stopifnot(length(value_cols) == nlyr(r))

zhandian_vec <- as.character(stations$zhandian)

long_df <- ext_vals %>%
  mutate(zhandian = zhandian_vec[ID]) %>%
  select(zhandian, cell, all_of(value_cols)) %>%
  tidyr::pivot_longer(cols = all_of(value_cols), names_to = "layer", values_to = "model_value") %>%
  group_by(zhandian) %>%
  mutate(date = dates[match(layer, value_cols)]) %>%
  ungroup() %>%
  mutate(
    Year  = lubridate::year(date),
    Month = lubridate::month(date),
    Day   = lubridate::day(date),
    model = model,
    variable = variable
  ) %>%
  select(zhandian, Year, Month, Day, model, variable, model_value, cell)

# ---- 6b. tas only: elevation / lapse-rate correction -----------------------------
# T_corrected = T_model - Gamma_d * (z_station - z_model), using the same
# nearest-cell approach (station points projected into the orography
# file's native CRS) rather than resampling anything.

if (variable == "tas") {
  
  rcm_family <- case_when(
    startsWith(model, "HCLIM") ~ "HCLIM",
    startsWith(model, "RACMO") ~ "RACMO",
    startsWith(model, "MetUM") ~ "MetUM",
    TRUE ~ NA_character_
  )
  orog_path <- orog_paths[[rcm_family]]
  
  if (is.null(orog_path) || !file.exists(orog_path)) {
    stop("Orography file not found for RCM family '", rcm_family, "': ",
         orog_path, "\nDownload it (static, shared across driving GCMs for ",
         "this RCM) and place it at that path before running tas jobs. ",
         "Wind jobs do not require this file.")
  }
  
  orog_r <- rast(orog_path)
  
  # Same CRS-fix logic as the main variable raster -- orography for a given
  # RCM sits on that RCM's native grid, so the same needs_gridfix /
  # HCLIM-template handling applies.
  if (is.na(crs(orog_r)) || crs(orog_r) == "") {
    if (needs_gridfix) {
      set.crs(orog_r, crsfix)
      set.ext(orog_r, exfix)
    } else {
      hclim_template <- rast(hclim_crs_template_path)
      extents_match <- isTRUE(all.equal(
        as.vector(ext(orog_r)), as.vector(ext(hclim_template)), tolerance = 1e-6
      ))
      if (!extents_match) {
        stop("Orography raster has no CRS and its extent does not match the ",
             "HCLIM CRS template -- refusing to guess. File: ", orog_path)
      }
      crs(orog_r) <- crs(hclim_template)
    }
  }
  
  orog_stations_native <- project(stations, crs(orog_r))
  orog_vals <- terra::extract(orog_r, orog_stations_native, method = "simple")[, 2]
  
  elev_lookup <- data.frame(
    zhandian = zhandian_vec,
    z_station = stations$elevation,
    z_model = orog_vals
  )
  
  long_df <- long_df %>%
    left_join(elev_lookup, by = "zhandian") %>%
    mutate(
      corrected_model_value = model_value - lapse_rate_K_per_m * (z_station - z_model)
    )
  
  n_missing_orog <- sum(is.na(elev_lookup$z_model))
  if (n_missing_orog > 0) {
    message("  ", n_missing_orog, " station(s) have no orography value at ",
            "their nearest cell (outside domain) -- corrected_model_value ",
            "will be NA for those.")
  }
}

write.csv(long_df, out_csv, row.names = FALSE)
message("  wrote ", out_csv, " (", nrow(long_df), " rows)")

# ---- 7. Per-station diagnostics: distance to the cell actually used ------------
# Cell centre in the model's native CRS, projected back to EPSG:3031 and
# compared against the station's true EPSG:3031 location. Constant per
# station (same cell every day), so this is one row per station, not per day.

cell_xy_native <- xyFromCell(r, ext_vals$cell)
cell_pts_native <- vect(cell_xy_native, crs = crs(r))
# NA cells (station outside the raster's real data domain) can't be
# projected/measured -- keep them as NA distance rather than dropping the
# station from the diagnostics table.
valid <- !is.na(ext_vals$cell)

dist_m <- rep(NA_real_, nrow(ext_vals))
if (any(valid)) {
  cell_pts_3031 <- project(cell_pts_native[valid], "EPSG:3031")
  station_pts_3031 <- project(stations_native[valid], "EPSG:3031")
  dist_m[valid] <- terra::distance(cell_pts_3031, station_pts_3031, pairwise = TRUE)
}

diagnostics <- data.frame(
  zhandian = zhandian_vec,
  model = model,
  variable = variable,
  cell = ext_vals$cell,
  dist_to_nearest_cell_m = dist_m,
  n_na_days = sapply(seq_len(nrow(ext_vals)), function(i) {
    if (is.na(ext_vals$cell[i])) NA_integer_ else sum(is.na(unlist(ext_vals[i, value_cols])))
  })
)

# tas jobs: append the elevation values used for the lapse-rate correction,
# so any station with an implausibly large model/station elevation gap
# (and therefore a large correction) is easy to spot.
if (variable == "tas" && exists("elev_lookup")) {
  diagnostics <- diagnostics %>%
    left_join(elev_lookup, by = "zhandian") %>%
    mutate(elevation_diff_m = z_station - z_model)
}

write.csv(diagnostics, diag_csv, row.names = FALSE)
message("  wrote ", diag_csv)

if (any(is.na(ext_vals$cell))) {
  dropped <- zhandian_vec[is.na(ext_vals$cell)]
  message("  station(s) outside the model's data domain (NA cell), excluded: ",
          paste(dropped, collapse = ", "))
}

message("Done: ", model, " / ", variable)
