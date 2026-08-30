# ==============================================================================
# PolarRES26 -- DAILY VALUES for AntAirICE and AWS evaluations
# ==============================================================================
# Produces one multi-layer daily raster per model x year (reprojected and
# masked to the ice-free domain) for point-based validation against
# AntAirICE / AWS observations.
#
# Uses the same dataset, file-finding, and loading approach as the
# climatology script (find_variable_files / load_variable_series, with the
# RACMO/MetUM grid fix): these netCDFs are already daily-resolution with a
# real time dimension, so there's no sub-daily aggregation step here --
# just load the full series, subset to the target year by date, convert
# units, reproject/mask, and save. No MAR, no manual layer-trimming, no
# filler for missing periods.

# ---- 0. Setup ------------------------------------------------------------------

lib_loc <- paste(getwd(), "/r_lib", sep = "")

library(terra)
library(here)
library(arrow)
library(lubridate)
library(stringr)

# ---- 1. Configuration -----------------------------------------------------------

base_dir <- "/mnt/hpccs01/home/patterc2/n11222026/AntarcticFutureHabitat/Data/PolarRes26"

models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5",
            "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5",
            "MetUM_ERA5")

years_hist <- seq(1994, 2014, by = 1)

# One job = one (model, year) combination.
args <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])

job_grid <- expand.grid(model = models, year = years_hist, stringsAsFactors = FALSE)
model <- job_grid$model[job_index]
year  <- job_grid$year[job_index]

variable <- "tas"

message("Model: ", model, " | Year: ", year, " | Variable: ", variable)

model_dir <- file.path(base_dir, model)
outpath <- here("Data/Environmental_predictors/PolarRes26", model, "Validation")

# ERA5-driven models (HCLIM_ERA5, RACMO_ERA5, MetUM_ERA5) don't have a
# "historical" folder -- their historical-equivalent files live under
# "evaluation" instead. Same logic as the climatology script.
is_era5  <- endsWith(model, "_ERA5")
scenario <- if (is_era5) "evaluation" else "historical"

# RACMO and ERA5-driven grids (i.e. everything that isn't HCLIM) carry
# incorrect/missing CRS + extent metadata in their netCDFs -- fixed in
# place inside load_variable_series(), same as the climatology script.
needs_gridfix <- !startsWith(model, "HCLIM")

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

# ---- 2. Domain (ice-free mask) ---------------------------------------------------

domain <- rast(here("Data/Environmental_predictors/ice_free_upsamp_1km.tif"))
domain <- ifel(not.na(domain), 1, NA)

tmp_dir <- tempdir()
clean_tmp <- function(tmp_dir) {
  tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
  file.remove(tmp_files)
}

# ---- 3. Generic engine: find files, load a continuous daily series -------------
# Same as the climatology script: the variable name must appear as an
# exact path component (so "tas" never matches "tasmax"), and the scenario
# folder must match.

find_variable_files <- function(model_dir, variable_name, scenario) {
  all_files <- list.files(model_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
  if (length(all_files) == 0) return(character(0))
  components_list <- str_split(all_files, "/")
  is_var  <- vapply(components_list, function(cmp) variable_name %in% cmp, logical(1))
  is_scen <- grepl(scenario, all_files, fixed = TRUE)
  sort(all_files[is_var & is_scen])
}

# Load every matching file into ONE time-ordered SpatRaster with real
# dates. DIAGNOSTICS ASSUMED: no spillover, no gaps, no duplicate dates,
# daily, standard calendar.
load_variable_series <- function(model_dir, variable_name, scenario) {
  files <- find_variable_files(model_dir, variable_name, scenario)
  if (length(files) == 0) {
    stop(sprintf("No files found for variable '%s', scenario '%s' under %s",
                 variable_name, scenario, model_dir))
  }
  r <- terra::rast(files, subds = variable_name)
  
  if (needs_gridfix) {
    terra::set.crs(r, crsfix)
    terra::set.ext(r, exfix)
  }
  
  dates <- as.Date(terra::time(r))
  if (length(dates) != terra::nlyr(r) || anyNA(dates)) {
    stop(sprintf("Could not read a clean daily time dimension for '%s' (%s).",
                 variable_name, scenario))
  }
  ord <- order(dates)
  list(r = r[[ord]], dates = dates[ord])
}

# ---- 4. Unit-conversion helper ---------------------------------------------------

to_celsius <- function(r) r - 273.15 # K to celsius

# ---- 5. Save helper ---------------------------------------------------------------

save_daily_raster <- function(r, variable, model, year, outpath) {
  label <- switch(variable, tas = "Temperature", wind = "Wind_Speed")
  out_file <- file.path(outpath, sprintf("Daily_%s_%s_%s_ICEFREE.tif", label, model, year))
  writeRaster(r, out_file, overwrite = TRUE)
  message("Saved: ", out_file)
}

# ==============================================================================
# 6. Load this model's full series and subset to the target year
# ==============================================================================

var_name <- switch(variable, tas = "tas", wind = "sfcWind")
series <- load_variable_series(model_dir, var_name, scenario)

keep_year <- lubridate::year(series$dates) == year
if (!any(keep_year)) stop(paste("No", var_name, "layers found for year", year))

r_daily <- series$r[[keep_year]]
if (variable == "tas") r_daily <- to_celsius(r_daily)

clean_tmp(tmp_dir)

# -------------------------------------------------------------------------
# Step 2: Reproject -> mask to ice-free domain
# -------------------------------------------------------------------------
r_daily <- terra::project(r_daily, domain, method = "near")
r_daily <- mask(r_daily, domain, maskvalue = NA)

# -------------------------------------------------------------------------
# Step 3: Save -- all daily layers for this year in one file
# -------------------------------------------------------------------------
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
save_daily_raster(r_daily, variable, model, year, outpath)