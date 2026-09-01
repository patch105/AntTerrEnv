# ==============================================================================
# PolarRES26 -- DAILY VALUES for model-to-model COMPARISON (historical only)
# ==============================================================================
# Produces one multi-layer daily raster per model x variable (tas, wind),
# covering the full 1995-2014 historical period, on each model's native
# grid (with the RACMO/ERA5-driven CRS+extent fix applied where needed).
#
# Unlike the validation version of this script, there is NO reprojection
# onto a common domain and NO ice-free masking here -- these are meant for
# direct model-to-model comparison on native grids, not point extraction.
#
# Uses the same dataset, file-finding, and loading approach as the
# climatology script (find_variable_files / load_variable_series, with the
# RACMO/MetUM grid fix): these netCDFs are already daily-resolution with a
# real time dimension, so there's no sub-daily aggregation step here --
# just load the full series, subset to the historical years, convert
# units where relevant, and save.

# ---- 0. Setup ------------------------------------------------------------------

lib_loc <- paste(getwd(), "/r_lib", sep = "")

library(terra)
library(here)
library(lubridate)
library(stringr)

# ---- 1. Configuration -----------------------------------------------------------

base_dir <- "/mnt/hpccs01/home/patterc2/n11222026/AntarcticFutureHabitat/Data/PolarRes26"

models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5",
            "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5",
            "MetUM_ERA5")

# Matches the historical period used in the climatology script (Script 2).
years_hist <- seq(1995, 2014, by = 1)
hist_range <- paste(min(years_hist), max(years_hist), sep = "_")

variables <- c("tas", "wind")

# One job = one (model, variable) combination.
args <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])

job_grid <- expand.grid(model = models, variable = variables, stringsAsFactors = FALSE)
model    <- job_grid$model[job_index]
variable <- job_grid$variable[job_index]

message("Model: ", model, " | Variable: ", variable, " | Period: ", hist_range)

model_dir <- file.path(base_dir, model)
outpath <- here("Data/Environmental_predictors/PolarRes26", model, "comparison")

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

tmp_dir <- tempdir()
clean_tmp <- function(tmp_dir) {
  tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
  file.remove(tmp_files)
}

# ---- 2. Generic engine: find files, load a continuous daily series -------------
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

# ---- 3. Unit-conversion helper ---------------------------------------------------

to_celsius <- function(r) r - 273.15 # K to celsius

# ---- 4. Save helper ---------------------------------------------------------------

save_daily_raster <- function(r, variable, model, range, outpath) {
  label <- switch(variable, tas = "Temperature", wind = "Wind_Speed")
  out_file <- file.path(outpath, sprintf("Daily_%s_%s_%s_historical_DAILY.tif", label, model, range))
  writeRaster(r, out_file, overwrite = TRUE)
  message("Saved: ", out_file)
}

# ==============================================================================
# 5. Load this model's full series and subset to the historical period
# ==============================================================================

var_name <- switch(variable, tas = "tas", wind = "sfcWind")
series <- load_variable_series(model_dir, var_name, scenario)

keep_years <- lubridate::year(series$dates) %in% years_hist
if (!any(keep_years)) stop(paste("No", var_name, "layers found for years", hist_range))

r_daily <- series$r[[keep_years]]
if (variable == "tas") r_daily <- to_celsius(r_daily)

clean_tmp(tmp_dir)

# -------------------------------------------------------------------------
# Save -- all daily layers for the historical period in one file, on the
# model's own (fixed, where needed) native grid. No reprojection, no
# domain masking.
# -------------------------------------------------------------------------
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
save_daily_raster(r_daily, variable, model, hist_range, outpath)