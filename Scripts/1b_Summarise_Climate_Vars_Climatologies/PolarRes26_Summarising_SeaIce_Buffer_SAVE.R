# ==============================================================================
# PolarRes26 -- STEP 2: SUMMARISE INTO CLIMATOLOGIES
# ==============================================================================

# ---- 0. Setup ------------------------------------------------------------------

lib_loc <- paste(getwd(), "/r_lib_new", sep = "")
# .libPaths(lib_loc)  

library(terra)
library(here)
library(lubridate)
library(stringr)
library(dplyr, lib.loc = lib_loc)
library(purrr, lib.loc = lib_loc)
library(ncdf4, lib.loc = lib_loc)


# ---- 1. Configuration -----------------------------------------------------------

base_dir <- "/mnt/hpccs01/home/patterc2/n11222026/AntarcticFutureHabitat/Data/PolarRes26"

# models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5", "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5", "MetUM_ERA5")
models <- c("RACMO_CESM2", "RACMO_MPI_ESM1")

ssp_scenario <- "ssp370"

years_hist   <- seq(1995, 2014, by = 1)
years_mid    <- seq(2041, 2060, by = 1)
years_future <- seq(2081, 2100, by = 1)

# One job = one (model, variable) combination, matching the two-argument
# job-array style. args[1] = model index (1-4), args[2] = variable index.
args <- commandArgs(trailingOnly = TRUE)
job_index    <- as.integer(args[1])

variables <- list("temp", "total_DD", "wind", "sea_ice",
                  "total_precip", "total_summer_precip", "mean_precip",
                  "mean_summer_precip", "solar_rad", "mean_melt", "total_melt",
                  "mean_snow", "mean_hurs", "vpd", "temp_min", "temp_max")
# temp_minmax added at the END (not alphabetically/logically placed) so it
# doesn't shift the job_index of every variable after it in job_table.



job_table <- expand.grid(model = models, variable = variables,
                         stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

model    <- job_table$model[job_index]
variable <- job_table$variable[job_index]
message("Job ", job_index, " -> Model: ", model, " | Variable: ", variable)

model_dir <- file.path(base_dir, model)

outpath <- here("Data/Environmental_predictors/PolarRes26", model)
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

# RACMO and ERA5-driven grids (i.e. everything that isn't HCLIM) carry
# incorrect/missing CRS + extent metadata in their netCDFs. `needs_gridfix`
# flags every non-HCLIM model, and `crsfix` / `exfix` are the corrected
# rotated-pole CRS and extent applied to every raster read for those models
# -- see the fix applied inside load_variable_series() in section 2, right
# after the raster is read in and before anything else touches it.
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

# RACMO's ssp series currently ends 2099-12-31 (one year short of HCLIM's
# 2100). For RACMO, the years actually used
# to compute the FUTURE period are clipped to 2099, but the output
# filenames keep the original "2081_2100" label (set below) so RACMO and
# HCLIM outputs stay directly comparable.
racmo_future_end_year <- 2099
years_future_for_model <- if (startsWith(model, "RACMO")) {
  years_future[years_future <= racmo_future_end_year]
} else {
  years_future
}

# The three climatology periods, each tagged with which scenario folder its
# source files live in. `range` always uses the ORIGINAL years_hist/mid/
# future vectors for the filename label, even when `years` (what's actually
# used to compute the climatology) has been clipped for a particular model.
#
# ERA5-driven models (HCLIM_ERA5, RACMO_ERA5, MetUM_ERA5) don't have an ssp
# series at all, and their historical-equivalent files live under an
# "evaluation" folder rather than "historical". So for these models we only
# build a HISTORICAL period, pointed at scenario "evaluation" -- the period
# label and range are kept identical to the other models' HISTORICAL entry
# so output filenames still read "..._HISTORICAL_1995_2014.tif" the same way.
is_era5 <- endsWith(model, "_ERA5")

periods <- if (is_era5) {
  list(
    HISTORICAL = list(years = years_hist, scenario = "evaluation",
                      range = paste(min(years_hist), max(years_hist), sep = "_"))
  )
} else {
  list(
    HISTORICAL = list(years = years_hist,   scenario = "historical",
                      range = paste(min(years_hist), max(years_hist), sep = "_")),
    MID        = list(years = years_mid,    scenario = ssp_scenario,
                      range = paste(min(years_mid), max(years_mid), sep = "_")),
    FUTURE     = list(years = years_future_for_model, scenario = ssp_scenario,
                      range = paste(min(years_future), max(years_future), sep = "_"))
  )
}

# ---- 2. Generic engine: find files, load a continuous daily series -------------

# Folder-based match, the variable name must appear as an exact path component (so "tas" never matches
# "tasmax"), and the scenario folder must match.

# This function finds all the files in the folder, matches the ones that are the right variable and scenario (historical or ssp) and just returns the matching ones
find_variable_files <- function(model_dir, variable_name, scenario) {
  
  all_files <- list.files(model_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
  
  if (length(all_files) == 0) return(character(0))
  components_list <- str_split(all_files, "/")
  
  is_var   <- vapply(components_list, function(cmp) variable_name %in% cmp, logical(1))
  
  is_scen  <- grepl(scenario, all_files, fixed = TRUE)
  
  sort(all_files[is_var & is_scen])
}

# Load every matching file into ONE time-ordered SpatRaster with real dates.
# DIAGNOSTICS ASSUMED: no spillover, no gaps, no duplicate dates, daily,
# standard calendar. If that's not true for some variable, this is the
# function to adjust (e.g. add a step to drop a known spillover day).
load_variable_series <- function(model_dir, variable_name, scenario) {
  files <- find_variable_files(model_dir, variable_name, scenario)
  if (length(files) == 0) {
    stop(sprintf("No files found for variable '%s', scenario '%s' under %s",
                 variable_name, scenario, model_dir))
  }
  r <- terra::rast(files, subds = variable_name)
  
  # RACMO / ERA5-driven grids: fix the CRS + extent metadata in place before
  # anything else (time parsing, unit conversion, etc.) touches this raster.
  if (needs_gridfix) {
    terra::set.crs(r, crsfix)
    terra::set.ext(r, exfix)
  }
  
  dates <- as.Date(terra::time(r))
  if (length(dates) != terra::nlyr(r) || anyNA(dates)) {
    stop(sprintf("Could not read a clean daily time dimension for '%s' (%s). ",
                 variable_name, scenario),
         "Re-check Script 1's diagnostics for this variable before trusting this run.")
  }
  ord <- order(dates)
  list(r = r[[ord]], dates = dates[ord])
}

# ---- 3. Generic engine: climatology maths ---------------------------------------

# Per-year mean for one (month, year_offset) combination, averaged across years.
# year_offset = -1 means "December belonging to the summer of `year`", i.e.
# pull December from year-1.
# NOTE: this is the plain/generic version used by temp, wind, precip, solar
# rad, melt, snow, hurs and vpd. Sea ice needs extra land-masking/South
# Georgia steps in between, so it has its own separate function instead of
# bolting that on here via a callback - see climatological_monthly_mean_sea_ice() in section 6.
climatological_monthly_mean <- function(r, dates, month, years, year_offset = 0) {
  yearly <- list()
  for (y in years) {
    target_year <- y + year_offset
    idx <- which(lubridate::month(dates) == month & lubridate::year(dates) == target_year)
    if (length(idx) == 0) {
      warning(sprintf("No data for month %d, year %d (climatology year %d) -- skipping",
                      month, target_year, y))
      next
    }
    sub <- r[[idx]]
    yearly[[as.character(y)]] <- app(sub, mean, na.rm = TRUE)
  }
  if (length(yearly) == 0) stop("No years had any data -- check the date range and scenario coverage.")
  app(rast(yearly), mean, na.rm = TRUE)
}

# All 12 climatological monthly means (offset 0, i.e. no year-crossing) --
# used for the "mean annual X" variables, which also save each month.
climatological_monthly_means_all <- function(r, dates, years) {
  out <- vector("list", 12)
  names(out) <- month.name
  for (m in 1:12) {
    out[[m]] <- climatological_monthly_mean(r, dates, m, years, year_offset = 0)
  }
  out
}

annual_mean_from_monthly <- function(monthly_list) {
  app(rast(monthly_list), mean, na.rm = TRUE)
}
# CHECKED: every "mean annual X" variable below (temp, wind, mean_precip,
# solar_rad, mean_melt, mean_snow, mean_hurs, vpd) builds its annual figure
# by calling this on its 12 climatological_monthly_mean() results -- i.e.
# the annual value is just the mean of the 12 monthly climatologies, not a
# separately-computed daily/annual mean. Nothing to change there.

# Sum (not mean) of a 12-month climatology list -- used for the "total"
# variables (total_DD, total_precip, total_melt) so their annual figure is
# a true Jan-Dec sum of the monthly climatological totals, rather than an
# average of them.
annual_total_from_monthly <- function(monthly_list) {
  app(rast(monthly_list), sum, na.rm = TRUE)
}

# Mean of climatological_monthly_mean() across a set of (month, offset) pairs.
climatological_multi_month_mean <- function(r, dates, months, offsets, years) {
  stopifnot(length(months) == length(offsets))
  parts <- Map(function(m, o) climatological_monthly_mean(r, dates, m, years, o),
               months, offsets)
  app(rast(parts), mean, na.rm = TRUE)
}

# DJF (Dec-Jan-Feb) austral summer: December comes from the year BEFORE.
summer_mean_climatology <- function(r, dates, years) {
  climatological_multi_month_mean(r, dates, c(12, 1, 2), c(-1, 0, 0), years)
}

# JJA austral winter -- no year-crossing, so this can reuse an already-
# computed monthly list (see the temp block below) instead of recomputing.
winter_mean_from_monthly <- function(monthly_list) {
  app(rast(monthly_list[c("June", "July", "August")]), mean, na.rm = TRUE)
}

# Sum-style: per-year (or per-season) total over an arbitrary set of
# (month, offset) windows, then averaged across years. `transform` is
# applied to the daily values before summing (e.g. a flux-to-mm conversion).
period_total_climatology <- function(r, dates, years, months, offsets, transform = NULL) {
  yearly <- list()
  for (y in years) {
    idx_all <- integer(0)
    for (i in seq_along(months)) {
      target_year <- y + offsets[i]
      idx_all <- c(idx_all, which(lubridate::month(dates) == months[i] & lubridate::year(dates) == target_year))
    }
    if (length(idx_all) == 0) { warning(sprintf("No data for year %d -- skipping", y)); next }
    sub <- r[[sort(idx_all)]]
    if (!is.null(transform)) sub <- transform(sub)
    yearly[[as.character(y)]] <- app(sub, sum, na.rm = TRUE)
  }
  if (length(yearly) == 0) stop("No years had any data -- check the date range and scenario coverage.")
  app(rast(yearly), mean, na.rm = TRUE)
}

annual_total_climatology <- function(r, dates, years, transform = NULL) {
  period_total_climatology(r, dates, years, 1:12, rep(0, 12), transform)
}
# NOTE: annual_total_climatology() above is now only used as the building
# block for summer_total_climatology() below (via period_total_climatology).
# The plain annual "total" variables (total_DD, total_precip, total_melt)
# no longer call it directly -- they now go through
# climatological_monthly_sum(s)_all() instead, so they save monthly
# climatology outputs the same way the "mean annual" variables do (see
# point 2 of the brief: daily -> monthly sum -> monthly climatology ->
# saved -> annual figure = mean of those 12 monthlies).

summer_total_climatology <- function(r, dates, years, transform = NULL) {
  period_total_climatology(r, dates, years, c(12, 1, 2), c(-1, 0, 0), transform)
}
# total_summer_precip (section 5b below) still uses this as-is -- left
# unchanged for now per the brief, but flagged there as an inconsistency:
# it's now the only "total" variable that doesn't produce monthly
# climatology outputs the way the others below do.

# Per-year monthly TOTAL for one (month, year_offset) combination, averaged
# across years -- the "total" counterpart to climatological_monthly_mean()
# above. This is what lets the accumulation/total-style variables
# (degree-days, melt, precip) be built the same way as the mean-style ones:
# 12 climatological monthly figures (here, monthly SUMS rather than monthly
# MEANS), saved individually, then combined with annual_mean_from_monthly()
# for the annual figure.
climatological_monthly_sum <- function(r, dates, month, years, year_offset = 0, transform = NULL) {
  yearly <- list()
  for (y in years) {
    target_year <- y + year_offset
    idx <- which(lubridate::month(dates) == month & lubridate::year(dates) == target_year)
    if (length(idx) == 0) {
      warning(sprintf("No data for month %d, year %d (climatology year %d) -- skipping",
                      month, target_year, y))
      next
    }
    sub <- r[[idx]]
    if (!is.null(transform)) sub <- transform(sub)
    yearly[[as.character(y)]] <- app(sub, sum, na.rm = TRUE)
  }
  if (length(yearly) == 0) stop("No years had any data -- check the date range and scenario coverage.")
  app(rast(yearly), mean, na.rm = TRUE)
}

# All 12 climatological monthly TOTALS (offset 0, no year-crossing) --
# mirrors climatological_monthly_means_all() but summing within each
# month instead of averaging.
climatological_monthly_sums_all <- function(r, dates, years, transform = NULL) {
  out <- vector("list", 12)
  names(out) <- month.name
  for (m in 1:12) {
    out[[m]] <- climatological_monthly_sum(r, dates, m, years, year_offset = 0, transform = transform)
  }
  out
}

# ---- 4. Unit-conversion helpers --------------------------------------------------
# DIAGNOSTICS ASSUMED: units match what Script 1's variable_summary CSV said
# matched the old scripts' expectations (K for temp, kg m-2 s-1 flux for pr).

to_celsius <- function(r) r - 273.15 # K to celcius
flux_to_mm <- function(r) r * 86400 # Daily flux to daily mm by multiplying by seconds

# ---- 5. Output queue -------------------------------------------------------------
# Every result is queued here and only actually written in the "FINALISE"
# section at the bottom of the script, once every product for this
# variable/model/period has been computed. RACMO and ERA5-driven models are
# written out in their own (fixed) native grid -- see section 1 for the
# CRS/extent fix applied to every raster as it's read in -- there is no
# longer any resampling onto a common HCLIM grid at the end.

pending_outputs <- list()

queue_raster <- function(r, filename) {
  pending_outputs[[length(pending_outputs) + 1]] <<- list(r = r, filename = filename)
  invisible(NULL)
}

# The actual disk write (only called from the FINALISE section at the end).
save_raster <- function(r, filename) {
  path <- file.path(outpath, filename)
  writeRaster(r, path, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  message("  wrote ", filename)
}

# ---- 6. Sea-ice: dedicated climatology function + spatial helpers --------------
# This does the same job as climatological_monthly_mean() in section 3, but with land-masking

# Per-year mean sea ice concentration for one (month, year_offset)
# combination, masking out land and correcting the South Georgia / South
# America artefact each year, then averaged across years.
climatological_monthly_mean_sea_ice <- function(r, dates, month, years, year_offset, coast, SG) {
  yearly <- list()
  for (y in years) {
    target_year <- y + year_offset
    idx <- which(lubridate::month(dates) == month & lubridate::year(dates) == target_year)
    if (length(idx) == 0) {
      warning(sprintf("No sea-ice data for month %d, year %d (climatology year %d) -- skipping",
                      month, target_year, y))
      next
    }
    r_month <- r[[idx]]
    
    # Mask out land/coast (keep everything OUTSIDE the coastline polygon)
    r_month <- mask(r_month, coast, inverse = TRUE, touches = FALSE)
    
    # Zero-out the South Georgia / South America cells that sit at 100% SIC
    # year-round -- this is a model-domain artefact, not real sea ice.
    IDs <- terra::extract(r_month, SG, cells = TRUE)
    IDs <- IDs[!is.na(IDs[, 2]) & IDs[, 2] == 100, ]
    if (nrow(IDs) > 0) r_month[IDs$cell] <- 0
    
    yearly[[as.character(y)]] <- app(r_month, mean, na.rm = TRUE)
  }
  if (length(yearly) == 0) stop("No years had any sea-ice data -- check the date range and scenario coverage.")
  app(rast(yearly), mean, na.rm = TRUE)
}

# Extract a concentration raster's mean value within a pre-built buffer
# around every ice-free-land cell, and place it back onto the domain grid.
extract_to_buffer <- function(conc_raster, domain, domain.pts, buffer_vect) {
  conc_raster <- ifel(is.na(conc_raster), 0, conc_raster)
  extracted <- terra::extract(conc_raster, buffer_vect, fun = mean, na.rm = TRUE)
  extracted[is.na(extracted)] <- 0
  out <- domain
  values(out) <- NA
  cell_ids <- cellFromXY(out, crds(domain.pts))
  out[cell_ids] <- extracted[, 2]
  out
}


# ==============================================================================
# 7. Create the variables 
# ==============================================================================

# ---- 1. TEMPERATURE (tas, K -> degC): annual, summer (DJF), winter (JJA) -------
if (variable == "temp") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- temp: ", period_name)
    
    series <- load_variable_series(model_dir, "tas", p$scenario) # result of this function is a list of rasters and their dates
    r_celsius <- to_celsius(series$r)
    
    # Annual: 12 monthly climatologies, saved individually, then averaged
    # Inputs are: all the rasters, their dates, and all the years in that period
    monthly <- climatological_monthly_means_all(r_celsius, series$dates, p$years)
    for (m in 1:12) {
      queue_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Temperature_%s_%s_%s.tif",
                                         month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    queue_raster(annual, sprintf("Mean_Annual_Temperature_%s_%s.tif", period_name, p$range))
    
    # Summer (DJF) -- needs its own December-from-previous-year sweep
    summer <- summer_mean_climatology(r_celsius, series$dates, p$years)
    queue_raster(summer, sprintf("Mean_Summer_Temperature_%s_%s.tif", period_name, p$range))
    
    # Winter (JJA) -- no year-crossing, reuse the annual monthly list
    winter <- winter_mean_from_monthly(monthly)
    queue_raster(winter, sprintf("Mean_Winter_Temperature_%s_%s.tif", period_name, p$range))
  }
}

# ---- 1b. TASMIN  (K -> degC): monthly climatology only ------------------
# Added for the separate bioclim-summarising script, which needs monthly
# min/max temperature climatologies (not the "mean" tas used by the "temp"
# variable above). No annual/summer/winter aggregate computed here since
# bioclim only needs the 12 monthly values each -- add one if useful
# elsewhere later.
if (variable == "temp_min") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- temp_min: ", period_name)
    
    tasmin_series <- load_variable_series(model_dir, "tasmin", p$scenario)
    
    tasmin_celsius <- to_celsius(tasmin_series$r)
    
    monthly_min <- climatological_monthly_means_all(tasmin_celsius, tasmin_series$dates, p$years)
    for (m in 1:12) {
      queue_raster(monthly_min[[m]], sprintf("Climatological_Monthly_Mean_TasMin_%s_%s_%s.tif",
                                             month.name[m], period_name, p$range))
    }
    
  }
}

# ---- 1c. TASMAX (K -> degC): monthly climatology only ------------------
# Added for the separate bioclim-summarising script, which needs monthly
# min/max temperature climatologies (not the "mean" tas used by the "temp"
# variable above). No annual/summer/winter aggregate computed here since
# bioclim only needs the 12 monthly values each -- add one if useful
# elsewhere later.
if (variable == "temp_max") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- temp_max: ", period_name)
    
    tasmax_series <- load_variable_series(model_dir, "tasmax", p$scenario)
    
    tasmax_celsius <- to_celsius(tasmax_series$r)
    
    monthly_min <- climatological_monthly_means_all(tasmax_celsius, tasmax_series$dates, p$years)
    for (m in 1:12) {
      queue_raster(monthly_min[[m]], sprintf("Climatological_Monthly_Mean_TasMax_%s_%s_%s.tif",
                                             month.name[m], period_name, p$range))
    }
    
  }
}

# ---- 2. TOTAL POSITIVE DEGREE-DAYS, two thresholds: -5 degC and 0 degC (tas) ----
if (variable == "total_DD") {
  
  # Both thresholds use the same "excess above the limit, floored at 0" rule
  dd_transform <- function(r_daily_celsius, limit) {
    ifel(r_daily_celsius > limit, r_daily_celsius - limit, 0)
  }
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- total_DD: ", period_name)
    
    series <- load_variable_series(model_dir, "tas", p$scenario)
    r_celsius <- to_celsius(series$r)
    
    # -5 degC threshold: daily -> monthly sum -> monthly climatology
    # (saved individually, same as the "mean annual X" variables), then the
    # annual figure is the mean of those 12 monthly climatologies
    monthly_minus5 <- climatological_monthly_sums_all(r_celsius, series$dates, p$years,
                                                      transform = function(r) dd_transform(r, -5))
    for (m in 1:12) {
      queue_raster(monthly_minus5[[m]], sprintf("Climatological_Monthly_Total_Degree_Days-5_%s_%s_%s.tif",
                                                month.name[m], period_name, p$range))
    }
    dd_minus5 <- annual_total_from_monthly(monthly_minus5)
    queue_raster(dd_minus5, sprintf("Mean_Total_Annual_Degree_Days-5_%s_%s.tif", period_name, p$range))
    
    # 0 degC threshold, same pattern
    monthly_0 <- climatological_monthly_sums_all(r_celsius, series$dates, p$years,
                                                 transform = function(r) dd_transform(r, 0))
    for (m in 1:12) {
      queue_raster(monthly_0[[m]], sprintf("Climatological_Monthly_Total_Degree_Days0_%s_%s_%s.tif",
                                           month.name[m], period_name, p$range))
    }
    dd_0 <- annual_total_from_monthly(monthly_0)
    queue_raster(dd_0, sprintf("Mean_Total_Annual_Degree_Days0_%s_%s.tif", period_name, p$range))
  }
}

# ---- 3. WIND SPEED (sfcWind, m s-1, no conversion): annual ----------------------
if (variable == "wind") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- wind: ", period_name)
    
    series <- load_variable_series(model_dir, "sfcWind", p$scenario)
    
    monthly <- climatological_monthly_means_all(series$r, series$dates, p$years)
    for (m in 1:12) {
      queue_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Wind_Speed_%s_%s_%s.tif",
                                         month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    queue_raster(annual, sprintf("Mean_Annual_Wind_Speed_%s_%s.tif", period_name, p$range))
  }
}

# ---- 4. SEA ICE (siconca, %): per-month climatologies, Oct-Feb -----------------
# Two products per month:
#   1. the concentration itself, cropped (not masked) to the 60 deg S extent
#      so every grid cell that overlaps that boundary is kept
#   2. the mean concentration within a set of buffer distances around every
#      ice-free-land cell (radii below -- easy to add/remove)
if (variable == "sea_ice") {
  
  # Reference layers that don't change by period
  Ant_extent <- vect(here("Data/PolarRes26/add_data_limit_v7.2.shp"))  # EPSG:3031, 60 deg S boundary
  coast <- vect(here("Data/PolarRes26/add_coastline_high_res_polygon_v7_12.shp"))
  SG <- vect(here("Data", "PolarRes26", "orkney.shp"))
  
  # Ice-free-land layers DO change by period: historical & mid use the
  # current domain, future uses the projected future domain
  domain_hist_mid <- rast(here("Data/PolarRes26/ice_free_domain.tif"))
  domain_future    <- rast(here("Data/PolarRes26/ice_free_future_domain.tif"))
  
  # Buffer radii (km) around ice-free-land cells -- malleable, just add/remove values
  sea_ice_buffer_km <- c(2, 5, 10, 50, 100)
  
  # The sea-ice season of interest: Oct, Nov, Dec, Jan, Feb. Oct/Nov/Dec are
  # drawn from the year BEFORE the target year (same austral-summer logic as
  # DJF elsewhere), Jan/Feb from the target year itself.
  sea_ice_months  <- c(10, 11, 12, 1, 2)
  sea_ice_offsets <- c(-1, -1, -1, 0, 0)
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- sea_ice: ", period_name)
    
    domain <- if (period_name == "FUTURE") domain_future else domain_hist_mid
    
    domain.pts <- as.points(domain, values = TRUE)
    domain.pts <- domain.pts[domain.pts$rock_union1 == 1, ]
    
    buffers <- lapply(sea_ice_buffer_km, function(km) terra::buffer(domain.pts, km * 1000))
    names(buffers) <- paste0(sea_ice_buffer_km, "km")
    
    series <- load_variable_series(model_dir, "siconca", p$scenario)
    
    for (i in seq_along(sea_ice_months)) {
      m <- sea_ice_months[i]
      offset <- sea_ice_offsets[i]
      month_label <- month.name[m]
      message("   month: ", month_label)
      
      conc <- climatological_monthly_mean_sea_ice(series$r, series$dates, m, p$years, offset, coast, SG)
      
      # Product 1: cropped to the 60S extent, keeping every overlapping cell
      conc_cropped <- terra::crop(conc, Ant_extent, snap = "out")
      queue_raster(conc_cropped, sprintf("Mean_%s_Sea_Ice_Concentration_%s_%s.tif",
                                         month_label, period_name, p$range))
      
      # Product 2: mean concentration within each buffer distance
      for (label in names(buffers)) {
        buffered <- extract_to_buffer(conc, domain, domain.pts, buffers[[label]])
        queue_raster(buffered, sprintf("Mean_%s_Sea_Ice_Concentration_%s_%s_%s.tif",
                                       month_label, label, period_name, p$range))
      }
    }
  }
}

# ---- 5a. TOTAL ANNUAL PRECIPITATION (pr, kg m-2 s-1 -> mm) -----------------------
if (variable == "total_precip") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- total_precip: ", period_name)
    
    series <- load_variable_series(model_dir, "pr", p$scenario)
    
    # daily -> monthly sum -> monthly climatology, saved individually
    monthly <- climatological_monthly_sums_all(series$r, series$dates, p$years, transform = flux_to_mm)
    for (m in 1:12) {
      queue_raster(monthly[[m]], sprintf("Climatological_Monthly_Total_Precipitation_%s_%s_%s.tif",
                                         month.name[m], period_name, p$range))
    }
    total <- annual_total_from_monthly(monthly)
    queue_raster(total, sprintf("Total_Annual_Precipitation_%s_%s.tif", period_name, p$range))
  }
}

# ---- 5b. TOTAL SUMMER (DJF) PRECIPITATION (pr, mm) -------------------------------
# NOTE / FLAGGED AS INCONSISTENT: left as a single per-season total (not run
# through the daily -> monthly-sum -> monthly-climatology pipeline used for
# the other total variables above and below). This means it's currently the
# only "total" variable that doesn't save monthly climatology outputs. Left
# as-is for now -- revisit if/when this should match the others.
if (variable == "total_summer_precip") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- total_summer_precip: ", period_name)
    
    series <- load_variable_series(model_dir, "pr", p$scenario)
    
    total <- summer_total_climatology(series$r, series$dates, p$years, transform = flux_to_mm)
    queue_raster(total, sprintf("Mean_Total_Summer_Precipitation_%s_%s.tif", period_name, p$range))
  }
}

# ---- 5c. MEAN ANNUAL PRECIPITATION (pr, mm/day): annual --------------------------
if (variable == "mean_precip") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- mean_precip: ", period_name)
    
    series <- load_variable_series(model_dir, "pr", p$scenario)
    r_mm <- flux_to_mm(series$r)
    
    monthly <- climatological_monthly_means_all(r_mm, series$dates, p$years)
    for (m in 1:12) {
      queue_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Precipitation_%s_%s_%s.tif",
                                         month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    queue_raster(annual, sprintf("Mean_Annual_Precipitation_%s_%s.tif", period_name, p$range))
  }
}

# ---- 5d. MEAN SUMMER (DJF) PRECIPITATION (pr, mm/day) ----------------------------
if (variable == "mean_summer_precip") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- mean_summer_precip: ", period_name)
    
    series <- load_variable_series(model_dir, "pr", p$scenario)
    r_mm <- flux_to_mm(series$r)
    
    summer <- summer_mean_climatology(r_mm, series$dates, p$years)
    queue_raster(summer, sprintf("Mean_Summer_Precipitation_%s_%s.tif", period_name, p$range))
  }
}

# ---- 6. SOLAR RADIATION (net SW = rsds - rsus, W m-2, no conversion) ------------
if (variable == "solar_rad") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- solar_rad: ", period_name)
    
    rsds <- load_variable_series(model_dir, "rsds", p$scenario)
    rsus <- load_variable_series(model_dir, "rsus", p$scenario)
    
    if (!identical(rsds$dates, rsus$dates)) {
      stop("rsds and rsus dates don't line up for ", model, " / ", p$scenario,
           " -- check Script 1's date-gap report for these two variables.")
    }
    swnet_r <- rsds$r - rsus$r
    
    monthly <- climatological_monthly_means_all(swnet_r, rsds$dates, p$years)
    for (m in 1:12) {
      queue_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Solar_Radiation_%s_%s_%s.tif",
                                         month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    queue_raster(annual, sprintf("Mean_Annual_Solar_Radiation_%s_%s.tif", period_name, p$range))
  }
}

# ---- 7. MEAN & TOTAL SNOW MELT (snm, converted to mm like precipitation) -------
# Mass of ice/snow melted at the surface per second (magnitude of melt)
if (variable == "mean_melt") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- mean_melt: ", period_name)
    
    series <- load_variable_series(model_dir, "snm", p$scenario)
    r_mm <- flux_to_mm(series$r)  # now converted to mm, same as precipitation
    
    monthly <- climatological_monthly_means_all(r_mm, series$dates, p$years)
    for (m in 1:12) {
      queue_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Melt_%s_%s_%s.tif",
                                         month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    queue_raster(annual, sprintf("Mean_Annual_Melt_%s_%s.tif", period_name, p$range))
  }
}

# ---- 7b. TOTAL ANNUAL SNOW MELT (snm, converted to mm like precipitation) -------
if (variable == "total_melt") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- total_melt: ", period_name)
    
    series <- load_variable_series(model_dir, "snm", p$scenario)
    
    monthly <- climatological_monthly_sums_all(series$r, series$dates, p$years, transform = flux_to_mm)
    for (m in 1:12) {
      queue_raster(monthly[[m]], sprintf("Climatological_Monthly_Total_Melt_%s_%s_%s.tif",
                                         month.name[m], period_name, p$range))
    }
    total <- annual_total_from_monthly(monthly)
    queue_raster(total, sprintf("Mean_Total_Annual_Melt_%s_%s.tif", period_name, p$range))
  }
}

# ---- 8. MEAN SNOW COVER (snc, %, no conversion): annual, monthly climatology only --
if (variable == "mean_snow") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- mean_snow: ", period_name)
    
    series <- load_variable_series(model_dir, "snc", p$scenario)
    
    monthly <- climatological_monthly_means_all(series$r, series$dates, p$years)
    for (m in 1:12) {
      queue_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Snow_Cover_%s_%s_%s.tif",
                                         month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    queue_raster(annual, sprintf("Mean_Annual_Snow_Cover_%s_%s.tif", period_name, p$range))
    
  }
}

# ---- 9. MEAN RELATIVE HUMIDITY (hurs, %): monthly, annual, summer, winter ------

if (variable == "mean_hurs") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- mean_hurs: ", period_name)
    
    series <- load_variable_series(model_dir, "hurs", p$scenario)
    
    monthly <- climatological_monthly_means_all(series$r, series$dates, p$years)
    for (m in 1:12) {
      queue_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Relative_Humidity_%s_%s_%s.tif",
                                         month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    queue_raster(annual, sprintf("Mean_Annual_Relative_Humidity_%s_%s.tif", period_name, p$range))
    
    summer <- summer_mean_climatology(series$r, series$dates, p$years)
    queue_raster(summer, sprintf("Mean_Summer_Relative_Humidity_%s_%s.tif", period_name, p$range))
    
    winter <- winter_mean_from_monthly(monthly)
    queue_raster(winter, sprintf("Mean_Winter_Relative_Humidity_%s_%s.tif", period_name, p$range))
  }
}

# ---- 10. VAPOUR PRESSURE DEFICIT (Pa): monthly, annual --------------------------
# Computed from hurs (relative humidity, %) and tas (air temperature), daily,
# then run through the same monthly-climatology + annual-mean pattern as
# temp/wind/etc.
#
# This uses the same Sonntag 1990 / WMO(2008) Magnus formula as
# bigleaf::rH.to.VPD() / Esat.slope() -- but reimplemented directly as raster
# arithmetic below (compute_vpd_Pa()) rather than calling those functions on
# the rasters directly. bigleaf's Esat.slope() ends with
# `data.frame(Esat, Delta)`, which doesn't work on a SpatRaster, so calling
# rH.to.VPD() on raster inputs would error. 
if (variable == "vpd") {
  
  # rH.to.VPD(rH, Tair) returns VPD in kPa; hurs is %, so rH = hurs/100
  compute_vpd_Pa <- function(hurs_pct, tas_celsius) {
    a <- 611.2; b <- 17.62; c <- 243.12  # Sonntag 1990 Magnus coefficients
    rH <- hurs_pct / 100
    esat_Pa <- a * exp((b * tas_celsius) / (c + tas_celsius))
    vpd_Pa <- esat_Pa * (1 - rH)
    vpd_Pa
  }
  
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- vpd: ", period_name)
    
    hurs_series <- load_variable_series(model_dir, "hurs", p$scenario)
    tas_series  <- load_variable_series(model_dir, "tas",  p$scenario)
    
    if (!identical(hurs_series$dates, tas_series$dates)) {
      stop("hurs and tas dates don't line up for ", model, " / ", p$scenario,
           " -- check Script 1's date-gap report for these two variables.")
    }
    
    tas_celsius <- to_celsius(tas_series$r)
    vpd_daily_Pa <- compute_vpd_Pa(hurs_series$r, tas_celsius)
    
    monthly <- climatological_monthly_means_all(vpd_daily_Pa, hurs_series$dates, p$years)
    for (m in 1:12) {
      queue_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_VPD_%s_%s_%s.tif",
                                         month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    queue_raster(annual, sprintf("Mean_Annual_VPD_%s_%s.tif", period_name, p$range))
  }
}

# ==============================================================================
# FINALISE: write every queued output to disk, in its own native grid.
# RACMO / ERA5-driven models were already fixed to the correct CRS + extent
# when each raster was read in (see load_variable_series() in section 2), so
# there is nothing left to do here except write -- no resampling onto the
# HCLIM grid.
# ==============================================================================

message("Finalising ", length(pending_outputs), " output(s) for ", model, " / ", variable, "...")
for (item in pending_outputs) {
  save_raster(item$r, item$filename)
}

message("Done: ", model, " / ", variable)