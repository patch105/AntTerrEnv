
# ==============================================================================
# PolarRes26 -- STEP 2b: SUMMARISE THE "ADDITIONAL" HCLIM_CESM2 VARIABLES
#   (siconca, rsus, snm, snc) INTO CLIMATOLOGIES
#
# WHY THIS SCRIPT EXISTS
# -----------------------------------------------------------------------------
# Script 2's load_variable_series() assumes every variable lives as a set of
# files that, once read with rast(files, subds = variable_name), form
# one clean continuous daily (or sub-daily) time series with no spillover
# layers. That's true for most PolarRes26 variables, but siconca / rsus /
# snm / snc for HCLIM_CESM2 are still delivered in Script 1's ORIGINAL
# format: one file per calendar month, and (as documented in Script 1) each
# monthly file's LAST layer is actually the first timestep of the *next*
# month, which has to be dropped before the file's layers can be trusted.
#
# Rather than reusing Script 1's hard-coded hourly/6-hourly/daily layer-count
# math to find and drop that spillover layer, load_variable_series_monthly_
# chunks() below does it generically: it reads each file's REAL per-layer
# dates via terra::time() (the same mechanism Script 2 already trusts) and
# keeps only the layers whose date falls in that file's own calendar month.
# This also means a corrupt/missing monthly file just results in that
# year-month combination being absent from the assembled series -- which
# climatological_monthly_mean()/_sum() (copied from Script 2, unmodified)
# already knows to skip with a warning. So the "insert an NA filler and
# explicitly skip that month" step Script 1 did by hand for known corrupt
# files (siconca: Apr/Jun 1994, snc: Nov 1995, rsus: May 1997) happens here
# automatically, for ANY missing month -- report_missing_year_months() below
# prints exactly what got skipped so you can cross-check it against Script
# 1's list of known corrupt files (and catch anything new).
#
# Everything else -- the climatology maths, period definitions, output
# naming, output folder -- is taken as-is from Script 2, so outputs from
# this script sit alongside Script 2's HCLIM_CESM2 outputs seamlessly.
#
# DIAGNOSTICS ASSUMED (same caveat as Script 2): each file's time dimension
# parses to real, non-NA dates via terra::time(). If that turns out not to
# be true for one of these variables, this is the function to adjust.
# ==============================================================================

# ---- 0. Setup ------------------------------------------------------------------

lib_loc <- paste(getwd(), "/r_lib_new", sep = "")
# .libPaths(lib_loc)

library(terra)
library(here)
library(lubridate)
library(stringr)

# ---- 1. Configuration -----------------------------------------------------------

model <- "HCLIM_CESM2"  # this script is only ever run for this one model

# Standard PolarRes26 folder (continuous-file variables, e.g. rsds) --
# same base_dir Script 2 uses.
base_dir  <- "/mnt/hpccs01/home/patterc2/n11222026/AntarcticFutureHabitat/Data/PolarRes26"
model_dir <- file.path(base_dir, model)

# The "additional" folder holding siconca / rsus / snm / snc in Script 1's
# original monthly-chunked-file format.
additional_dir <- here("Data/PolarRes26/HCLIM_CESM2_additional")

ssp_scenario <- "ssp370"

years_hist   <- seq(1995, 2014, by = 1)
years_mid    <- seq(2041, 2060, by = 1)
years_future <- seq(2081, 2100, by = 1)

# One job = one variable.
args <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])

variables <- list("sea_ice", "solar_rad_net", "mean_melt", "total_melt", "mean_snow")
variable  <- variables[[job_index]]
message("Job ", job_index, " -> Model: ", model, " | Variable: ", variable)

# Output goes to the SAME place as Script 2's HCLIM_CESM2 outputs -- no
# "_additional" in the output path, per the naming convention requested.
outpath <- here("Data/Environmental_predictors/PolarRes26", model)
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

# HCLIM_CESM2 needs no CRS/extent gridfix (that's only for RACMO/ERA5-driven
# grids in Script 2), and it isn't ERA5-driven, so it gets the full
# HISTORICAL/MID/FUTURE period set with no RACMO-style future-year clipping.
periods <- list(
  HISTORICAL = list(years = years_hist,   scenario = "historical",
                    range = paste(min(years_hist),   max(years_hist),   sep = "_")),
  MID        = list(years = years_mid,    scenario = ssp_scenario,
                    range = paste(min(years_mid),    max(years_mid),    sep = "_")),
  FUTURE     = list(years = years_future, scenario = ssp_scenario,
                    range = paste(min(years_future), max(years_future), sep = "_"))
)

# ---- 2a. Loader for STANDARD continuous files (identical to Script 2) ----------
# Needed here only for rsds, which is NOT one of the "additional" variables --
# it already lives in the normal continuous-file format under model_dir, so
# solar_rad_net can load it exactly the way Script 2 does.

find_variable_files <- function(model_dir, variable_name, scenario) {
  all_files <- list.files(model_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
  if (length(all_files) == 0) return(character(0))
  components_list <- str_split(all_files, "/")
  is_var  <- vapply(components_list, function(cmp) variable_name %in% cmp, logical(1))
  is_scen <- grepl(scenario, all_files, fixed = TRUE)
  sort(all_files[is_var & is_scen])
}

load_variable_series <- function(model_dir, variable_name, scenario) {
  files <- find_variable_files(model_dir, variable_name, scenario)
  if (length(files) == 0) {
    stop(sprintf("No files found for variable '%s', scenario '%s' under %s",
                 variable_name, scenario, model_dir))
  }
  r <- rast(files, subds = variable_name)
  dates <- as.Date(terra::time(r))
  if (length(dates) != terra::nlyr(r) || anyNA(dates)) {
    stop(sprintf("Could not read a clean daily time dimension for '%s' (%s).",
                 variable_name, scenario))
  }
  ord <- order(dates)
  list(r = r[[ord]], dates = dates[ord])
}

# ---- 2b. Loader for the "additional" MONTHLY-CHUNK files (siconca/rsus/snm/snc) --

# Matches files by a leading "variable_" prefix on the filename (e.g.
# "rsus_fp_..."), same substring-on-scenario match Script 2 uses.
find_monthly_files <- function(dir, variable_name, scenario) {
  all_files <- list.files(dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
  if (length(all_files) == 0) return(character(0))
  is_var  <- grepl(paste0("^", variable_name, "_"), basename(all_files))
  is_scen <- grepl(scenario, all_files, fixed = TRUE)
  sort(all_files[is_var & is_scen])
}

# See the header comment block at the top of this script for the reasoning
# behind this function. In short: read each monthly file's real per-layer
# dates, keep only the layers that fall in that file's own calendar month
# (this is what drops the spillover layer Script 1 handled by hand), then
# concatenate everything into one continuous, date-sorted series.
load_variable_series_monthly_chunks <- function(dir, variable_name, scenario) {
  
  files <- find_monthly_files(dir, variable_name, scenario)
  if (length(files) == 0) {
    stop(sprintf("No monthly-chunk files found for variable '%s', scenario '%s' under %s",
                 variable_name, scenario, dir))
  }
  
  all_r     <- list()
  all_dates <- as.Date(character(0))
  
  for (f in files) {
    
    r <- rast(f)
    dates <- as.Date(terra::time(r))
    
    if (length(dates) != terra::nlyr(r) || anyNA(dates)) {
      stop(sprintf("Could not read a clean time dimension from file:\n  %s\n", f),
           "Re-check this file's CF time metadata before trusting this run.")
    }
    
    # This file's "own" calendar month = the month of its FIRST layer. Any
    # trailing layer(s) dated into the following month are the spillover
    # and get dropped here.
    own_year  <- lubridate::year(dates[1])
    own_month <- lubridate::month(dates[1])
    keep <- lubridate::year(dates) == own_year & lubridate::month(dates) == own_month
    
    if (!any(keep)) {
      warning(sprintf("File contributed no layers within its own month -- check:\n  %s", f))
      next
    }
    
    all_r[[length(all_r) + 1]] <- r[[which(keep)]]
    all_dates <- c(all_dates, dates[keep])
  }
  
  if (length(all_r) == 0) {
    stop(sprintf("No usable layers assembled for variable '%s', scenario '%s'.",
                 variable_name, scenario))
  }
  
  r_all <- rast(all_r)
  
  # De-duplicate on date (keep first occurrence) in case two files' own-month
  # windows ever overlap, then sort chronologically.
  dup <- duplicated(all_dates)
  if (any(dup)) {
    warning(sprintf("%d duplicate date(s) found for '%s' / '%s' -- keeping first occurrence.",
                    sum(dup), variable_name, scenario))
    r_all     <- r_all[[!dup]]
    all_dates <- all_dates[!dup]
  }
  
  ord <- order(all_dates)
  list(r = r_all[[ord]], dates = all_dates[ord])
}

# Prints which (year, month) combinations are absent from an assembled
# series, for the years actually needed by this period. This is the
# generic stand-in for Script 1's hand-maintained list of known corrupt
# files (siconca: Apr/Jun 1994; snc: Nov 1995; rsus: May 1997) -- anything
# it reports here is a month climatological_monthly_mean()/_sum() will skip
# automatically. Cross-check the printed list against Script 1's comments
# to confirm nothing new has gone missing.
report_missing_year_months <- function(dates, years, variable_name, scenario) {
  present  <- unique(format(dates, "%Y-%m"))
  expected <- as.vector(outer(years, sprintf("%02d", 1:12), paste, sep = "-"))
  missing  <- setdiff(expected, present)
  if (length(missing) > 0) {
    message(sprintf("  [%s / %s] missing %d year-month(s) (auto-skipped by the climatology functions): %s",
                    variable_name, scenario, length(missing), paste(sort(missing), collapse = ", ")))
  } else {
    message(sprintf("  [%s / %s] no missing year-months detected for %d-%d.",
                    variable_name, scenario, min(years), max(years)))
  }
}

# ---- 3. Climatology maths (identical to Script 2) --------------------------------

climatological_monthly_mean <- function(r, dates, month, years) {
  yearly <- list()
  for (y in years) {
    idx <- which(lubridate::month(dates) == month & lubridate::year(dates) == y)
    if (length(idx) == 0) {
      warning(sprintf("No data for month %d, year %d -- skipping", month, y))
      next
    }
    sub <- r[[idx]]
    yearly[[as.character(y)]] <- app(sub, mean, na.rm = TRUE)
  }
  if (length(yearly) == 0) stop("No years had any data -- check the date range and scenario coverage.")
  app(rast(yearly), mean, na.rm = TRUE)
}

climatological_monthly_means_all <- function(r, dates, years) {
  out <- vector("list", 12)
  names(out) <- month.name
  for (m in 1:12) {
    out[[m]] <- climatological_monthly_mean(r, dates, m, years)
  }
  out
}

annual_mean_from_monthly <- function(monthly_list) {
  app(rast(monthly_list), mean, na.rm = TRUE)
}

annual_total_from_monthly <- function(monthly_list) {
  app(rast(monthly_list), sum, na.rm = TRUE)
}

seasonal_mean_from_monthly <- function(monthly_list, months) {
  app(rast(monthly_list[months]), mean, na.rm = TRUE)
}

seasonal_total_from_monthly <- function(monthly_list, months) {
  app(rast(monthly_list[months]), sum, na.rm = TRUE)
}

climatological_monthly_sum <- function(r, dates, month, years, transform = NULL) {
  yearly <- list()
  for (y in years) {
    idx <- which(lubridate::month(dates) == month & lubridate::year(dates) == y)
    if (length(idx) == 0) {
      warning(sprintf("No data for month %d, year %d -- skipping", month, y))
      next
    }
    sub <- r[[idx]]
    if (!is.null(transform)) sub <- transform(sub)
    yearly[[as.character(y)]] <- app(sub, sum, na.rm = TRUE)
  }
  if (length(yearly) == 0) stop("No years had any data -- check the date range and scenario coverage.")
  app(rast(yearly), mean, na.rm = TRUE)
}

climatological_monthly_sums_all <- function(r, dates, years, transform = NULL) {
  out <- vector("list", 12)
  names(out) <- month.name
  for (m in 1:12) {
    out[[m]] <- climatological_monthly_sum(r, dates, m, years, transform = transform)
  }
  out
}

# ---- 4. Unit-conversion helpers (identical to Script 2) --------------------------

flux_to_mm <- function(r) r * 86400  # daily flux (kg m-2 s-1) to daily mm

# ---- 5. Output writing (identical to Script 2) ------------------------------------

save_raster <- function(r, filename) {
  path <- file.path(outpath, filename)
  writeRaster(r, path, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  message("  wrote ", filename)
}

# ==============================================================================
# 6. Create the variables
# ==============================================================================

# ---- SEA ICE (siconca, %): monthly climatology, full grid, no masking ----------
# Same simplified treatment as Script 2's sea_ice block (12 monthly
# climatologies only, no land mask / buffer / South Georgia correction),
# loaded via the monthly-chunk loader to handle the spillover layer.
if (variable == "sea_ice") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- sea_ice: ", period_name)
    
    series <- load_variable_series_monthly_chunks(additional_dir, "siconca", p$scenario)
    report_missing_year_months(series$dates, p$years, "siconca", p$scenario)
    
    monthly <- climatological_monthly_means_all(series$r, series$dates, p$years)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Sea_Ice_Concentration_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
  }
}

# ---- SOLAR RADIATION NET (net SW = rsds - rsus, W m-2) ---------------------------
# rsds loads the normal (continuous-file) way, exactly as in Script 2.
# rsus loads via the monthly-chunk loader, which will also transparently
# absorb the known corrupt May-1997 rsus file (see report_missing_year_
# months() output below) the same way Script 1's manual NA-fill-then-skip
# did -- that month is simply excluded from the days used to compute the
# climatology for years containing it.
if (variable == "solar_rad_net") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- solar_rad_net: ", period_name)
    
    rsds <- load_variable_series(model_dir, "rsds", p$scenario)
    rsus <- load_variable_series_monthly_chunks(additional_dir, "rsus", p$scenario)
    report_missing_year_months(rsus$dates, p$years, "rsus", p$scenario)
    
    common_dates <- as.Date(intersect(rsds$dates, rsus$dates), origin = "1970-01-01")
    if (length(common_dates) == 0) {
      stop("rsds and rsus share no common dates for ", model, " / ", p$scenario)
    }
    
    rsds_r <- rsds$r[[match(common_dates, rsds$dates)]]
    rsus_r <- rsus$r[[match(common_dates, rsus$dates)]]
    
    swnet_r <- rsds_r - rsus_r
    
    monthly <- climatological_monthly_means_all(swnet_r, common_dates, p$years)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Net_Solar_Radiation_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    save_raster(annual, sprintf("Mean_Annual_Net_Solar_Radiation_%s_%s.tif", period_name, p$range))
  }
}

# ---- MEAN & TOTAL SNOW MELT (snm, converted to mm like precipitation) ------------
# Same maths as Script 2's mean_melt/total_melt blocks, loaded via the
# monthly-chunk loader for the spillover-layer idiosyncrasy.
if (variable == "mean_melt") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- mean_melt: ", period_name)
    
    series <- load_variable_series_monthly_chunks(additional_dir, "snm", p$scenario)
    report_missing_year_months(series$dates, p$years, "snm", p$scenario)
    r_mm <- flux_to_mm(series$r)
    
    monthly <- climatological_monthly_means_all(r_mm, series$dates, p$years)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Melt_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    save_raster(annual, sprintf("Mean_Annual_Melt_%s_%s.tif", period_name, p$range))
  }
}

if (variable == "total_melt") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- total_melt: ", period_name)
    
    series <- load_variable_series_monthly_chunks(additional_dir, "snm", p$scenario)
    report_missing_year_months(series$dates, p$years, "snm", p$scenario)
    
    monthly <- climatological_monthly_sums_all(series$r, series$dates, p$years, transform = flux_to_mm)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Total_Melt_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    total <- annual_total_from_monthly(monthly)
    save_raster(total, sprintf("Mean_Total_Annual_Melt_%s_%s.tif", period_name, p$range))
  }
}

# ---- MEAN SNOW COVER (snc, %): annual + monthly climatology ----------------------
# Same maths as Script 2's mean_snow block. The known corrupt Nov-1995 snc
# file (if it's still missing from the additional folder) is auto-skipped
# by the loader/climatology functions -- check the report_missing_year_
# months() message below to confirm it's the only gap.
if (variable == "mean_snow") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- mean_snow: ", period_name)
    
    series <- load_variable_series_monthly_chunks(additional_dir, "snc", p$scenario)
    report_missing_year_months(series$dates, p$years, "snc", p$scenario)
    
    monthly <- climatological_monthly_means_all(series$r, series$dates, p$years)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Snow_Cover_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    save_raster(annual, sprintf("Mean_Annual_Snow_Cover_%s_%s.tif", period_name, p$range))
  }
}

message("Done: ", model, " / ", variable)