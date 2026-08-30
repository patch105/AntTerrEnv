# ==============================================================================
# PolarRes26 -- STEP 2b: SUMMARISE INTO BIOCLIM VARIABLES
# ==============================================================================
# Reads the monthly climatology .tif files produced by PolarRes26_Summarising.R
# (12 months each of TasMin, TasMax, and Total Precipitation, per model and
# period) and runs them through predicts::bcvars() to get the standard 19
# WorldClim-style bioclimatic variables, saved as one .tif per BIO layer
# (BIO1, BIO2, ... BIO19) per model/period.
#
# One job = one model (all periods for that model are done in the same job).
# Run with no argument (or an out-of-range one) to print the model lookup
# table instead of failing silently.
#
# DEPENDENCY: this needs PolarRes26_Summarising.R's "temp_minmax" variable to
# have been run (it produces the TasMin/TasMax monthly climatologies) and
# "total_precip" to have been run (monthly precipitation totals). If either
# is missing for a given model/period, this stops with a message telling you
# exactly which file(s) are missing.
#
# bcvars() expects prec/tmin/tmax as three 12-layer SpatRasters, with layers
# in calendar order (Jan first). Filenames are built explicitly in
# month.name order (Jan-Dec) below to guarantee that, rather than relying on
# whatever order list.files() would return.
#
# tmin/tmax are expected in degC (which is what "temp_minmax" produces) and
# prec in mm (which is what "total_precip" produces).
# ==============================================================================

# ---- 0. Setup ------------------------------------------------------------------

lib_loc <- paste(getwd(), "/r_lib_new", sep = "")

library(terra)
library(here)
library(predicts, lib.loc = lib_loc)

# ---- 1. Configuration -----------------------------------------------------------

# models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5", "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5", "MetUM_ERA5")
models <- c("RACMO_CESM2", "RACMO_MPI_ESM1")

years_hist   <- seq(1995, 2014, by = 1)
years_mid    <- seq(2041, 2060, by = 1)
years_future <- seq(2081, 2100, by = 1)

input_base  <- here("Data/Environmental_predictors/PolarRes26")
output_base <- here("Data/Environmental_predictors/PolarRes26_Bioclim")

# One job = one model. job_index (1-7) picks a row from `models`.
args <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])


model <- models[job_index]
message("Job ", job_index, "/", length(models), " -> Model: ", model)

# Just the period labels + filename range strings needed to find the
# already-computed monthly climatology outputs -- doesn't need
# years/scenario/RACMO-clipping info the way the main summarising script
# does, since no daily data is read here.
#
# ERA5-driven models (HCLIM_ERA5, RACMO_ERA5, MetUM_ERA5) only ever have a
# HISTORICAL set of monthly climatology files -- PolarRes26_Summarising.R
# doesn't produce MID/FUTURE outputs for them (no ssp series exists), so
# looking for those here would just fail with a missing-file error.
is_era5 <- endsWith(model, "_ERA5")

periods <- if (is_era5) {
  list(HISTORICAL = paste(min(years_hist), max(years_hist), sep = "_"))
} else {
  list(
    HISTORICAL = paste(min(years_hist),   max(years_hist),   sep = "_"),
    MID        = paste(min(years_mid),    max(years_mid),    sep = "_"),
    FUTURE     = paste(min(years_future), max(years_future), sep = "_")
  )
}

model_in_dir <- file.path(input_base, model)
if (!dir.exists(model_in_dir)) {
  stop("Model directory not found: ", model_in_dir)
}

model_out_dir <- file.path(output_base, model)
dir.create(model_out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Helper: load one variable's 12 monthly climatology files, in order -----

load_monthly_climatology <- function(model_dir, name_pattern, period_name, range_label) {
  filenames <- sprintf(name_pattern, month.name, period_name, range_label)
  paths <- file.path(model_dir, filenames)
  
  missing <- !file.exists(paths)
  if (any(missing)) {
    stop("Missing monthly climatology file(s) -- has the matching variable ",
         "been run for this model/period yet?\n",
         paste(paths[missing], collapse = "\n"))
  }
  
  # Loaded in month.name (Jan-Dec) order, matching the order `filenames`
  # was built in -- this is what bcvars() requires.
  rast(paths)
}

# ---- 3. Build bioclim variables for this model, every period -------------------

for (period_name in names(periods)) {
  range_label <- periods[[period_name]]
  message("-- bioclim: ", model, " / ", period_name)
  
  tmin <- load_monthly_climatology(model_in_dir, "Climatological_Monthly_Mean_TasMin_%s_%s_%s.tif",
                                   period_name, range_label)
  tmax <- load_monthly_climatology(model_in_dir, "Climatological_Monthly_Mean_TasMax_%s_%s_%s.tif",
                                   period_name, range_label)
  prec <- load_monthly_climatology(model_in_dir, "Climatological_Monthly_Total_Precipitation_%s_%s_%s.tif",
                                   period_name, range_label)
  
  bio <- predicts::bcvars(prec, tmin, tmax, filename = "")
  
  # Save each BIO layer separately (BIO1, BIO2, ... BIO19) rather than one
  # multi-layer file
  for (i in seq_len(nlyr(bio))) {
    out_path <- file.path(model_out_dir, sprintf("BIO%d_%s_%s.tif", i, period_name, range_label))
    writeRaster(bio[[i]], out_path, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  }
  message("  wrote BIO1-BIO", nlyr(bio), " for ", period_name, " to ", model_out_dir)
}

message("Done: ", model)