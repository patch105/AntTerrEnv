# ==============================================================================
# PolarRes26 -- STEP 3: REPROJECT, RESAMPLE, AND MASK EVERY CLIMATOLOGY OUTPUT
# ==============================================================================
# Takes every .tif Script 2 wrote (for every model) and produces masked
# versions of each, on a common grid:
#   - <name>_COASTLINE.tif       -- masked to the coastline domain
#   - <name>_ICEFREE.tif         -- masked to the (historical/mid) ice-free domain
#   - <name>_ICEFREE_FUTURE.tif  -- masked to the future ice-free domain
#                                    (FUTURE-period files only, saved as a
#                                    separate, additional output)
#
# Deliberately does NOT hard-code the list of variables/periods/products --
# it just finds every .tif file Script 2 actually produced, for every model,
# and processes each one. That way this script automatically covers
# whatever Script 2 currently outputs (including anything added, renamed,
# or removed there later) without the two scripts needing to be kept in
# sync by hand.
#
# One job = one input file, selected via a single job_index (same pattern as
# Script 2) -- run with no argument (or an out-of-range one) to print the
# full table of every file that will be processed, and how many jobs that is.
# ==============================================================================

# ---- 0. Setup ------------------------------------------------------------------


lib_loc <- paste(getwd(), "/r_lib_new", sep = "")

library(dplyr, lib.loc = lib_loc)
library(purrr, lib.loc = lib_loc)
library(terra)
library(here)

# ---- 1. Configuration -----------------------------------------------------------

models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5", "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5", "MetUM_ERA5")

input_base  <- here("Data/Environmental_predictors/PolarRes26_Bioclim")
output_base <- here("Data/Environmental_predictors/PolarRes26_Bioclim/Regridded")

coast_domain           <- rast(here("Data/PolarRes26/coast_domain.tif"))
ice_free_domain        <- rast(here("Data/PolarRes26/ice_free_domain.tif"))
ice_free_future_domain <- rast(here("Data/PolarRes26/ice_free_future_domain.tif"))

# Template file used to recover a CRS for HCLIM inputs that come out of
# Script 1/2 with an empty/missing CRS (a known HCLIM quirk). Only ever
# used as a source of CRS metadata -- never as data -- and only after the
# extent of the broken input has been checked against this template's
# extent (see section 5).
hclim_crs_template_path <- here(
  "Data/PolarRes26/HCLIM_CESM2/historical/r11i1p1f1/HCLIM43-ALADIN/v1-r1/day/hurs/v20251130/",
  "hurs_ANT-12_CESM2_historical_r11i1p1f1_HCLIMcom-DMI_HCLIM43-ALADIN_v1-r1_day_19860101-19901231.nc"
)


# ---- 2. Build the job table: every .tif Script 2 produced, across every model --

job_table <- map_dfr(models, function(model) {
  model_dir <- file.path(input_base, model)
  if (!dir.exists(model_dir)) {
    message("Model directory not found, skipping: ", model_dir)
    return(tibble())
  }
  files <- list.files(model_dir, pattern = "\\.tif$", full.names = TRUE)
  tibble(model = model, input_path = files, filename = basename(files))
})

if (nrow(job_table) == 0) {
  stop("No .tif files found under ", input_base, " -- has Script 2 been run yet?")
}

message(nrow(job_table), " output file(s) found across ",
        length(unique(job_table$model)), " model(s).")

# ---- 3. job_index selects ONE input file (same pattern as Script 2) ------------

args <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])


this_job    <- job_table[job_index, ]
model       <- this_job$model
input_path  <- this_job$input_path
filename    <- this_job$filename
message("Job ", job_index, "/", nrow(job_table), " -> ", model, " / ", filename)

# ---- 4. Output paths -------------------------------------------------------------

out_dir <- file.path(output_base, model)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stem <- tools::file_path_sans_ext(filename)
coastline_path     <- file.path(out_dir, paste0(stem, "_COASTLINE.tif"))
icefree_path       <- file.path(out_dir, paste0(stem, "_ICEFREE.tif"))
icefree_futurepath <- file.path(out_dir, paste0(stem, "_ICEFREE_FUTURE.tif"))

# Period is read off the filename. FUTURE-period files additionally get an
# ICEFREE_FUTURE output; historical/mid files only get the two generic
# COASTLINE/ICEFREE outputs.
is_future <- grepl("_FUTURE_", filename, fixed = TRUE)

# ---- 5. Project, mask, save --------------------------------------------------

r <- rast(input_path)

if (is.na(crs(r)) || crs(r) == "") {
  
  # Known HCLIM quirk: some HCLIM outputs come out of Script 1/2 with no
  # CRS attached, even though their grid/extent is fine. For HCLIM only,
  # recover the CRS from a known-good HCLIM template file -- but only after
  # confirming the extents actually line up, so we never silently stamp a
  # wrong CRS onto a raster that doesn't actually match the template grid.
  if (grepl("^HCLIM", model)) {
    
    hclim_template <- rast(hclim_crs_template_path)
    
    extents_match <- isTRUE(all.equal(
      as.vector(ext(r)), as.vector(ext(hclim_template)),
      tolerance = 1e-6
    ))
    
    if (!extents_match) {
      stop("Input raster has no CRS and its extent does not match the HCLIM ",
           "CRS template (", hclim_crs_template_path, ") -- refusing to guess ",
           "a CRS for a grid that doesn't line up. Input: ", input_path,
           "\n  input extent:    ", paste(round(as.vector(ext(r)), 4), collapse = ", "),
           "\n  template extent: ", paste(round(as.vector(ext(hclim_template)), 4), collapse = ", "))
    }
    
    crs(r) <- crs(hclim_template)
    message("  input had no CRS -- extent matched the HCLIM template, so CRS was ",
            "copied from it: ", input_path)
    
  } else {
    # RACMO / MetUM (or anything else): an empty CRS here is not a known,
    # safe-to-patch quirk -- treat it as the hard error it always has been.
    stop("Input raster has no CRS: ", input_path,
         " -- reprojecting without a source CRS would silently produce garbage.",
         " Check this file (and Script 1/2's terra version, if this is a RACMO file).")
  }
}

# Step 1: Make a domain for target CRS with matching resolution
target_grid <- rast(extent = ext(coast_domain), crs = crs(coast_domain),
                    resolution = res(coast_domain))

# Sanity check: extent must be a clean integer multiple of the resolution,
# otherwise project()/rast() will silently produce a grid that doesn't
# tile evenly and downstream resample() will misalign.
stopifnot(
  (xmax(target_grid) - xmin(target_grid)) %% res(target_grid) == 0,
  (ymax(target_grid) - ymin(target_grid)) %% res(target_grid) == 0
)


#extent values are clean multiples of 0 and the resolution

# Step 2: Reproject into the domain's CRS with bilinear interpolation
r <- project(r, target_grid, method = "bilinear")

# Mask to the coastline domain
r_coast <- mask(r, coast_domain, maskvalue = NA)
writeRaster(r_coast, coastline_path, overwrite = TRUE)
message("  wrote ", coastline_path)

# Mask to the ice-free domain
r_icefree <- mask(r, ice_free_domain, maskvalue = NA)
writeRaster(r_icefree, icefree_path, overwrite = TRUE)
message("  wrote ", icefree_path)

# Mask to the future ice-free domain (FUTURE-period files only)
if (is_future) {
  r_icefree_future <- mask(r, ice_free_future_domain, maskvalue = NA)
  writeRaster(r_icefree_future, icefree_futurepath, overwrite = TRUE)
  message("  wrote ", icefree_futurepath)
}

message("Done: ", model, " / ", filename)