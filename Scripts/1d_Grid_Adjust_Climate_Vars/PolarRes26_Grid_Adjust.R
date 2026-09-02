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
# SEA ICE is a special case (see section 5b below): instead of the generic
# COASTLINE/ICEFREE/ICEFREE_FUTURE masks, every monthly sea-ice climatology
# file gets three other products, produced on the same reprojected/resampled
# grid, using the SAME historical/future ice-free-land inputs as every other
# variable:
#   - <name>_CONCENTRATION.tif   -- coast/SG masked, cropped to 60 deg S
#   - <name>_BUFFER_<radius>.tif -- mean concentration within each buffer
#                                   distance of every ice-free-land cell
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

input_base  <- here("Data/Environmental_predictors/PolarRes26")
output_base <- here("Data/Environmental_predictors/PolarRes26/Regridded")

coast_domain <- rast(here("Data/coast_domain.tif"))
ice_free_domain <- rast(here("Data/ice_free_domain.tif"))
ice_free_future_domain <- rast(here("Data/ice_free_future_domain.tif"))

# Template file used to recover a CRS for HCLIM inputs that come out of
# Script 1/2 with an empty/missing CRS (a known HCLIM quirk). Only ever
# used as a source of CRS metadata -- never as data -- and only after the
# extent of the broken input has been checked against this template's
# extent (see section 5).
hclim_crs_template_path <- here(
  "Data/PolarRes26/HCLIM_CESM2/historical/r11i1p1f1/HCLIM43-ALADIN/v1-r1/day/hurs/v20251130/",
  "hurs_ANT-12_CESM2_historical_r11i1p1f1_HCLIMcom-DMI_HCLIM43-ALADIN_v1-r1_day_19860101-19901231.nc"
)

# If TRUE, skip an input file entirely (no reprojection/resampling work at
# all) when every output it would produce already exists on disk. Set to
# FALSE to always reprocess and overwrite.
SKIP_EXISTING <- FALSE

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
job_index <- suppressWarnings(as.integer(args[1]))

# No argument, a non-numeric argument, or an out-of-range job_index: print
# the full job table (as documented above) instead of erroring out on a
# malformed subscript.
if (length(args) == 0 || is.na(job_index) || job_index < 1 || job_index > nrow(job_table)) {
  print(job_table, n = Inf)
  message(nrow(job_table), " job(s) total. Pass a job_index between 1 and ",
          nrow(job_table), " to run a single job.")
  quit(save = "no", status = 0)
}

this_job    <- job_table[job_index, ]
model       <- this_job$model
input_path  <- this_job$input_path
filename    <- this_job$filename
message("Job ", job_index, "/", nrow(job_table), " -> ", model, " / ", filename)

# ---- 4. Output paths -------------------------------------------------------------

out_dir <- file.path(output_base, model)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stem <- tools::file_path_sans_ext(filename)
coastline_path      <- file.path(out_dir, paste0(stem, "_COASTLINE.tif"))
icefree_path        <- file.path(out_dir, paste0(stem, "_ICEFREE.tif"))
icefree_futurepath  <- file.path(out_dir, paste0(stem, "_ICEFREE_FUTURE.tif"))

# Sea-ice detection only needs the filename, so it's done here (before any
# raster is read) -- both to route to the right branch below and to know
# which output paths to check for the skip-existing logic. Every monthly
# sea-ice climatology file gets the concentration/buffer treatment, not
# just a subset of months.
is_sea_ice <- grepl("Sea_Ice_Concentration", filename, fixed = TRUE)

# Period is also read off the filename. FUTURE-period files additionally
# get an ICEFREE_FUTURE output; historical/mid files only get the two
# generic COASTLINE/ICEFREE outputs. This applies uniformly to every
# variable, sea ice included.
is_future <- grepl("_FUTURE_", filename, fixed = TRUE)

sea_ice_buffer_km <- c(2, 5, 10, 50, 100)  # must match the buffer radii used in section 5a

# Every output path this job would produce -- used only for the
# skip-existing check. COASTLINE/ICEFREE are always produced; ICEFREE_FUTURE
# only for FUTURE-period files; CONCENTRATION/BUFFER only for sea ice (in
# addition to, not instead of, the COASTLINE/ICEFREE/ICEFREE_FUTURE set).
domain_mask_outputs <- c(coastline_path, icefree_path)
if (is_future) domain_mask_outputs <- c(domain_mask_outputs, icefree_futurepath)

expected_outputs <- if (is_sea_ice) {
  c(file.path(out_dir, paste0(stem, "_CONCENTRATION.tif")),
    file.path(out_dir, paste0(stem, "_BUFFER_", sea_ice_buffer_km, "km.tif")))
} else {
  domain_mask_outputs
}

# ---- 4b. Skip-existing check ------------------------------------------------------
# If every output this job would produce already exists, skip the whole job
# (no reprojection/resampling work at all).
if (SKIP_EXISTING && all(file.exists(expected_outputs))) {
  message("  all expected output(s) already exist -- skipping: ",
          paste(basename(expected_outputs), collapse = ", "))
  quit(save = "no", status = 0)
}

# ---- 5. Project, resample, mask, save --------------------------------------------

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
                    resolution = res(coast_domain)) # 10 km domain

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

# ---- 5a-helper. Shared COASTLINE / ICEFREE / ICEFREE_FUTURE writer -------------
# Used by every variable, including sea ice, on the already-reprojected/
# resampled raster `r`. ICEFREE_FUTURE is only written for FUTURE-period
# files and is saved as its own, separate output alongside COASTLINE/ICEFREE.
write_domain_masks <- function(r, is_future) {
  
  r_coast <- mask(r, coast_domain, maskvalue = NA)
  writeRaster(r_coast, coastline_path, overwrite = TRUE)
  message("  wrote ", coastline_path)
  
  r_icefree <- mask(r, ice_free_domain, maskvalue = NA)
  writeRaster(r_icefree, icefree_path, overwrite = TRUE)
  message("  wrote ", icefree_path)
  
  if (is_future) {
    r_icefree_future <- mask(r, ice_free_future_domain, maskvalue = NA)
    writeRaster(r_icefree_future, icefree_futurepath, overwrite = TRUE)
    message("  wrote ", icefree_futurepath)
  }
}

# ---- 5b. Sea-ice concentration files: special-case handling --------------------
# These skip the generic COASTLINE/ICEFREE/ICEFREE_FUTURE masking below and
# instead get the three products that used to be produced directly inside
# Script 2's old sea-ice block, now run on the already-reprojected/resampled
# (common-grid) raster `r` instead of on the raw daily model data:
#   1. coast/land masking + South Georgia / South America artefact removal
#   2. crop (not mask) to the 60 deg S extent -> concentration product
#   3. mean concentration within buffer distances of every ice-free-land
#      cell -> buffer products

if (is_sea_ice) {
  
  message("  sea-ice file detected -- running coast/SG masking, 60S crop, and buffer products")
  
  # Reference layers
  Ant_extent <- vect(here("Data/PolarRes26/add_data_limit_v7.2.shp"))  # EPSG:3031, 60 deg S boundary
  coast <- vect(here("Data/PolarRes26/add_coastline_high_res_polygon_v7_12.shp"))
  SG <- vect(here("Data", "PolarRes26", "orkney.shp"))
  
  # Ice-free-land layers DO change by period: historical & mid use the
  # current domain, future uses the projected future domain. Which one
  # applies is determined from the filename (it carries the period label);
  # `is_future` was already computed in section 4. Uses the SAME
  # ice_free_domain / ice_free_future_domain objects loaded once in section 1
  # for every other variable, rather than re-reading separate copies from a
  # different path -- so sea ice and every other variable are guaranteed to
  # be using identical ice-free-land inputs.
  sea_ice_domain <- if (is_future) ice_free_future_domain else ice_free_domain
  
  domain.pts <- as.points(sea_ice_domain, values = TRUE)
  domain.pts <- domain.pts[domain.pts$rock_union1 == 1, ]
  
  buffers <- lapply(sea_ice_buffer_km, function(km) terra::buffer(domain.pts, km * 1000))
  names(buffers) <- paste0(sea_ice_buffer_km, "km")
  
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
  
  # 1. Mask out land/coast (keep everything OUTSIDE the coastline polygon)
  conc <- mask(r, coast, inverse = TRUE, touches = FALSE)
  
  # Zero-out the South Georgia / South America cells that sit at 100% SIC
  # year-round -- this is a model-domain artefact, not real sea ice.
  IDs <- terra::extract(conc, SG, cells = TRUE)
  IDs <- IDs[!is.na(IDs[, 2]) & IDs[, 2] == 100, ]
  if (nrow(IDs) > 0) conc[IDs$cell] <- 0
  
  # 2. Cropped to the 60S extent, keeping every overlapping cell
  conc_cropped <- terra::crop(conc, Ant_extent, snap = "out")
  concentration_path <- file.path(out_dir, paste0(stem, "_CONCENTRATION.tif"))
  writeRaster(conc_cropped, concentration_path, overwrite = TRUE)
  message("  wrote ", concentration_path)
  
  # 3. Mean concentration within each buffer distance
  for (label in names(buffers)) {
    buffered <- extract_to_buffer(conc, sea_ice_domain, domain.pts, buffers[[label]])
    buffer_path <- file.path(out_dir, paste0(stem, "_BUFFER_", label, ".tif"))
    writeRaster(buffered, buffer_path, overwrite = TRUE)
    message("  wrote ", buffer_path)
  }
  
  message("Done: ", model, " / ", filename)
  quit(save = "no", status = 0)
}

# ---- 5c. Every other variable: generic COASTLINE / ICEFREE / ICEFREE_FUTURE ----

write_domain_masks(r, is_future)

message("Done: ", model, " / ", filename)