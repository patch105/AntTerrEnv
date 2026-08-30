# ==============================================================================
# PolarRes26 -- STEP 3: REPROJECT, RESAMPLE, AND MASK EVERY CLIMATOLOGY OUTPUT
# ==============================================================================
# Takes every .tif Script 2 wrote (for every model) and produces two masked
# versions of each, on a common grid:
#   - <name>_COASTLINE.tif -- masked to the coastline domain
#   - <name>_ICEFREE.tif   -- masked to the ice-free domain
#
# SEA ICE is a special case (see section 5a below): instead of the generic
# COASTLINE/ICEFREE masks, EVERY monthly sea-ice climatology file gets three
# other products produced on the same reprojected/resampled grid:
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

# models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5", "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5", "MetUM_ERA5")
models <- c("RACMO_MPI_ESM1", "RACMO_CESM2")

input_base  <- here("Data/Environmental_predictors/PolarRes26")
output_base <- here("Data/Environmental_predictors/PolarRes26/Regridded")

coast_domain    <- rast(here("Data/PolarRes26/coast_domain.tif"))
ice_free_domain <- rast(here("Data/PolarRes26/ice_free_domain.tif"))

# If TRUE, skip an input file entirely (no reprojection/resampling work at
# all) when every output it would produce already exists on disk. Set to
# FALSE to always reprocess and overwrite.
SKIP_EXISTING <- TRUE

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
coastline_path <- file.path(out_dir, paste0(stem, "_COASTLINE.tif"))
icefree_path   <- file.path(out_dir, paste0(stem, "_ICEFREE.tif"))

# Sea-ice detection only needs the filename, so it's done here (before any
# raster is read) -- both to route to the right branch below and to know
# which output paths to check for the skip-existing logic. Every monthly
# sea-ice climatology file gets the concentration/buffer treatment, not
# just a subset of months.
is_sea_ice <- grepl("Sea_Ice_Concentration", filename, fixed = TRUE)

sea_ice_buffer_km <- c(2, 5, 10, 50, 100)  # must match the buffer radii used in section 5a

# Every output path this job would produce, depending on which branch it
# takes below -- used only for the skip-existing check.
expected_outputs <- if (is_sea_ice) {
  c(file.path(out_dir, paste0(stem, "_CONCENTRATION.tif")),
    file.path(out_dir, paste0(stem, "_BUFFER_", sea_ice_buffer_km, "km.tif")))
} else {
  c(coastline_path, icefree_path)
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

# For now: Hardcode the resolution to 11km
reso <- 11000

if (is.na(crs(r)) || crs(r) == "") {
  stop("Input raster has no CRS: ", input_path,
       " -- reprojecting without a source CRS would silently produce garbage.",
       " Check this file (and Script 1/2's terra version, if this is a RACMO file).")
}

# Step 1: Make a domain for target CRS with matching resolution
target_grid <- rast(extent = ext(coast_domain), crs = crs(coast_domain),
                    resolution = reso)

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

# Step 3: Resample to 1km domain grid (already in the right CRS & extent, so
# "near" here is just chopping up 10km to 1km grids
r <- resample(r, coast_domain, method = "near")

# ---- 5a. Sea-ice concentration files: special-case handling --------------------
# These skip the generic COASTLINE/ICEFREE masking below and instead get the
# three products that used to be produced directly inside Script 2's old
# sea-ice block, now run on the already-reprojected/resampled (common-grid)
# raster `r` instead of on the raw daily model data:
#   1. coast/land masking + South Georgia / South America artefact removal
#   2. crop (not mask) to the 60 deg S extent -> concentration product
#   3. mean concentration within buffer distances of every ice-free-land
#      cell -> buffer products
# Applies to every monthly sea-ice climatology file (is_sea_ice and
# sea_ice_buffer_km were already determined in section 4, ahead of the
# skip-existing check).

if (is_sea_ice) {
  
  message("  sea-ice file detected -- running coast/SG masking, 60S crop, and buffer products")
  
  # Reference layers
  Ant_extent <- vect(here("Data/PolarRes26/add_data_limit_v7.2.shp"))  # EPSG:3031, 60 deg S boundary
  coast <- vect(here("Data/PolarRes26/add_coastline_high_res_polygon_v7_12.shp"))
  SG <- vect(here("Data", "PolarRes26", "orkney.shp"))
  
  # Ice-free-land layers DO change by period: historical & mid use the
  # current domain, future uses the projected future domain. Which one
  # applies is determined from the filename (it carries the period label).
  domain_hist_mid <- rast(here("Data/PolarRes26/ice_free_domain.tif"))
  domain_future    <- rast(here("Data/PolarRes26/ice_free_future_domain.tif"))
  domain <- if (grepl("_FUTURE_", filename, fixed = TRUE)) domain_future else domain_hist_mid
  
  domain.pts <- as.points(domain, values = TRUE)
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
    buffered <- extract_to_buffer(conc, domain, domain.pts, buffers[[label]])
    buffer_path <- file.path(out_dir, paste0(stem, "_BUFFER_", label, ".tif"))
    writeRaster(buffered, buffer_path, overwrite = TRUE)
    message("  wrote ", buffer_path)
  }
  
  message("Done: ", model, " / ", filename)
  quit(save = "no", status = 0)
}

# Mask to the coastline domain
r_coast <- mask(r, coast_domain, maskvalue = NA)
writeRaster(r_coast, coastline_path, overwrite = TRUE)
message("  wrote ", coastline_path)

# Mask to the ice-free domain
r_icefree <- mask(r, ice_free_domain, maskvalue = NA)
writeRaster(r_icefree, icefree_path, overwrite = TRUE)
message("  wrote ", icefree_path)

message("Done: ", model, " / ", filename)