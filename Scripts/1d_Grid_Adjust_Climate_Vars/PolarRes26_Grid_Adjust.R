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
# SEA ICE is a special case (see sections 5b-5e below): instead of the
# generic COASTLINE/ICEFREE/ICEFREE_FUTURE masks, every monthly sea-ice
# climatology file gets three other products, produced on the same
# reprojected/resampled grid, using the SAME historical/future ice-free-land
# inputs as every other variable:
#   - <name>_CONCENTRATION.tif   -- coast/SG masked, cropped to 60 deg S
#   - <name>_BUFFER_<radius>.tif -- mean concentration within each buffer
#                                   distance of every ice-free-land cell
#
# For RACMO sea-ice files specifically, an extra landmask fix runs first
# (section 5b): RACMO sea-ice output sits on a rotated-pole grid that comes
# out of Script 1/2 without a usable CRS, and the RACMO landmask
# (ANT11_masks.nc) is offset from it in index space. A fixed template sic
# file (always the same one, regardless of which RACMO job is running) is
# used to recover the rotated-pole CRS and the index-space offset against
# the landmask; that CRS + landmask are then applied to THIS job's actual
# raster before it's reprojected. HCLIM and MetUM sea ice (and every other
# RACMO variable) are untouched by this and go through the normal
# CRS-handling in section 5a.
#
# Sea-ice files also get a non-interactive coastal gap-fill after
# reprojection (section 5d): small NA gaps left by bilinear reprojection
# near the coastline are patched via patch-labelling + focal mean, with a
# fixed size_thresh of 100 cells for every model -- no histogram
# inspection, no manual tuning. Patches touching the raster edge are never
# filled, regardless of size, since those are out-of-domain rather than
# local artifacts.
#
# Deliberately does NOT hard-code the list of variables/periods/products --
# it just finds every .tif file Script 2 actually produced, for every
# model, and processes each one. That way this script automatically covers
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
# extent (see section 5a).
hclim_crs_template_path <- here(
  "Data/PolarRes26/HCLIM_CESM2/historical/r11i1p1f1/HCLIM43-ALADIN/v1-r1/day/hurs/v20251130/",
  "hurs_ANT-12_CESM2_historical_r11i1p1f1_HCLIMcom-DMI_HCLIM43-ALADIN_v1-r1_day_19860101-19901231.nc"
)

# Fixed RACMO sic template + landmask, used ONLY to derive the rotated-pole
# CRS ("crsfix") and the index-space offset ("off") between the RACMO
# native grid and the RACMO landmask. Always the same template file,
# regardless of which RACMO sea-ice job is actually being processed
# (see section 5b).
racmo_sic_template_path <- here(
  "Data/PolarRes26/RACMO_CESM2/historical/r11i1p1f1/RACMO24P-NN/v1-r1/day/siconca/v20250116/",
  "siconca_ANT-12_CESM2_historical_r11i1p1f1_UU-IMAU_RACMO24P-NN_v1-r1_day_19850101-19851231.nc"
)
racmo_lsm_path <- here("Data/ANT11_masks.nc")

# Index-space offset of the RACMO landmask relative to the RACMO sic
# template grid. A fixed property of the RACMO ANT-12 grid.
racmo_lsm_offset <- c(34, 30)

# Fixed bounding box (rotated-pole degrees) used to correct the RACMO sic
# template's extent before deriving the landmask offset. A fixed property
# of the RACMO ANT-12 grid, not of any particular file.
racmo_sic_template_ext <- ext(c(144, 210, -28.1, 25))

racmo_crsfix <- "GEOGCRS[\"Rotated_pole\",
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

# Reference layers for the sea-ice coast/SG masking + 60S crop + buffer step.
Ant_extent <- vect(here("Data/PolarRes26/add_data_limit_v7.2.shp"))  # EPSG:3031, 60 deg S boundary
coast      <- vect(here("Data/PolarRes26/add_coastline_high_res_polygon_v7_12.shp"))
SG         <- vect(here("Data/PolarRes26/orkney.shp"))

sea_ice_buffer_km <- c(2, 5, 10, 50, 100, 200)  # must match the buffer radii used in section 5e

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

# Every output path this job would produce -- used only for the
# skip-existing check. COASTLINE/ICEFREE are always produced; ICEFREE_FUTURE
# only for FUTURE-period files; CONCENTRATION/BUFFER only for sea ice (in
# addition to, not instead of, the COASTLINE/ICEFREE/ICEFREE_FUTURE set).
domain_mask_outputs <- c(coastline_path, icefree_path)
if (is_future) domain_mask_outputs <- c(domain_mask_outputs, icefree_futurepath)

concentration_path <- file.path(out_dir, paste0(stem, "_CONCENTRATION.tif"))
buffer_paths        <- file.path(out_dir, paste0(stem, "_BUFFER_", sea_ice_buffer_km, "km.tif"))

expected_outputs <- if (is_sea_ice) {
  c(concentration_path, buffer_paths)
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

# ---- 5a. CRS handling -------------------------------------------------------------

if (is.na(crs(r)) || crs(r) == "") {
  
  if (grepl("^HCLIM", model)) {
    
    # Known HCLIM quirk: some HCLIM outputs come out of Script 1/2 with no
    # CRS attached, even though their grid/extent is fine. For HCLIM only,
    # recover the CRS from a known-good HCLIM template file -- but only
    # after confirming the extents actually line up, so we never silently
    # stamp a wrong CRS onto a raster that doesn't actually match the
    # template grid.
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
    # RACMO (non-sea-ice) / MetUM / anything else: an empty CRS here is not
    # a known, safe-to-patch quirk -- treat it as the hard error it always
    # has been.
    stop("Input raster has no CRS: ", input_path,
         " -- reprojecting without a source CRS would silently produce garbage.",
         " Check this file (and Script 1/2's terra version, if this is a RACMO file).")
  }
}

if (grepl("^RACMO", model) && is_sea_ice) {
  
  # ---- 5b. RACMO sea-ice landmask fix ------------------------------------------
  # The template establishes crsfix/off once (always the same template,
  # regardless of which RACMO sea-ice job is running); the actual masking
  # is applied to THIS job's `r`, not the template.
  message("  RACMO sea-ice file -- applying rotated-pole CRS and landmask fix")
  
  sic_template <- rast(racmo_sic_template_path)
  set.ext(sic_template, racmo_sic_template_ext)
  set.crs(sic_template, racmo_crsfix)
  
  rs <- res(sic_template)
  
  lsm <- rast(racmo_lsm_path, "LSM", guessCRS = FALSE)
  x0 <- xmin(sic_template) - racmo_lsm_offset[1] * rs[1]
  y0 <- ymin(sic_template) - racmo_lsm_offset[2] * rs[2]
  ext(lsm) <- ext(x0, x0 + ncol(lsm) * rs[1], y0, y0 + nrow(lsm) * rs[2])
  crs(lsm) <- racmo_crsfix
  
  lsm_c <- crop(lsm, sic_template)
  if (!isTRUE(compareGeom(lsm_c, sic_template, stopOnError = FALSE))) {
    stop("RACMO landmask does not align with the sic template grid after ",
         "applying the fixed offset -- check racmo_lsm_offset / the template ",
         "file. Template: ", racmo_sic_template_path)
  }
  
  # Apply the fix to the ACTUAL job raster, not the template.
  crs(r) <- racmo_crsfix
  if (!isTRUE(compareGeom(r, lsm_c, stopOnError = FALSE))) {
    stop("This job's sea-ice raster does not share the RACMO sic template's ",
         "grid, so the landmask can't be applied directly -- an explicit ",
         "ext() fix (like the template's) may be needed for this file. ",
         "Input: ", input_path)
  }
  
  # Where lsm == 1 is land -- mask it out.
  r <- mask(r, lsm_c, maskvalues = 1)
  
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

# ---- 5c. Every non-sea-ice variable: generic COASTLINE / ICEFREE / ICEFREE_FUTURE --

if (!is_sea_ice) {
  
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
  
  message("Done: ", model, " / ", filename)
  quit(save = "no", status = 0)
}

# ---- 5d. Sea ice only: coastal gap-fill (non-interactive) -------------------------
# Bilinear reprojection can leave small NA gaps right at the coastline.
# Fixed threshold, no histogram, no manual tuning -- run the same way for
# every model.

message("  sea-ice file detected -- running gap-fill, coast/SG masking, 60S crop, and buffer products")

size_thresh <- 100

r <- subst(r, NA, -999)
r <- mask(r, coast_domain, inverse = TRUE)

missing <- r == -999
miss_r <- ifel(missing, 1, NA)
p <- patches(miss_r, directions = 8, zeroAsNA = TRUE)
patch_sizes <- freq(p)

# Patches touching the raster edge are a second, independent signal that a
# patch is out-of-domain rather than a local coastline artifact -- never
# fill these, however small.
edge_ids <- unique(na.omit(c(
  p[1, ],                 # top row
  p[nrow(p), ],           # bottom row
  p[, 1],                 # left col
  p[, ncol(p)]            # right col
)))

# Only small, interior patches are treated as genuinely fixable.
good_ids <- patch_sizes$value[
  !(patch_sizes$value %in% edge_ids) & patch_sizes$count <= size_thresh
]
fillable <- p %in% good_ids

# Build a working raster where ONLY the fillable cells are NA; true
# out-of-domain -999 stays excluded from the averaging so it can't leak
# into interpolated values.
work <- r
work[fillable] <- NA
work[missing & !fillable] <- NA

# Gap-fill iteratively with a focal mean that only touches NA cells.
filled <- work
w <- matrix(1, 3, 3)
for (i in 1:10) {
  still_na <- is.na(filled) & fillable
  if (!any(as.logical(values(still_na)), na.rm = TRUE)) break
  filled <- focal(filled, w = w, fun = "mean", na.rm = TRUE, na.policy = "only")
}

# Stitch back: replace only the fillable cells, leave everything else
# (including the true -999 domain edges) exactly as it was.
r_final <- r
r_final[fillable] <- filled[fillable]
r <- ifel(r_final == -999, NA, r_final)

# ---- 5e. Sea ice only: coast/SG masking, 60S crop, and buffer products ------------
# Ice-free-land layers DO change by period: historical & mid use the
# current domain, future uses the projected future domain. Uses the SAME
# ice_free_domain / ice_free_future_domain objects loaded once in section 1
# for every other variable, rather than re-reading separate copies from a
# different path -- so sea ice and every other variable are guaranteed to
# be using identical ice-free-land inputs.

sea_ice_domain <- ice_free_domain

domain.pts <- as.points(sea_ice_domain, values = TRUE)

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


# Zero-out the South Georgia / South America cells that sit at 100% SIC
# year-round -- this is a model-domain artefact, not real sea ice.
IDs <- terra::extract(conc, SG, cells = TRUE)
IDs <- IDs[!is.na(IDs[, 2]) & IDs[, 2] == 100, ]
if (nrow(IDs) > 0) conc[IDs$cell] <- 0

# 2. Cropped to the 60S extent, keeping every overlapping cell
conc_cropped <- terra::crop(conc, Ant_extent, snap = "out")
writeRaster(conc_cropped, concentration_path, overwrite = TRUE)
message("  wrote ", concentration_path)

# 3. Mean concentration within each buffer distance
for (i in seq_along(buffers)) {
  label <- names(buffers)[i]
  buffered <- extract_to_buffer(conc, sea_ice_domain, domain.pts, buffers[[label]])
  writeRaster(buffered, buffer_paths[i], overwrite = TRUE)
  message("  wrote ", buffer_paths[i])
}

message("Done: ", model, " / ", filename)