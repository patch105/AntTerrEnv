# ==============================================================================
# PolarRes26 -- STEP 3: REPROJECT, RESAMPLE, AND MASK EVERY CLIMATOLOGY OUTPUT
# ==============================================================================
# Takes every .tif Script 2 wrote (for every model) and produces two masked
# versions of each, on a common grid:
#   - <name>_COASTLINE.tif -- masked to the coastline domain
#   - <name>_ICEFREE.tif   -- masked to the ice-free domain
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

# Mask to the coastline domain
r_coast <- mask(r, coast_domain, maskvalue = NA)
writeRaster(r_coast, coastline_path, overwrite = TRUE)
message("  wrote ", coastline_path)

# Mask to the ice-free domain
r_icefree <- mask(r, ice_free_domain, maskvalue = NA)
writeRaster(r_icefree, icefree_path, overwrite = TRUE)
message("  wrote ", icefree_path)

message("Done: ", model, " / ", filename)