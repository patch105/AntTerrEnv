# ==============================================================================
# PolarRes26 -- STEP 3 (COMPARISON SUBSET): REPROJECT, RESAMPLE, AND MASK
# ==============================================================================
# Same grid-adjustment logic as the full Step 3 script, but scoped to the
# temperature-only "comparison" outputs (HISTORICAL/FUTURE mean annual,
# monthly climatologies, summer, winter) produced by the trimmed evaluation
# version of Script 2. Reads from each model's
# Data/Environmental_predictors/PolarRes26/<model>/comparison folder, and
# writes to a matching .../Regridded/<model>/comparison folder.
#
# There's no sea-ice handling here -- these are temperature files only, so
# every input just gets the generic COASTLINE/ICEFREE masked outputs, same
# as the non-sea-ice branch of the full Step 3 script.
#
# One job = one input file, selected via a single job_index (same pattern
# as the full script) -- run with no argument (or an out-of-range one) to
# see the full job table.
# ==============================================================================

# ---- 0. Setup ------------------------------------------------------------------

lib_loc <- paste(getwd(), "/r_lib_new", sep = "")

library(dplyr, lib.loc = lib_loc)
library(purrr, lib.loc = lib_loc)
library(terra)
library(here)

# ---- 1. Configuration -----------------------------------------------------------

models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5",
            "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5",
            "MetUM_ERA5")

input_base  <- here("Data/Environmental_predictors/PolarRes26")
output_base <- here("Data/Environmental_predictors/PolarRes26/Regridded")

coast_domain    <- rast(here("Data/PolarRes26/coast_domain.tif"))
ice_free_domain <- rast(here("Data/PolarRes26/ice_free_domain.tif"))

# If TRUE, skip an input file entirely (no reprojection/resampling work at
# all) when every output it would produce already exists on disk. Set to
# FALSE to always reprocess and overwrite.
SKIP_EXISTING <- TRUE

# ---- 2. Build the job table: every .tif in each model's "comparison" folder ----

job_table <- map_dfr(models, function(model) {
  model_dir <- file.path(input_base, model, "comparison")
  if (!dir.exists(model_dir)) {
    message("Comparison directory not found, skipping: ", model_dir)
    return(tibble())
  }
  files <- list.files(model_dir, pattern = "\\.tif$", full.names = TRUE)
  tibble(model = model, input_path = files, filename = basename(files))
})

if (nrow(job_table) == 0) {
  stop("No .tif files found under any <model>/comparison folder under ", input_base,
       " -- has the comparison version of Script 2 been run yet?")
}

message(nrow(job_table), " output file(s) found across ",
        length(unique(job_table$model)), " model(s).")

# ---- 3. job_index selects ONE input file -----------------------------------------

args <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])

this_job   <- job_table[job_index, ]
model      <- this_job$model
input_path <- this_job$input_path
filename   <- this_job$filename
message("Job ", job_index, "/", nrow(job_table), " -> ", model, " / ", filename)

# ---- 4. Output paths -------------------------------------------------------------
# Same "comparison" subfolder, but under Regridded/<model>/ instead of
# <model>/ directly.

out_dir <- file.path(output_base, model, "comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stem <- tools::file_path_sans_ext(filename)
coastline_path <- file.path(out_dir, paste0(stem, "_COASTLINE.tif"))
icefree_path   <- file.path(out_dir, paste0(stem, "_ICEFREE.tif"))

expected_outputs <- c(coastline_path, icefree_path)

# ---- 4b. Skip-existing check ------------------------------------------------------

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

# Step 2: Reproject into the domain's CRS with bilinear interpolation
r <- project(r, target_grid, method = "bilinear")

# Step 3: Resample to the 1km domain grid (already in the right CRS &
# extent, so "near" here is just chopping up 10km to 1km grids)
r <- resample(r, coast_domain, method = "near")

# ---- 6. Mask, save ----------------------------------------------------------------
# No sea-ice special case here -- these are temperature-only comparison
# outputs, so every file just gets the two generic masked products.

r_coast <- mask(r, coast_domain, maskvalue = NA)
writeRaster(r_coast, coastline_path, overwrite = TRUE)
message("  wrote ", coastline_path)

r_icefree <- mask(r, ice_free_domain, maskvalue = NA)
writeRaster(r_icefree, icefree_path, overwrite = TRUE)
message("  wrote ", icefree_path)

message("Done: ", model, " / ", filename)