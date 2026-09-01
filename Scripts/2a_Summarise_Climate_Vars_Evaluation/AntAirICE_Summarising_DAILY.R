# Making the 2003 to 2014 daily ice-free values for AntAirICE
# HPC - Single job per raster

lib_loc <- paste(getwd(),"/r_lib_new",sep="")

library(terra)
library(here)
library(arrow)
library(lubridate)

# Extract the arguments provided in the command line
args <- commandArgs(trailingOnly = TRUE)

# The first argument is the job index (1-based for SLURM arrays)
job_index <- as.integer(args[1])

# -------------------------------------------------------------------------
# File discovery
# -------------------------------------------------------------------------
variable_names <- list.files(
  "/mnt/hpccs01/home/patterc2/n11222026/AntarcticFutureHabitat/Data/AntAirICE",
  pattern = "\\.tif$"
)

variable_paths <- list.files(
  "/mnt/hpccs01/home/patterc2/n11222026/AntarcticFutureHabitat/Data/AntAirICE",
  pattern   = "\\.tif$",
  full.names = TRUE,
  recursive  = TRUE
)

# Use job_index to select the file for this job
file_path <- variable_paths[job_index]
file_name <- variable_names[job_index]

# Set the output directory
outpath <- here("Data/AntAirICE/Summarised")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load domain (ice-free mask) — adjust path as needed
# -------------------------------------------------------------------------
domain <- rast(here("Data/Environmental_predictors/ice_free_upsamp_1km.tif"))

# Set domain values
domain <- ifel(not.na(domain), 1, NA)

# -------------------------------------------------------------------------
# Daily values: 2003–2014 ice-free Antarctic air temperature
# -------------------------------------------------------------------------

# 1. Load — get only layer 1 per raster
r <- rast(file_path)[[1]]

# 2. Reproject to Antarctic Polar Stereographic (EPSG:3031)
r <- project(r, "epsg:3031")

# 3. Reproject/resample to match the ice-free domain extent & resolution
#    (bilinear interpolation to correct for new domain)
r <- project(r, domain, method = "bilinear")

# 4. Crop and mask to ice-free areas
r <- mask(r, domain, maskvalue = NA)

# 5. Scale: multiply by 0.1
r <- r * 0.1

# -------------------------------------------------------------------------
# Save output
# -------------------------------------------------------------------------
# Strip .tif, append _ICEFREE.tif
out_name <- paste0(tools::file_path_sans_ext(file_name), "_ICEFREE.tif")
out_file  <- file.path(outpath, out_name)

writeRaster(r, out_file, overwrite = TRUE)

message("Saved: ", out_file)