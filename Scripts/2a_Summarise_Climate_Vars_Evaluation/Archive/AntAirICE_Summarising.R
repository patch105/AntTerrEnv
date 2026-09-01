
# Making the 2003 to 2014 monthly climatology for AntAirICE

# Annual, Summer (DJF), and Winter (JJA)

# HPC - Single job per raster
lib_loc <- paste(getwd(),"/r_lib",sep="")

library(terra)
library(here)
library(arrow)
library(lubridate)

# Extract the arguments provided in the command line
args <- commandArgs(trailingOnly = TRUE)

# The first argument is the job index (1-based for SLURM arrays)
job_index <- as.integer(args[1])

data_type <- "daily"
data_type <- "monthly"

##############################
# DAILY VALUES ------------------------------------------------------------
##############################

if(data_type == "daily") {
  
  # -------------------------------------------------------------------------
  # File discovery
  # -------------------------------------------------------------------------
  variable_names <- list.files(
    "/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/AntAirICE",
    pattern = "\\.tif$"
  )
  
  variable_paths <- list.files(
    "/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/AntAirICE",
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
  
  
}




##############################
# MONTHLY CLIMATOLOGY VALUES -----------------------------------------------
##############################

if(data_type == "monthly") {
  
  
  
  
}

# 
# # Making all monthly climatologies ----------------------------------------
# 
# tmp_dir <- tempdir()
# 
# years_hist <- seq(2003, 2014, by = 1)
# months <- seq(1, 12, by = 1)
# 
# # Get the days range for each month (what day index is in that month)
# get_doy_range <- function(year, month) {
#   first_day <- ymd(paste(year, month, "01", sep = "-"))
#   last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
#   
#   doy_start <- yday(first_day)
#   doy_end <- yday(last_day)
#   
#   return(seq(doy_start, doy_end))
# }
# 
# # STEP 1: Calculate monthly means for each year
# # Store all monthly means organized by month across years
# monthly_means_by_month <- vector("list", 12)
# names(monthly_means_by_month) <- month.name # month.name is a built in constant
# 
# for(y in seq_along(years_hist)) {
#   
#   # Get the raster for the year
#   year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_hist[y]),
#                                      x = variable_paths, perl = TRUE)]
#   
#   print(years_hist[y])
#   print(length(year_files))
#   
#   r <- terra::rast(year_files)
#   r
#   
#   for(m in seq_along(months)) {
#     
#     Doy <- get_doy_range(years_hist[y], months[m])
#     
#     # Subset for the month of interest (1 raster per day)
#     r_month <- r[[Doy]]
#     
#     # Convert from Kelvin to Celsius
#     r_month <- r_month - 273.15
#     
#     # Take the monthly average from daily values
#     r_month_mean <- app(r_month, mean, na.rm = TRUE)
#     
#     # Store in the appropriate month's list
#     if(is.null(monthly_means_by_month[[m]])) {
#       monthly_means_by_month[[m]] <- list()
#     }
#     monthly_means_by_month[[m]][[y]] <- r_month_mean
#     
#   }
#   
#   # Clean up temp files after each year
#   tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
#   file.remove(tmp_files)
#   
# }
# 
# # STEP 2: Calculate climatological monthly means (average across all years for each month)
# climatological_monthly_means <- list()
# 
# for(m in seq_along(months)) {
#   
#   # Stack all years for this month
#   month_stack <- rast(monthly_means_by_month[[m]])
#   
#   # Calculate mean across all years
#   climatological_mean <- app(month_stack, mean, na.rm = TRUE)
#   
#   # Save the climatological monthly mean
#   month_name <- sprintf("%02d", m)
#   name <- paste0(outpath, "/Climatological_Monthly_Mean_Temperature_",
#                  month.name[m], "_1995_2014.tif")
#   writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
#   
#   # Store in list for annual calculation
#   climatological_monthly_means[[m]] <- climatological_mean
#   
#   # Clean up temp files
#   tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
#   file.remove(tmp_files)
#   
# }
# 
# # Mean ANNUAL summary -----------------------------------------------------
# 
# # STEP 3: Calculate mean annual temperature from the 12 climatological monthly means
# annual_stack <- rast(climatological_monthly_means)
# mean_annual_temp <- app(annual_stack, mean, na.rm = TRUE)
# 
# # Save the final climatological mean annual temperature
# writeRaster(mean_annual_temp,
#             paste0(outpath, "/Mean_Annual_Temperature_HISTORICAL_1995_2014.tif"),
#             gdal = c("COMPRESS=NONE"), overwrite = TRUE)
# 
# 
# 
# # Mean SUMMER summary  ----------------------------------------------------
# 
# 
# 
# 
# # Mean WINTER summary -----------------------------------------------------


