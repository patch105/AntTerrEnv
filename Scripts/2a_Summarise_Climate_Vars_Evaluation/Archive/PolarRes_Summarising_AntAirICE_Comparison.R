# HPC
lib_loc <- paste(getwd(),"/r_lib",sep="")

library(terra)
library(here)
library(arrow)
library(lubridate)


# extract the arguments provided in the command line
args <- commandArgs(trailingOnly = TRUE)
# The first argument is now the job index
job_index <- as.integer(args[1])

models <- list("MAR_MPI_ESM1", "MAR_CESM2", "HCLIM_MPI_ESM1", "HCLIM_CESM2")
# models <- list("HCLIM_CESM2")

model = models[[job_index]]
print(paste0("Model is: ", model))


if(model == "MAR_MPI_ESM1"){
  
  ##########################################
  # MAR MPI-ESM1 ------------------------------------------------------------
  ##########################################
  
  files <- list.files("/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes/MAR_MPI_ESM1_tif")
  
  file_paths <- list.files("/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes/MAR_MPI_ESM1_tif",
                           full.names = TRUE, recursive = TRUE)
  
  # Set the output directory
  outpath <- here("Data/Environmental_predictors/PolarRes/MAR_MPI_ESM1/Validation")
  
  
  tmp_dir <- tempdir()
  
  # Find TT all years
  variable_names <- files[grepl(pattern = "TT", files)]
  variable_paths <- file_paths[grepl(pattern = "TT", file_paths)]
  
  years_hist <- seq(2003, 2014, by = 1)
  months <- seq(1, 12, by = 1)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    # Calculate 8-hourly layer indices
    index_start <- (doy_start - 1) * 3 + 1
    index_end <- doy_end * 3
    
    return(seq(index_start, index_end))
  }
  
  # STEP 1: Calculate monthly means for each year
  # Store all monthly means organized by month across years
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name # month.name is a built in constant
  
  # For every year
  for(y in seq_along(years_hist)) {
    
    r <- terra::rast(variable_paths[grepl(variable_paths, pattern = years_hist[y])])
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_hist[y], months[m])
      
      # Subset for the month of interest (3 rasters per day)
      r_3hr <- r[[Doy]]
      
      #First, need to integrate over the time period to get daily sum
      
      # Number of 8-layer (daily) blocks
      n_days <- nlyr(r_3hr) / 3
      
      # Error check: should be divisible by 3
      if (nlyr(r_3hr) %% 3 != 0) stop("Layer count not divisible by 3")
      
      # Create an index that repeats each group of 3 layers
      index <- rep(1:n_days, each = 3)
      
      # Group every 3 layers and compute daily mean
      r_month <- terra::tapp(r_3hr, index = index, fun = "mean", na.rm = TRUE)
      
      # Take the monthly average from daily values
      r_month_mean <- app(r_month, mean, na.rm = TRUE)
      
      # Store in the appropriate month's list
      if(is.null(monthly_means_by_month[[m]])) {
        monthly_means_by_month[[m]] <- list()
      }
      monthly_means_by_month[[m]][[y]] <- r_month_mean
      
    }
    
    # Clean up temp files after each year
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 2: Calculate climatological monthly means (average across all years for each month)
  climatological_monthly_means <- list()
  
  for(m in seq_along(months)) {
    
    # Stack all years for this month
    month_stack <- rast(monthly_means_by_month[[m]])
    
    # Calculate mean across all years
    climatological_mean <- app(month_stack, mean, na.rm = TRUE)
    
    # Save the climatological monthly mean
    month_name <- sprintf("%02d", m)
    name <- paste0(outpath, "/Climatological_Monthly_Mean_Temperature_",
                   month.name[m], "_2003_2014.tif")
    writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
    
    # Store in list for annual calculation
    climatological_monthly_means[[m]] <- climatological_mean
    
    # Clean up temp files
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 3: Calculate mean annual temperature from the 12 climatological monthly means
  annual_stack <- rast(climatological_monthly_means)
  mean_annual_temp <- app(annual_stack, mean, na.rm = TRUE)
  
  # Save the final climatological mean annual temperature
  writeRaster(mean_annual_temp,
              paste0(outpath, "/Mean_Annual_Temperature_HISTORICAL_2003_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # 2. Summer temperature ---------------------------------------------------
  
  # Find TT all years
  variable_names <- files[grepl(pattern = "TT", files)]
  variable_paths <- file_paths[grepl(pattern = "TT", file_paths)]
  
  # MONTHLY VALUES HISTORICAL  ------------------------------------------------
  
  months = seq(1, 12, by=1)
  years_hist <- seq(2003, 2014, by = 1) # Include 1994 to get Dec for 1995
  
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    # Calculate 8-hourly layer indices
    index_start <- (doy_start - 1) * 3 + 1
    index_end <- doy_end * 3
    
    return(seq(index_start, index_end))
  }
  
  # STEP 1: Calculate monthly means for each year from daily values
  # Store all monthly means organized by month across years
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name
  
  for(y in seq_along(years_hist)) {
    
    r <- terra::rast(variable_paths[grepl(variable_paths, pattern = years_hist[y])])
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_hist[y], months[m])
      
      # Subset for the month of interest (3 rasters per day)
      r_3hr <- r[[Doy]]
      
      #First, need to integrate over the time period to get daily sum
      
      # Number of 8-layer (daily) blocks
      n_days <- nlyr(r_3hr) / 3
      
      # Error check: should be divisible by 3
      if (nlyr(r_3hr) %% 3 != 0) stop("Layer count not divisible by 3")
      
      # Create an index that repeats each group of 3 layers
      index <- rep(1:n_days, each = 3)
      
      # Group every 3 layers and compute daily mean
      r_month <- terra::tapp(r_3hr, index = index, fun = "mean", na.rm = TRUE)
      
      # Take the monthly average from daily values
      r_month_mean <- app(r_month, mean, na.rm = TRUE)
      
      # Store in the appropriate month's list
      if(is.null(monthly_means_by_month[[m]])) {
        monthly_means_by_month[[m]] <- list()
      }
      monthly_means_by_month[[m]][[y]] <- r_month_mean
      
    }
    
    # Clean up temp files after each year
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 2: Calculate climatological monthly means (average across all years for each month)
  
  # December 2003 - 2013 (for summers 2004 - 2014)
  # January 2004 - 2014 (for summers 2004 - 2014)
  # February 2004 - 2014 (for summers 2004 - 2014)
  
  ## JANUARY
  
  m = 1
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove first layer (January 2003)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_JAN <- app(month_stack, mean, na.rm = TRUE)
  
  ## FEBRUARY
  
  m = 2
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove first layer (February 2003)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_FEB <- app(month_stack, mean, na.rm = TRUE)
  
  ## DECEMBER
  
  m = 12
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove final layer (December *2014*)
  month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  
  # Calculate mean across all years
  climatological_mean_DEC <- app(month_stack, mean, na.rm = TRUE)
  
  
  # STEP 3: Calculate mean summer temperature from the 3 climatological monthly means
  # Summer = DJF (December, January, February)
  
  summer_stack <- c(climatological_mean_DEC, climatological_mean_JAN, climatological_mean_FEB)
  mean_summer_temp <- app(summer_stack, mean, na.rm = TRUE)
  
  # Save the climatological mean summer temperature
  writeRaster(mean_summer_temp,
              paste0(outpath, "/Mean_Summer_Temperature_HISTORICAL_2004_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # 3. Winter temperature ---------------------------------------------------
  
  # Find TT all years
  variable_names <- files[grepl(pattern = "TT", files)]
  variable_paths <- file_paths[grepl(pattern = "TT", file_paths)]
  
  # HISTORICAL WINTER TEMPERATURE (2003-2014) ------------------------------
  
  # Load the three climatological monthly means needed for winter(JJA)
  jun_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_June_2003_2014.tif"))
  jul_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_July_2003_2014.tif"))
  aug_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_August_2003_2014.tif"))
  
  # Calculate mean summer temperature from the 3 climatological monthly means
  winter_stack <- c(jun_clim, jul_clim, aug_clim)
  mean_winter_temp <- app(winter_stack, mean, na.rm = TRUE)
  
  writeRaster(mean_winter_temp,
              paste0(outpath, "/Mean_Winter_Temperature_HISTORICAL_2003_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
}

if(model == "MAR_CESM2"){ 
  
  ##########################################
  # MAR CESM2 ------------------------------------------------------------
  ##########################################
  
  files <- list.files("/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes/MAR_CESM2_tif")
  
  file_paths <- list.files("/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes/MAR_CESM2_tif",
                           full.names = TRUE, recursive = TRUE)
  
  # Set the output directory
  outpath <- here("Data/Environmental_predictors/PolarRes/MAR_CESM2/Validation")
  
  
  tmp_dir <- tempdir()
  
  # Find TT all years
  variable_names <- files[grepl(pattern = "TT", files)]
  variable_paths <- file_paths[grepl(pattern = "TT", file_paths)]
  
  years_hist <- seq(2003, 2014, by = 1)
  months <- seq(1, 12, by = 1)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    # Calculate 8-hourly layer indices
    index_start <- (doy_start - 1) * 3 + 1
    index_end <- doy_end * 3
    
    return(seq(index_start, index_end))
  }
  
  # STEP 1: Calculate monthly means for each year
  # Store all monthly means organized by month across years
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name # month.name is a built in constant
  
  # For every year
  for(y in seq_along(years_hist)) {
    
    r <- terra::rast(variable_paths[grepl(variable_paths, pattern = years_hist[y])])
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_hist[y], months[m])
      
      # Subset for the month of interest (3 rasters per day)
      r_3hr <- r[[Doy]]
      
      #First, need to integrate over the time period to get daily sum
      
      # Number of 8-layer (daily) blocks
      n_days <- nlyr(r_3hr) / 3
      
      # Error check: should be divisible by 3
      if (nlyr(r_3hr) %% 3 != 0) stop("Layer count not divisible by 3")
      
      # Create an index that repeats each group of 3 layers
      index <- rep(1:n_days, each = 3)
      
      # Group every 3 layers and compute daily mean
      r_month <- terra::tapp(r_3hr, index = index, fun = "mean", na.rm = TRUE)
      
      # Take the monthly average from daily values
      r_month_mean <- app(r_month, mean, na.rm = TRUE)
      
      # Store in the appropriate month's list
      if(is.null(monthly_means_by_month[[m]])) {
        monthly_means_by_month[[m]] <- list()
      }
      monthly_means_by_month[[m]][[y]] <- r_month_mean
      
    }
    
    # Clean up temp files after each year
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 2: Calculate climatological monthly means (average across all years for each month)
  climatological_monthly_means <- list()
  
  for(m in seq_along(months)) {
    
    # Stack all years for this month
    month_stack <- rast(monthly_means_by_month[[m]])
    
    # Calculate mean across all years
    climatological_mean <- app(month_stack, mean, na.rm = TRUE)
    
    # Save the climatological monthly mean
    month_name <- sprintf("%02d", m)
    name <- paste0(outpath, "/Climatological_Monthly_Mean_Temperature_",
                   month.name[m], "_2003_2014.tif")
    writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
    
    # Store in list for annual calculation
    climatological_monthly_means[[m]] <- climatological_mean
    
    # Clean up temp files
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 3: Calculate mean annual temperature from the 12 climatological monthly means
  annual_stack <- rast(climatological_monthly_means)
  mean_annual_temp <- app(annual_stack, mean, na.rm = TRUE)
  
  # Save the final climatological mean annual temperature
  writeRaster(mean_annual_temp,
              paste0(outpath, "/Mean_Annual_Temperature_HISTORICAL_2003_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # 2. Summer temperature ---------------------------------------------------
  
  # Find TT all years
  variable_names <- files[grepl(pattern = "TT", files)]
  variable_paths <- file_paths[grepl(pattern = "TT", file_paths)]
  
  # MONTHLY VALUES HISTORICAL  ------------------------------------------------
  
  months = seq(1, 12, by=1)
  years_hist <- seq(2003, 2014, by = 1) # Include 1994 to get Dec for 1995
  
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    # Calculate 8-hourly layer indices
    index_start <- (doy_start - 1) * 3 + 1
    index_end <- doy_end * 3
    
    return(seq(index_start, index_end))
  }
  
  # STEP 1: Calculate monthly means for each year from daily values
  # Store all monthly means organized by month across years
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name
  
  for(y in seq_along(years_hist)) {
    
    r <- terra::rast(variable_paths[grepl(variable_paths, pattern = years_hist[y])])
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_hist[y], months[m])
      
      # Subset for the month of interest (3 rasters per day)
      r_3hr <- r[[Doy]]
      
      #First, need to integrate over the time period to get daily sum
      
      # Number of 8-layer (daily) blocks
      n_days <- nlyr(r_3hr) / 3
      
      # Error check: should be divisible by 3
      if (nlyr(r_3hr) %% 3 != 0) stop("Layer count not divisible by 3")
      
      # Create an index that repeats each group of 3 layers
      index <- rep(1:n_days, each = 3)
      
      # Group every 3 layers and compute daily mean
      r_month <- terra::tapp(r_3hr, index = index, fun = "mean", na.rm = TRUE)
      
      # Take the monthly average from daily values
      r_month_mean <- app(r_month, mean, na.rm = TRUE)
      
      # Store in the appropriate month's list
      if(is.null(monthly_means_by_month[[m]])) {
        monthly_means_by_month[[m]] <- list()
      }
      monthly_means_by_month[[m]][[y]] <- r_month_mean
      
    }
    
    # Clean up temp files after each year
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 2: Calculate climatological monthly means (average across all years for each month)
  
  # December 2003 - 2013 (for summers 2004 - 2014)
  # January 2004 - 2014 (for summers 2004 - 2014)
  # February 2004 - 2014 (for summers 2004 - 2014)
  
  ## JANUARY
  
  m = 1
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove first layer (January 2003)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_JAN <- app(month_stack, mean, na.rm = TRUE)
  
  ## FEBRUARY
  
  m = 2
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove first layer (February 2003)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_FEB <- app(month_stack, mean, na.rm = TRUE)
  
  ## DECEMBER
  
  m = 12
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove final layer (December *2014*)
  month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  
  # Calculate mean across all years
  climatological_mean_DEC <- app(month_stack, mean, na.rm = TRUE)
  
  
  # STEP 3: Calculate mean summer temperature from the 3 climatological monthly means
  # Summer = DJF (December, January, February)
  
  summer_stack <- c(climatological_mean_DEC, climatological_mean_JAN, climatological_mean_FEB)
  mean_summer_temp <- app(summer_stack, mean, na.rm = TRUE)
  
  # Save the climatological mean summer temperature
  writeRaster(mean_summer_temp,
              paste0(outpath, "/Mean_Summer_Temperature_HISTORICAL_2004_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # 3. Winter temperature ---------------------------------------------------
  
  # Find TT all years
  variable_names <- files[grepl(pattern = "TT", files)]
  variable_paths <- file_paths[grepl(pattern = "TT", file_paths)]
  
  # HISTORICAL WINTER TEMPERATURE (2003-2014) ------------------------------
  
  # Load the three climatological monthly means needed for winter(JJA)
  jun_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_June_2003_2014.tif"))
  jul_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_July_2003_2014.tif"))
  aug_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_August_2003_2014.tif"))
  
  # Calculate mean summer temperature from the 3 climatological monthly means
  winter_stack <- c(jun_clim, jul_clim, aug_clim)
  mean_winter_temp <- app(winter_stack, mean, na.rm = TRUE)
  
  writeRaster(mean_winter_temp,
              paste0(outpath, "/Mean_Winter_Temperature_HISTORICAL_2003_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
}


if(model == "HCLIM_MPI_ESM1"){
  
  ##########################################
  # HCLIM MPI-ESM1 ------------------------------------------------------------
  ##########################################
  
  files <- list.files("/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes/HCLIM_MPI_ESM1_tif")
  
  file_paths <- list.files("/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes/HCLIM_MPI_ESM1_tif",
                           full.names = TRUE, recursive = TRUE)
  
  outpath <- here("Data/Environmental_predictors/PolarRes/HCLIM_MPI_ESM1/Validation")
  
  tmp_dir <- tempdir()
  
  # Find tas all years
  variable_names <- files[grepl(pattern = "tas", files)]
  variable_paths <- file_paths[grepl(pattern = "tas", file_paths)]
  
  # HISTORICAL --------------------------------------------------------------
  
  library(lubridate)
  
  years_hist <- seq(2003, 2014, by = 1)
  months <- seq(1, 12, by = 1)
  
  # Get the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    return(seq(doy_start, doy_end))
  }
  
  # STEP 1: Calculate monthly means for each year
  # Store all monthly means organized by month across years
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name # month.name is a built in constant
  
  for(y in seq_along(years_hist)) {
    
    # Get the raster for the year
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_hist[y]),
                                       x = variable_paths, perl = TRUE)]
    
    print(years_hist[y])
    print(length(year_files))
    
    r <- terra::rast(year_files)
    r
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_hist[y], months[m])
      
      # Subset for the month of interest (1 raster per day)
      r_month <- r[[Doy]]
      
      # Convert from Kelvin to Celsius
      r_month <- r_month - 273.15
      
      # Take the monthly average from daily values
      r_month_mean <- app(r_month, mean, na.rm = TRUE)
      
      # Store in the appropriate month's list
      if(is.null(monthly_means_by_month[[m]])) {
        monthly_means_by_month[[m]] <- list()
      }
      monthly_means_by_month[[m]][[y]] <- r_month_mean
      
    }
    
    # Clean up temp files after each year
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 2: Calculate climatological monthly means (average across all years for each month)
  climatological_monthly_means <- list()
  
  for(m in seq_along(months)) {
    
    # Stack all years for this month
    month_stack <- rast(monthly_means_by_month[[m]])
    
    # Calculate mean across all years
    climatological_mean <- app(month_stack, mean, na.rm = TRUE)
    
    # Save the climatological monthly mean
    month_name <- sprintf("%02d", m)
    name <- paste0(outpath, "/Climatological_Monthly_Mean_Temperature_",
                   month.name[m], "_2003_2014.tif")
    writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
    
    # Store in list for annual calculation
    climatological_monthly_means[[m]] <- climatological_mean
    
    # Clean up temp files
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 3: Calculate mean annual temperature from the 12 climatological monthly means
  annual_stack <- rast(climatological_monthly_means)
  mean_annual_temp <- app(annual_stack, mean, na.rm = TRUE)
  
  # Save the final climatological mean annual temperature
  writeRaster(mean_annual_temp,
              paste0(outpath, "/Mean_Annual_Temperature_HISTORICAL_2003_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # 2. Summer temperature ---------------------------------------------------
  
  # Find TT all years
  variable_names <- files[grepl(pattern = "tas", files)]
  variable_paths <- file_paths[grepl(pattern = "tas", file_paths)]
  
  library(lubridate)
  
  months <- seq(1, 12, by = 1)
  
  # Get the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    return(seq(doy_start, doy_end))
  }
  
  
  # HISTORICAL SUMMER TEMPERATURE (2004-2014) ------------------------------
  
  
  years_hist <- seq(2003, 2014, by = 1)  # Include 1994 to get Dec for 1995
  months = seq(1, 12, by=1)
  
  # STEP 1: Calculate monthly means for each year from daily values
  # Store all monthly means organized by month across years
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name
  
  for(y in seq_along(years_hist)) {
    
    # Get the raster for the year
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_hist[y]),
                                       x = variable_paths, perl = TRUE)]
    
    print(years_hist[y])
    print(length(year_files))
    
    r <- terra::rast(year_files)
    r
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_hist[y], months[m])
      
      # Subset for the month of interest (1 raster per day)
      r_month <- r[[Doy]]
      
      # Convert from Kelvin to Celsius
      r_month <- r_month - 273.15
      
      # Take the monthly average from daily values
      r_month_mean <- app(r_month, mean, na.rm = TRUE)
      
      # Store in the appropriate month's list
      if(is.null(monthly_means_by_month[[m]])) {
        monthly_means_by_month[[m]] <- list()
      }
      monthly_means_by_month[[m]][[y]] <- r_month_mean
      
    }
    
    # Clean up temp files after each year
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 2: Calculate climatological monthly means (average across all years for each month)
  
  # December 2003 - 2013 (for summers 2004 - 2014)
  # January 2004 - 2014 (for summers 2004 - 2014)
  # February 2004 - 2014 (for summers 2004 - 2014)
  
  ## JANUARY
  
  m = 1
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove first layer (January 2003)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_JAN <- app(month_stack, mean, na.rm = TRUE)
  
  ## FEBRUARY
  
  m = 2
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove first layer (February 2003)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_FEB <- app(month_stack, mean, na.rm = TRUE)
  
  ## DECEMBER
  
  m = 12
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove final layer (December *2014*)
  month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  
  # Calculate mean across all years
  climatological_mean_DEC <- app(month_stack, mean, na.rm = TRUE)
  
  
  
  # STEP 3: Calculate mean summer temperature from the 3 climatological monthly means
  # Summer = DJF (December, January, February)
  
  summer_stack <- c(climatological_mean_DEC, climatological_mean_JAN, climatological_mean_FEB)
  mean_summer_temp <- app(summer_stack, mean, na.rm = TRUE)
  
  # Save the climatological mean summer temperature
  writeRaster(mean_summer_temp,
              paste0(outpath, "/Mean_Summer_Temperature_HISTORICAL_2004_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # 3. Winter temperature ---------------------------------------------------
  
  
  # HISTORICAL WINTER TEMPERATURE (2003-2014) ------------------------------
  
  # Load the three climatological monthly means needed for winter(JJA)
  jun_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_June_2003_2014.tif"))
  jul_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_July_2003_2014.tif"))
  aug_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_August_2003_2014.tif"))
  
  # Calculate mean summer temperature from the 3 climatological monthly means
  winter_stack <- c(jun_clim, jul_clim, aug_clim)
  mean_winter_temp <- app(winter_stack, mean, na.rm = TRUE)
  
  writeRaster(mean_winter_temp,
              paste0(outpath, "/Mean_Winter_Temperature_HISTORICAL_2003_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
}


if(model == "HCLIM_CESM2"){
  
  ##########################################
  # HCLIM CESM2 ------------------------------------------------------------
  ##########################################
  
  files <- list.files("/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes/HCLIM_CESM2_tif")
  
  file_paths <- list.files("/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes/HCLIM_CESM2_tif",
                           full.names = TRUE, recursive = TRUE)
  
  outpath <- here("Data/Environmental_predictors/PolarRes/HCLIM_CESM2/Validation")
  
  tmp_dir <- tempdir()
  
  # Find tas all years
  variable_names <- files[grepl(pattern = "tas", files)]
  variable_paths <- file_paths[grepl(pattern = "tas", file_paths)]
  
  
  # HISTORICAL --------------------------------------------------------------
  
  library(lubridate)
  
  years_hist <- seq(2003, 2014, by = 1)
  months <- seq(1, 12, by = 1)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    # Calculate Hourly layer indices
    index_start <- (doy_start - 1) * 24 + 1
    index_end <- doy_end * 24
    
    return(seq(index_start, index_end))
  }
  
  # STEP 1: Calculate monthly means for each year
  # Store all monthly means organized by month across years
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name # month.name is a built in constant
  
  for(y in seq_along(years_hist)) {
    
    # Get the rasters for the year, then remove the final day
    
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_hist[y]),
                                       x = variable_paths, perl = TRUE)]
    year_files <- year_files[!grepl(pattern = paste0((years_hist[y] - 1), "12"), x = year_files)]
    
    print(years_hist[y])
    print(length(year_files))
    
    # Load each file, remove its last layer, then combine
    r_list <- lapply(year_files, function(file) {
      r_temp <- terra::rast(file)
      r_temp[[1:(terra::nlyr(r_temp) - 1)]]  # Remove last layer
    })
    
    # Combine all into one raster
    r <- terra::rast(r_list)
    r
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_hist[y],months[m])
      
      # Subset the relevant layers for the month
      r_1hr <- r[[Doy]]
      
      #First, need to integrate over the time period to get daily sum
      
      # Number of 24-layer (daily) blocks
      n_days <- nlyr(r_1hr) / 24
      
      # Error check: should be divisible by 24
      if (nlyr(r_1hr) %% 24 != 0) stop("Layer count not divisible by 24")
      
      # Create an index that repeats each group of 24 layers
      index <- rep(1:n_days, each = 24)
      
      # Group every 24 layers and compute daily mean
      r_month <- terra::tapp(r_1hr, index = index, fun = "mean", na.rm = TRUE)
      
      # Convert from Kelvin to Celsius
      r_month <- r_month - 273.15
      
      # Take the monthly average from daily values
      r_month_mean <- app(r_month, mean, na.rm = TRUE)
      
      # Store in the appropriate month's list
      if(is.null(monthly_means_by_month[[m]])) {
        monthly_means_by_month[[m]] <- list()
      }
      monthly_means_by_month[[m]][[y]] <- r_month_mean
      
      
    }
    
    # Clean up temp files after each year
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  
  # STEP 2: Calculate climatological monthly means (average across all years for each month)
  climatological_monthly_means <- list()
  
  for(m in seq_along(months)) {
    
    # Stack all years for this month
    month_stack <- rast(monthly_means_by_month[[m]])
    
    # Calculate mean across all years
    climatological_mean <- app(month_stack, mean, na.rm = TRUE)
    
    # Save the climatological monthly mean
    month_name <- sprintf("%02d", m)
    name <- paste0(outpath, "/Climatological_Monthly_Mean_Temperature_",
                   month.name[m], "_2003_2014.tif")
    writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
    
    # Store in list for annual calculation
    climatological_monthly_means[[m]] <- climatological_mean
    
    # Clean up temp files
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 3: Calculate mean annual temperature from the 12 climatological monthly means
  annual_stack <- rast(climatological_monthly_means)
  mean_annual_temp <- app(annual_stack, mean, na.rm = TRUE)
  
  # Save the final climatological mean annual temperature
  writeRaster(mean_annual_temp,
              paste0(outpath, "/Mean_Annual_Temperature_HISTORICAL_2003_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  
  # 2. Summer temperature ---------------------------------------------------
  
  years_hist <- seq(2003, 2014, by = 1)  # Include 1994 to get Dec for 1995
  months = seq(1, 12, by=1)
  
  # STEP 1: Calculate monthly means for each year from daily values
  # Store all monthly means organized by month across years
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name
  
  for(y in seq_along(years_hist)) {
    
    # Get the rasters for the year, then remove the final day
    
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_hist[y]),
                                       x = variable_paths, perl = TRUE)]
    year_files <- year_files[!grepl(pattern = paste0((years_hist[y] - 1), "12"), x = year_files)]
    
    print(years_hist[y])
    print(length(year_files))
    
    # Load each file, remove its last layer, then combine
    r_list <- lapply(year_files, function(file) {
      r_temp <- terra::rast(file)
      r_temp[[1:(terra::nlyr(r_temp) - 1)]]  # Remove last layer
    })
    
    # Combine all into one raster
    r <- terra::rast(r_list)
    r
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_hist[y], months[m])
      
      # Subset the relevant layers for the month
      r_1hr <- r[[Doy]]
      
      #First, need to integrate over the time period to get daily sum
      
      # Number of 24-layer (daily) blocks
      n_days <- nlyr(r_1hr) / 24
      
      # Error check: should be divisible by 24
      if (nlyr(r_1hr) %% 24 != 0) stop("Layer count not divisible by 24")
      
      # Create an index that repeats each group of 24 layers
      index <- rep(1:n_days, each = 24)
      
      # Group every 24 layers and compute daily mean
      r_month <- terra::tapp(r_1hr, index = index, fun = "mean", na.rm = TRUE)
      
      # Convert from Kelvin to Celsius
      r_month <- r_month - 273.15
      
      # Take the monthly average from daily values
      r_month_mean <- app(r_month, mean, na.rm = TRUE)
      
      # Store in the appropriate month's list
      if(is.null(monthly_means_by_month[[m]])) {
        monthly_means_by_month[[m]] <- list()
      }
      monthly_means_by_month[[m]][[y]] <- r_month_mean
      
    }
    
    # Clean up temp files after each year
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 2: Calculate climatological monthly means (average across all years for each month)
  
  # December 2003 - 2013 (for summers 2004 - 2014)
  # January 2004 - 2014 (for summers 2004 - 2014)
  # February 2004 - 2014 (for summers 2004 - 2014)
  
  ## JANUARY
  
  m = 1
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove first layer (January 2003)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_JAN <- app(month_stack, mean, na.rm = TRUE)
  
  ## FEBRUARY
  
  m = 2
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove first layer (February 2003)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_FEB <- app(month_stack, mean, na.rm = TRUE)
  
  ## DECEMBER
  
  m = 12
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove final layer (December *2014*)
  month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  
  # Calculate mean across all years
  climatological_mean_DEC <- app(month_stack, mean, na.rm = TRUE)
  
  
  
  # STEP 3: Calculate mean summer temperature from the 3 climatological monthly means
  # Summer = DJF (December, January, February)
  
  summer_stack <- c(climatological_mean_DEC, climatological_mean_JAN, climatological_mean_FEB)
  mean_summer_temp <- app(summer_stack, mean, na.rm = TRUE)
  
  # Save the climatological mean summer temperature
  writeRaster(mean_summer_temp,
              paste0(outpath, "/Mean_Summer_Temperature_HISTORICAL_2004_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # 3. Winter temperature ---------------------------------------------------
  
  
  # Load the three climatological monthly means needed for winter(JJA)
  jun_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_June_2003_2014.tif"))
  jul_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_July_2003_2014.tif"))
  aug_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_August_2003_2014.tif"))
  
  # Calculate mean summer temperature from the 3 climatological monthly means
  winter_stack <- c(jun_clim, jul_clim, aug_clim)
  mean_winter_temp <- app(winter_stack, mean, na.rm = TRUE)
  
  writeRaster(mean_winter_temp,
              paste0(outpath, "/Mean_Winter_Temperature_HISTORICAL_2003_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
}


