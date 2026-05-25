
## NEW CODE TO ADD IN: 
# year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_future[y]), 
#                                    x = variable_paths, perl = TRUE)]

# Non-HPC
# lib_loc <- .libPaths() 

# HPC
lib_loc <- paste(getwd(),"/r_lib",sep="")

library(terra)
library(here)
library(arrow)
library(lubridate)


#files <- list.files("Z:/AntarcticFutureHabitat/Data/PolarRes/HCLIM_CESM2_tif")
files <- list.files("/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes/HCLIM_CESM2_tif")

#file_paths <- list.files("Z:/AntarcticFutureHabitat/Data/PolarRes/HCLIM_CESM2_tif",
#full.names = TRUE, recursive = TRUE)
file_paths <- list.files("/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes/HCLIM_CESM2_tif",
                         full.names = TRUE, recursive = TRUE)

# Set the output directory
# outpath <- "Z:/AntarcticFutureHabitat/Data/Environmental_predictors/PolarRes/HCLIM_CESM2"
outpath <- here("Data/Environmental_predictors/PolarRes/HCLIM_CESM2")


tmp_dir <- tempdir()

# extract the arguments provided in the command line
args <- commandArgs(trailingOnly = TRUE)
# The first argument is now the job index
job_index <- as.integer(args[1])

# Set the variables
variables <- list("temp", "total_DD", "wind", "sea_ice_nov", "sea_ice_oct_feb", "total_precip", "total_summer_precip", "mean_precip", "mean_summer_precip", "solar_rad", "mean_melt", "total_melt", "mean_snow", "summer_snow")

# variables <- list("temp", "wind", "mean_snow", "summer_snow")

variable = variables[[job_index]]
print(paste0("Variable is: ", variable))


# 1. TEMPERATURE  (originally in Kelvin) -----------------------------------

if(variable == "temp"){
  
  # Find tas all years
  variable_names <- files[grepl(pattern = "tas", files)]
  variable_paths <- file_paths[grepl(pattern = "tas", file_paths)]
  
  
  # HISTORICAL --------------------------------------------------------------
  
  library(lubridate)
  
  years_hist <- seq(1995, 2014, by = 1)
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
                   month.name[m], "_1995_2014.tif")
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
              paste0(outpath, "/Mean_Annual_Temperature_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # FUTURE CLIMATOLOGY (2081-2100) -----------------------------------------
  
  years_future <- seq(2081, 2100, by = 1)
  months <- seq(1, 12, by = 1)
  library(lubridate)
  
  # STEP 1: Calculate monthly means for each year
  # Store all monthly means organized by month across years
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name # month.name is a built in constant
  
  for(y in seq_along(years_future)) {
    
    # Get the rasters for the year, then remove the final day
    
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_future[y]),
                                       x = variable_paths, perl = TRUE)]
    
    print(years_future[y])
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
      
      Doy <- get_doy_range(years_future[y],months[m])
      
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
                   month.name[m], "_2081_2100.tif")
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
              paste0(outpath, "/Mean_Annual_Temperature_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # 2. MEAN SUMMER TEMPERATURE ----------------------------------------------
  # **ALSO** there was two corrupt file which weren't able to be downloaded: tas_fp_ANT11_ANT11_ALADIN43_v1_CESM2_r11i1p1f1_historical_1hr_199406010000-199407010000.nc
  # JUNE 1994 (30 days)

  # Find TT all years
  variable_names <- files[grepl(pattern = "tas", files)]
  variable_paths <- file_paths[grepl(pattern = "tas", file_paths)]

  library(lubridate)

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


  # HISTORICAL SUMMER TEMPERATURE (1995-2014) ------------------------------
  # Have to re-calculate to include December 1994

  years_hist <- seq(1994, 2014, by = 1)  # Include 1994 to get Dec for 1995
  months = seq(1, 12, by=1)

  # STEP 1: Calculate monthly means for each year from daily values
  # Store all monthly means organized by month across years
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name

  for(y in seq_along(years_hist)) {

    # Get the rasters for the year, then remove the final day

    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_hist[y]),
                                       x = variable_paths, perl = TRUE)]
    
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
    
    # ---- INSERT NA filler for missing rsus June 1994 (31 days) ----
    # This creates an NA raster for the missing month to preserve the right number of days
    if(years_hist[y] == 1994) {
      june_doys <- get_doy_range(1994, 6)           # days 152-181 (30 days)
      
      filler_june <- rast(replicate(length(june_doys),  # 30 NA layers matching grid
                                    init(r[[1]], fun = NA)))
      
      # Split r at the June insertion point, then reassemble
      r <- rast(list(
        r[[1:(june_doys[1] - 1)]],   # days before June
        filler_june,                            # NA placeholder for June
        r[[june_doys[1]:nlyr(r)]] # days after June
      ))
    }
    # ---------------------------------------------------------------
    
    for(m in seq_along(months)) {

      # ---- SKIP June 1994 entirely (data missing for that month) ----
      
      if(years_hist[y] == 1994 && m == 6) {
        message("Skipping June 1994 - corrupt file")
        monthly_means_by_month[[6]][[y]] <- NA  # NA preserves the slot
        next
      }
      # -------------------------------------------------------------------
      
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

  # December 1994 - 2013 (for summers 1995 - 2014)
  # January 1995 - 2014 (for summers 1995 - 2014)
  # February 1995 - 2014 (for summers 1995 - 2014)

  ## JANUARY

  m = 1
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])

  # Remove first layer (January 1994)
  month_stack <- month_stack[[2:nlyr(month_stack)]]

  # Calculate mean across all years
  climatological_mean_JAN <- app(month_stack, mean, na.rm = TRUE)

  ## FEBRUARY

  m = 2
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])

  # Remove first layer (February 1994)
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
              paste0(outpath, "/Mean_Summer_Temperature_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)


  # FUTURE SUMMER TEMPERATURE (2081-2099) ----------------------------------
  # Have to re-calculate to include December 2080

  years_future <- seq(2080, 2100, by = 1)  # Include 2080 to get Dec for 2081

  # STEP 1: Calculate monthly means for each year from daily values
  # Store all monthly means organized by month across years
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name

  for(y in seq_along(years_future)) {

    # Get the rasters for the year, then remove the final day

    # Get the raster for the year
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_future[y]),
                                       x = variable_paths, perl = TRUE)]
    
    print(years_future[y])
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

      Doy <- get_doy_range(years_future[y], months[m])

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
  # December 2080 - 2099 (for summers 2081 - 2100)
  # January 2081 - 2100 (for summers 2081 - 2100)
  # February 2081 - 2100 (for summers 2081 - 2100)


  ## JANUARY

  m = 1
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])

  # # Remove final layer (January 2100)
  # month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  # Remove first layer (January 2080)
  month_stack <- month_stack[[2:nlyr(month_stack)]]

  # Calculate mean across all years
  climatological_mean_JAN <- app(month_stack, mean, na.rm = TRUE)

  ## FEBRUARY

  m = 2
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])

  # # Remove final layer (February 2100)
  # month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  # Remove first layer (February 2080)
  month_stack <- month_stack[[2:nlyr(month_stack)]]

  # Calculate mean across all years
  climatological_mean_FEB <- app(month_stack, mean, na.rm = TRUE)

  ## DECEMBER

  m = 12
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])

  # Remove final layer (December *2100*)
  month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  # # Remove other final layer (December *2013*)
  # month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]

  # Calculate mean across all years
  climatological_mean_DEC <- app(month_stack, mean, na.rm = TRUE)


  # STEP 3: Calculate mean summer snow cover from the 3 climatological monthly means
  # Summer = DJF (December, January, February)

  summer_stack <- c(climatological_mean_DEC, climatological_mean_JAN, climatological_mean_FEB)
  mean_summer_temp_future <- app(summer_stack, mean, na.rm = TRUE)


  # Save the climatological mean summer temperature
  writeRaster(mean_summer_temp_future,
              paste0(outpath, "/Mean_Summer_Temperature_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  
  # 2B. MEAN WINTER TEMPERATURE ----------------------------------------------

  # HISTORICAL WINTER TEMPERATURE (1995-2014) ------------------------------

  # Load the three climatological monthly means needed for winter(JJA)
  jun_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_June_1995_2014.tif"))
  jul_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_July_1995_2014.tif"))
  aug_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_August_1995_2014.tif"))

  # Calculate mean summer temperature from the 3 climatological monthly means
  winter_stack <- c(jun_clim, jul_clim, aug_clim)
  mean_winter_temp <- app(winter_stack, mean, na.rm = TRUE)

  writeRaster(mean_winter_temp,
              paste0(outpath, "/Mean_Winter_Temperature_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)


  # FUTURE WINTER TEMPERATURE (2081-2100) ------------------------------

  # Load the three climatological monthly means needed for winter(JJA)
  jun_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_June_2081_2100.tif"))
  jul_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_July_2081_2100.tif"))
  aug_clim <- rast(paste0(outpath, "/Climatological_Monthly_Mean_Temperature_August_2081_2100.tif"))

  # Calculate mean summer temperature from the 3 climatological monthly means
  winter_stack <- c(jun_clim, jul_clim, aug_clim)
  mean_winter_temp <- app(winter_stack, mean, na.rm = TRUE)

  writeRaster(mean_winter_temp,
              paste0(outpath, "/Mean_Winter_Temperature_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
}




# 3. TOTAL POSITIVE DEGREE DAYS above -5 ----------------------------------------------------

# Find TT all years

if(variable == "total_DD"){
  
  variable_names <- files[grepl(pattern = "tas", files)]
  variable_paths <- file_paths[grepl(pattern = "tas", file_paths)]
  
  library(lubridate)
  
  months <- seq(1, 12, by = 1)
  
  
  # HISTORICAL --------------------------------------------------------------
  
  
  # For every year, calculate the mean
  
  years_hist <- seq(2009, 2014, by = 1)
  
  # Positive degree-days sum (above -5) for the year
  
  # Calculated from the daily values
  
  for(y in seq_along(years_hist)) {
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month
    
    # Get the rasters for the year, but exclude the file that spans from previous year (because of naming)
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_hist[y]),
                                       x = variable_paths, perl = TRUE)]
    year_files <- year_files[!grepl(pattern = paste0((years_hist[y] - 1), "12"), x = year_files)]
    
    print(years_hist[y])
    print(length(year_files))
    
    
    # Load each file, remove its last layer (first hour of next month), then combine
    r_list <- lapply(year_files, function(file) {
      r_temp <- terra::rast(file)
      r_temp[[1:(terra::nlyr(r_temp) - 1)]]  # Remove last layer
    })
    
    
    # Combine all into one raster
    r <- terra::rast(r_list)
    r
    
    # Error check: total layers should be divisible by 24
    if (terra::nlyr(r) %% 24 != 0) stop(paste("Year", years_hist[y], ": layer count not divisible by 24"))
    
    # Total number of days in the raster
    n_days <- terra::nlyr(r) / 24
    
    # Index grouping every 24 hourly layers into one day
    index <- rep(1:n_days, each = 24)
    
    # Collapse hourly -> daily by taking the mean across each 24-layer block
    r_daily <- terra::tapp(r, index = index, fun = "mean", na.rm = TRUE)
    
    # Convert from Kelvin to Celcius
    r_daily <- r_daily - 273.15
    
    limit <- -5 # -5 degrees celsius
    
    # For every cell for every day, was it above -5? If so, set to its temperature above the limit
    r_daily <- ifel(r_daily > limit, r_daily - limit, 0)
    
    # Sum these to get the yearly positive degree-day sum:
    yearly_pdd_sum <- app(r_daily, sum, na.rm = TRUE)
    
    # Save the yearly total
    name <- paste0(outpath, "/Total_Annual_Degree_Days-5_Year_", years_hist[y],".tif")
    
    writeRaster(yearly_pdd_sum, name, gdal=c("COMPRESS=NONE"), overwrite = T)
    
    tmp_files <- list.files(tmp_dir, full.names = T, pattern = "^file")
    
    file.remove(tmp_files)
    
  }
  
  
  # ANNUAL MEAN OF TOTAL DEGREE DAYS FOR ALL YEARS (Historical) -----------------
  
  # Update years hist
  years_hist <- seq(1995, 2014, by = 1)
  
  # List all the yearly mean rasters you just saved
  annual_mean_files <- list.files(outpath, pattern = "^Total_Annual_Degree_Days-5_Year_.*\\.tif$", full.names = TRUE)
  
  # Subset just the relevant ones based on year:
  # Collapse the years into a regex pattern: "1995|1996|...|2014"
  year_pattern <- paste0("(", paste(years_hist, collapse = "|"), ")")
  
  # Subset the files
  annual_mean_files <- annual_mean_files[grepl(year_pattern, annual_mean_files)]
  
  # Stack them together
  annual_means <- rast(annual_mean_files)
  
  # Calculate the overall mean
  final_mean <- app(annual_means, mean, na.rm = TRUE)
  
  # Save the final mean raster
  writeRaster(final_mean, paste0(outpath, "/Mean_Annual_Total_Degree_Days-5_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # FUTURE --------------------------------------------------------------

  # Find TT all years
  variable_names <- files[grepl(pattern = "tas", files)]
  variable_paths <- file_paths[grepl(pattern = "tas", file_paths)]

  library(lubridate)

  months <- seq(1, 12, by = 1)


  # For every year, calculate the mean

  years_future <- seq(2081, 2100, by = 1)


  # Positive degree-days sum (above -5) for the year

  # Calculated from the daily values

  for(y in seq_along(years_future)) {

    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month

    # Get the rasters for the year, but exclude the file that spans from previous year (because of naming)
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_future[y]),
                                       x = variable_paths, perl = TRUE)]
    year_files <- year_files[!grepl(pattern = paste0((years_future[y] - 1), "12"), x = year_files)]


    print(years_future[y])
    print(length(year_files))
    
    # Load each file, remove its last layer (first hour of next month), then combine
    r_list <- lapply(year_files, function(file) {
      r_temp <- terra::rast(file)
      r_temp[[1:(terra::nlyr(r_temp) - 1)]]  # Remove last layer
    })


    # Combine all into one raster
    r <- terra::rast(r_list)
    r
    
    # Error check: total layers should be divisible by 24
    if (terra::nlyr(r) %% 24 != 0) stop(paste("Year", years_future[y], ": layer count not divisible by 24"))

    # Total number of days in the raster
    n_days <- terra::nlyr(r) / 24

    # Index grouping every 24 hourly layers into one day
    index <- rep(1:n_days, each = 24)

    # Collapse hourly -> daily by taking the mean across each 24-layer block
    r_daily <- terra::tapp(r, index = index, fun = "mean", na.rm = TRUE)

    # Convert from Kelvin to Celcius
    r_daily <- r_daily - 273.15

    limit <- -5 # -5 degrees celsius

    # For every cell for every day, was it above -5? If so, set to its temperature above the limit
    r_daily <- ifel(r_daily > limit, r_daily - limit, 0)

    # Sum these to get the yearly positive degree-day sum:
    yearly_pdd_sum <- app(r_daily, sum, na.rm = TRUE)

    # Save the yearly total
    name <- paste0(outpath, "/Total_Annual_Degree_Days-5_Year_", years_future[y],".tif")

    writeRaster(yearly_pdd_sum, name, gdal=c("COMPRESS=NONE"), overwrite = T)

    tmp_files <- list.files(tmp_dir, full.names = T, pattern = "^file")

    file.remove(tmp_files)

  }


  # ANNUAL MEAN OF TOTAL DEGREE DAYS FOR ALL YEARS (Future) -----------------

  # List all the yearly mean rasters you just saved
  annual_mean_files <- list.files(outpath, pattern = "^Total_Annual_Degree_Days-5_Year_.*\\.tif$", full.names = TRUE)

  # Subset just the relevant ones based on year:
  # Collapse the years into a regex pattern: "2081|2082|...|2100"
  year_pattern <- paste0("(", paste(years_future, collapse = "|"), ")")

  # Subset the files
  annual_mean_files <- annual_mean_files[grepl(year_pattern, annual_mean_files)]

  # Stack them together
  annual_means <- rast(annual_mean_files)

  # Calculate the overall mean
  final_mean <- app(annual_means, mean, na.rm = TRUE)

  # Save the final mean raster
  writeRaster(final_mean, paste0(outpath, "/Mean_Annual_Total_Degree_Days-5_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
}



# 4. WIND (WIND SPEED) (m s-1)  -------------------------------------------------------------
if(variable == "wind"){
  
  variable_names <- files[grepl(pattern = "sfcWind", files)]
  variable_paths <- file_paths[grepl(pattern = "sfcWind", file_paths)]
  
  # These are stored in a MONTHLY format, with *hourly* values per month
  # Each monthly raster includes the first *hour* of the next month, which I deal with
  # So there's for e.g. 745 hourly values rather than 744 for January 1994
  
  # Since all months have downloaded, I don't need to do any other fiddling around for corrupt files
  
  
  # HISTORICAL --------------------------------------------------------------
  
  # For every year, pull out the monthly rasters, calculate their mean, then summarise for entire year
  
  years_hist <- seq(1995, 2014, by = 1)
  months = seq(1, 12, by=1)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    # Calculate 6-hourly layer indices
    index_start <- (doy_start - 1) * 24 + 1
    index_end <- doy_end * 24
    
    return(seq(index_start, index_end))
  }
  
  # STEP 1: Calculate monthly means for each year
  # Store all monthly means organized by month across years
  
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name # month.name is a built in constant
  
  
  for(y in seq_along(years_hist)) {
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month
    
    # Get the rasters for the year, but exclude the file that spans from previous year (because of naming)
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_hist[y]),
                                       x = variable_paths, perl = TRUE)]
    year_files <- year_files[!grepl(pattern = paste0((years_hist[y] - 1), "12"), x = year_files)]
    
    print(years_hist[y])
    print(length(year_files))
    
    
    # Load each file, remove its last layer (first hour of next month), then combine
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
    name <- paste0(outpath, "/Climatological_Monthly_Mean_Wind_Speed_",
                   month.name[m], "_1995_2014.tif")
    writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
    
    # Store in list for annual calculation
    climatological_monthly_means[[m]] <- climatological_mean
    
    # Clean up temp files
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 3: Calculate mean annual wind speed from the 12 climatological monthly means
  annual_stack <- rast(climatological_monthly_means)
  mean_annual_temp <- app(annual_stack, mean, na.rm = TRUE)
  
  # Save the final climatological mean annual temperature
  writeRaster(mean_annual_temp,
              paste0(outpath, "/Mean_Annual_Wind_Speed_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # FUTURE CLIMATOLOGY (2081 - 2100) -------------------------------------------
  
  # For every year, pull out the monthly rasters, calculate their mean, then summarise for entire year
  
  years_future <- seq(2081, 2100, by = 1)
  months = seq(1, 12, by=1)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    # Calculate 6-hourly layer indices
    index_start <- (doy_start - 1) * 24 + 1
    index_end <- doy_end * 24
    
    return(seq(index_start, index_end))
  }
  
  # STEP 1: Calculate monthly means for each year
  # Store all monthly means organized by month across years
  
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name # month.name is a built in constant
  
  
  for(y in seq_along(years_future)) {
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month
    
    # Get the rasters for the year, but exclude the file that spans from previous year (because of naming)
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_future[y]),
                                       x = variable_paths, perl = TRUE)]
    year_files <- year_files[!grepl(pattern = paste0((years_future[y] - 1), "12"), x = year_files)]
    
    print(years_future[y])
    print(length(year_files))
    
    # Load each file, remove its last layer (first hour of next month), then combine
    r_list <- lapply(year_files, function(file) {
      r_temp <- terra::rast(file)
      r_temp[[1:(terra::nlyr(r_temp) - 1)]]  # Remove last layer
    })
    
    
    # Combine all into one raster
    r <- terra::rast(r_list)
    r
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_future[y],months[m])
      
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
    name <- paste0(outpath, "/Climatological_Monthly_Mean_Wind_Speed_",
                   month.name[m], "_2081_2100.tif")
    writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
    
    # Store in list for annual calculation
    climatological_monthly_means[[m]] <- climatological_mean
    
    # Clean up temp files
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 3: Calculate mean annual wind speed from the 12 climatological monthly means
  annual_stack <- rast(climatological_monthly_means)
  mean_annual_temp <- app(annual_stack, mean, na.rm = TRUE)
  
  # Save the final climatological mean annual temperature
  writeRaster(mean_annual_temp,
              paste0(outpath, "/Mean_Annual_Wind_Speed_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
}




# 5. SEA ICE SICONCA (%) --------------------------------------------------
# These are stored in a MONTHLY format, with daily values per month
# Each monthly raster includes the *first day* of the next month, which I deal with

if(variable == "sea_ice_nov"){
  
  # Find sea ice
  variable_names <- files[grepl(pattern = "siconca", files)]
  variable_paths <- file_paths[grepl(pattern = "siconca", file_paths)]
  
  # HISTORICAL --------------------------------------------------------------
  
  # For every year, pull out the monthly rasters, calculate their mean, then summarise for entire year
  
  years_hist <- seq(1995, 2014, by = 1)
  
  months = seq(1, 12, by=1)
  
  library(lubridate)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    return(seq(doy_start, doy_end))
  }
  
  
  ### Preparation for sea ice calculation
  
  # Load the ice-free union layer
  domain <- rast(here("Data/Environmental_predictors/ice_free_upsamp_1km.tif"))
  
  # Set domain values
  domain <- ifel(not.na(domain), 1, NA)
  
  # domain <- aggregate(domain, fact = 100, fun = mean, na.rm = T, na.omit = T)
  
  # Turn cell centres into points
  domain.pts <- as.points(domain, values = T)
  domain.pts <- domain.pts[domain.pts$rock_union1 == 1, ]
  
  # For every cell centre, buffer by 100km
  domain.pts_buffer <- terra::buffer(domain.pts, 100000)
  
  # Load the coast to crop out land ice
  coast <- vect(here("Data/Environmental_predictors/add_coastline_high_res_polygon_v7_10.shp"))
  
  # Load the polygon delimiting South Georgia and South America (currently 100 % SIC)
  
  SG <- vect(here("Data", "PolarRes", "orkney.shp"))
  
  ### Extract the sea-ice concentration for every day, then average for November
  
  for(y in seq_along(years_hist)) {
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the first day of the next month
    
    # Get the rasters for the year, but exclude the file that spans from previous year (DEC)
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
    
    # Now select just November
    Doy <- get_doy_range(years_hist[y], 11)
    
    r <- r[[Doy]]
    
    # Crop out the coast
    r <- mask(r, coast, inverse = TRUE, touches = F)
    
    # Find South Georgia and South America cells
    IDs <- extract(r, SG, cells = TRUE)
    IDs <- IDs[!is.na(IDs[,2]) & IDs[,2] == 100, ]
    # Set those cells to zero in the original raster
    r[IDs$cell] <- 0
    
    # Get the mean November sea ice concentration
    r <- app(r, mean, na.rm = T)
    
    # Save the yearly mean November sea ice conc
    name <- paste0(outpath, "/Mean_November_Sea_Ice_Concentration_Year_", years_hist[y],".tif")
    
    writeRaster(r, name, gdal=c("COMPRESS=NONE"), overwrite = T)
    
    # Make all areas that are not sea ice 0 not NA
    # So that they're counted in mean sea ice
    r <- ifel(is.na(r), 0, r)
    
    # Pull out the mean sea ice concentration value within the buffer for every cell
    extracted_mean <- terra::extract(r, domain.pts_buffer, fun = mean, na.rm = T)
    
    # Classify NAs as 0
    extracted_mean[is.na(extracted_mean)] <- 0
    
    # Create a new empty raster matching domain
    sea_ice_NOV <- domain
    values(sea_ice_NOV) <- NA
    
    # Find the cell numbers of the domain pts
    cell_ids <- cellFromXY(sea_ice_NOV, crds(domain.pts))
    
    # Assign extracted sea ice values to these cells
    sea_ice_NOV[cell_ids] <- extracted_mean[, 2]
    
    # Save mean for that month per year
    name <- paste0(outpath, "/Mean_November_Sea_Ice_Concentration_100km_Year_", years_hist[y], ".tif")
    
    writeRaster(sea_ice_NOV, name, gdal=c("COMPRESS=NONE"), overwrite = T)
    
    tmp_files <- list.files(tmp_dir, full.names = T, pattern = "^file")
    
    file.remove(tmp_files)
    
    
  }
  
  
  # NOVEMBER MEAN SEA ICE FOR ALL YEARS -------------------------------------
  
  # 1. Concentration in buffer from ice-free land
  
  # List all the yearly mean rasters you just saved
  november_mean_files <- list.files(outpath, pattern = "^Mean_November_Sea_Ice_Concentration_100km_Year_.*\\.tif$", full.names = TRUE)
  
  # Subset just the relevant ones based on year:
  # Collapse the years into a regex pattern: "1995|1996|...|2014"
  year_pattern <- paste0("(", paste(years_hist, collapse = "|"), ")")
  
  # Subset the files
  november_mean_files <- november_mean_files[grepl(year_pattern, november_mean_files)]
  
  
  # Stack them together
  november_means <- rast(november_mean_files)
  
  # Calculate the overall mean
  final_mean <- app(november_means, mean, na.rm = TRUE)
  
  # Save the final mean raster
  writeRaster(final_mean, paste0(outpath, "/Mean_November_Sea_Ice_Concentration_100km_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  # 2. Now just concentration
  
  # List all the yearly mean rasters you just saved
  november_mean_files <- list.files(outpath, pattern = "^Mean_November_Sea_Ice_Concentration_Year_.*\\.tif$", full.names = TRUE)
  
  # Subset just the relevant ones based on year:
  # Collapse the years into a regex pattern: "1995|1996|...|2014"
  year_pattern <- paste0("(", paste(years_hist, collapse = "|"), ")")
  
  # Subset the files
  november_mean_files <- november_mean_files[grepl(year_pattern, november_mean_files)]
  
  
  # Stack them together
  november_means <- rast(november_mean_files)
  
  # Calculate the overall mean
  final_mean <- app(november_means, mean, na.rm = TRUE)
  
  # Save the final mean raster
  writeRaster(final_mean, paste0(outpath, "/Mean_November_Sea_Ice_Concentration_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # FUTURE --------------------------------------------------------------
  
  years_future <- seq(2081, 2100, by = 1)
  
  months = seq(1, 12, by=1)
  
  ### Extract the sea-ice concentration for every day, then average for November
  
  for(y in seq_along(years_future)) {
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the first day of the next month
    
    # Get the rasters for the year, but exclude the file that spans from previous year
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_future[y]),
                                       x = variable_paths, perl = TRUE)]
    year_files <- year_files[!grepl(pattern = paste0((years_future[y] - 1), "12"), x = year_files)]
    
    print(years_future[y])
    print(length(year_files))
    
    # Load each file, remove its last layer, then combine
    r_list <- lapply(year_files, function(file) {
      r_temp <- terra::rast(file)
      r_temp[[1:(terra::nlyr(r_temp) - 1)]]  # Remove last layer
    })
    
    # Combine all into one raster
    r <- terra::rast(r_list)
    r
    
    # Now select just November
    Doy <- get_doy_range(years_future[y], 11)
    
    r <- r[[Doy]]
    
    # Crop out the coast
    r <- mask(r, coast, inverse = TRUE, touches = F)
    
    # Find South Georgia and South America cells
    IDs <- extract(r, SG, cells = TRUE)
    IDs <- IDs[!is.na(IDs[,2]) & IDs[,2] == 100, ]
    # Set those cells to zero in the original raster
    r[IDs$cell] <- 0
    
    # Get the mean November sea ice concentration
    r <- app(r, mean, na.rm = T)
    
    # Save the yearly mean November sea ice conc
    name <- paste0(outpath, "/Mean_November_Sea_Ice_Concentration_Year_", years_future[y],".tif")
    
    writeRaster(r, name, gdal=c("COMPRESS=NONE"), overwrite = T)
    
    # Make all areas that are not sea ice 0 not NA
    # So that they're counted in mean sea ice
    r <- ifel(is.na(r), 0, r)
    
    # Pull out the mean sea ice concentration value within the buffer for every cell
    extracted_mean <- terra::extract(r, domain.pts_buffer, fun = mean, na.rm = T)
    
    # Classify NAs as 0
    extracted_mean[is.na(extracted_mean)] <- 0
    
    # Create a new empty raster matching domain
    sea_ice_NOV <- domain
    values(sea_ice_NOV) <- NA
    
    # Find the cell numbers of the domain pts
    cell_ids <- cellFromXY(sea_ice_NOV, crds(domain.pts))
    
    # Assign extracted sea ice values to these cells
    sea_ice_NOV[cell_ids] <- extracted_mean[, 2]
    
    # Save mean for that month per year
    name <- paste0(outpath, "/Mean_November_Sea_Ice_Concentration_100km_Year_", years_future[y], ".tif")
    
    writeRaster(sea_ice_NOV, name, gdal=c("COMPRESS=NONE"), overwrite = T)
    
    tmp_files <- list.files(tmp_dir, full.names = T, pattern = "^file")
    
    file.remove(tmp_files)
    
    
  }
  
  
  # NOVEMBER MEAN SEA ICE FOR ALL YEARS (FUTURE) -------------------------------------
  
  # 1. Concentration in buffer from ice-free land
  
  # List all the yearly mean rasters you just saved
  november_mean_files <- list.files(outpath, pattern = "^Mean_November_Sea_Ice_Concentration_100km_Year_.*\\.tif$", full.names = TRUE)
  
  # Subset just the relevant ones based on year:
  # Collapse the years into a regex pattern: "1995|1996|...|2014"
  year_pattern <- paste0("(", paste(years_future, collapse = "|"), ")")
  
  # Subset the files
  november_mean_files <- november_mean_files[grepl(year_pattern, november_mean_files)]
  
  
  # Stack them together
  november_means <- rast(november_mean_files)
  
  # Calculate the overall mean
  final_mean <- app(november_means, mean, na.rm = TRUE)
  
  # Save the final mean raster
  writeRaster(final_mean, paste0(outpath, "/Mean_November_Sea_Ice_Concentration_100km_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  # 2. Now just concentration
  
  # List all the yearly mean rasters you just saved
  november_mean_files <- list.files(outpath, pattern = "^Mean_November_Sea_Ice_Concentration_Year_.*\\.tif$", full.names = TRUE)
  
  # Subset just the relevant ones based on year:
  # Collapse the years into a regex pattern: "1995|1996|...|2014"
  year_pattern <- paste0("(", paste(years_future, collapse = "|"), ")")
  
  # Subset the files
  november_mean_files <- november_mean_files[grepl(year_pattern, november_mean_files)]
  
  
  # Stack them together
  november_means <- rast(november_mean_files)
  
  # Calculate the overall mean
  final_mean <- app(november_means, mean, na.rm = TRUE)
  
  # Save the final mean raster
  writeRaster(final_mean, paste0(outpath, "/Mean_November_Sea_Ice_Concentration_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
}



# 5. SEA ICE (OCTOBER TO FEB)  -------------------------------------------------------------
# These are stored in a MONTHLY format, with daily values per month
# Each monthly raster includes the first day of the next month, which I deal with

# **ALSO** there were two corrupt file which weren't able to be downloaded: siconca_sfx_ANT11_ANT11_ALADIN43_v1_CESM2_r11i1p1f1_historical_day_199404010000-199405010000.nc
# APRIL 1994 (30 days)

# siconca_sfx_ANT11_ANT11_ALADIN43_v1_CESM2_r11i1p1f1_historical_day_199406010000-199407010000.nc
# JUNE 1994 (30 days)

# Find sea ice
if(variable == "sea_ice_oct_feb"){
  
  variable_names <- files[grepl(pattern = "siconca", files)]
  variable_paths <- file_paths[grepl(pattern = "siconca", file_paths)]
  
  
  library(lubridate)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    return(seq(doy_start, doy_end))
  }
  
  
  
  ### Preparation for sea ice calculation
  
  # Load the ice-free union layer
  domain <- rast(here("Data/Environmental_predictors/ice_free_upsamp_1km.tif"))
  
  # Set domain values
  domain <- ifel(not.na(domain), 1, NA)
  
  # domain <- aggregate(domain, fact = 100, fun = mean, na.rm = T, na.omit = T)
  
  # Turn cell centres into points
  domain.pts <- as.points(domain, values = T)
  domain.pts <- domain.pts[domain.pts$rock_union1 == 1, ]
  
  # For every cell centre, buffer by 100km
  domain.pts_buffer32 <- terra::buffer(domain.pts, 32000)
  domain.pts_buffer36 <- terra::buffer(domain.pts, 36000)
  domain.pts_buffer42 <- terra::buffer(domain.pts, 42000)
  
  # Load the coast to crop out land ice
  coast <- vect(here("Data/Environmental_predictors/add_coastline_high_res_polygon_v7_10.shp"))
  
  # Load the polygon delimiting South Georgia and South America (currently 100 % SIC)
  
  SG <- vect(here("Data", "PolarRes", "orkney.shp"))
  
  
  # HISTORICAL --------------------------------------------------------------
  
  
  ### Extract the sea-ice concentration for every month
  
  months = seq(1, 12, by=1)
  years_hist <- seq(1994, 2014, by = 1) # DO 1994 to not have to remove
  
  # Store all monthly means organized by month across years
  
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name # month.name is a built in constant
  
  
  for(y in seq_along(years_hist)) {
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the first day of the next month
    
    # Get the rasters for the year, but exclude the file that spans from previous year
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
    
    # ---- INSERT NA filler for missing rsus April/June 1994 (31 days) ----
    # This creates an NA raster for the missing month to preserve the right number of days
    if(years_hist[y] == 1994) {
      april_doys <- get_doy_range(1994, 4)          # days 91-120 (30 days)
      june_doys <- get_doy_range(1994, 6)           # days 152-181 (30 days)
      filler_april <- rast(replicate(length(april_doys),  # 30 NA layers matching grid
                               init(r[[1]], fun = NA)))
      
      filler_june <- rast(replicate(length(june_doys),  # 30 NA layers matching grid
                                     init(r[[1]], fun = NA)))
      
      # Split r at the April insertion point, then reassemble
      r <- rast(list(
        r[[1:(april_doys[1] - 1)]],   # days before April
        filler_april,                            # NA placeholder for April
        r[[april_doys[1]:nlyr(r)]] # days after April
      ))
      
      # Split r at the June insertion point, then reassemble
      june_split <- june_doys[1] + length(april_doys)
      r <- rast(list(
        r[[1:(june_split - 1)]],
        filler_june,
        r[[june_split:nlyr(r)]]
      ))
    }
    # ---------------------------------------------------------------
    
    for(m in seq_along(months)) {
      
      # ---- SKIP April & June 1994 entirely (data missing for that month) ----
      if(years_hist[y] == 1994 && m == 4) {
        message("Skipping April 1994 - corrupt file")
        monthly_means_by_month[[4]][[y]] <- NA  # NA preserves the slot
        next
      }
      
      if(years_hist[y] == 1994 && m == 6) {
        message("Skipping April 1994 - corrupt file")
        monthly_means_by_month[[6]][[y]] <- NA  # NA preserves the slot
        next
      }
      # -------------------------------------------------------------------
      
      Doy <- get_doy_range(years_hist[y],months[m])
      
      # Subset the relevant layers for the month
      r_month <- r[[Doy]]
      
      # Crop out the coast
      r_month <- mask(r_month, coast, inverse = TRUE, touches = F)
      
      # Find South Georgia and South America cells
      IDs <- extract(r_month, SG, cells = TRUE)
      IDs <- IDs[!is.na(IDs[,2]) & IDs[,2] == 100, ]
      # Set those cells to zero in the original raster
      r_month[IDs$cell] <- 0
      
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
  
  # October 1994 - 2013 (for summers 1995 - 2014)
  # November 1994 - 2013 (for summers 1995 - 2014)
  # December 1994 - 2013 (for summers 1995 - 2014)
  # January 1995 - 2014 (for summers 1995 - 2014)
  # February 1995 - 2014 (for summers 1995 - 2014)
  
  
  # CONCENTRATION FIRST -----------------------------------------------------
  
  # OCTOBER
  
  m = 10
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove final layer (October *2014*)
  month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  
  # Calculate mean across all years
  climatological_mean_OCT <- app(month_stack, mean, na.rm = TRUE)
  
  
  # NOVEMBER
  
  m = 11
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove final layer (November *2014*)
  month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  
  # Calculate mean across all years
  climatological_mean_NOV <- app(month_stack, mean, na.rm = TRUE)
  
  
  # DECEMBER
  
  m = 12
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove final layer (December *2014*)
  month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  
  # Calculate mean across all years
  climatological_mean_DEC <- app(month_stack, mean, na.rm = TRUE)
  
  ## JANUARY
  
  m = 1
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove first layer (January 1994)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_JAN <- app(month_stack, mean, na.rm = TRUE)
  
  ## FEBRUARY
  
  m = 2
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove first layer (February 1994)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_FEB <- app(month_stack, mean, na.rm = TRUE)
  
  
  # STEP 3: Calculate mean sea ice from the 5 climatological monthly means
  # Summer = ONDJF (October, November, December, January, February)
  
  summer_stack <- c(climatological_mean_OCT, climatological_mean_NOV, climatological_mean_DEC, climatological_mean_JAN, climatological_mean_FEB)
  mean_summer_sea_ice <- app(summer_stack, mean, na.rm = TRUE)
  
  
  # Save the climatological mean summer sea ice concentration
  writeRaster(mean_summer_sea_ice,
              paste0(outpath, "/Mean_Oct_Feb_Sea_Ice_Concentration_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # SEA ICE CONCENTRATION IN BUFFER NOW  ------------------------------------
  
  # Make all areas that are not sea ice 0, not NA (COASTLINE)
  # So that they're counted in mean sea ice in buffer
  mean_summer_sea_ice <- ifel(is.na(mean_summer_sea_ice), 0, mean_summer_sea_ice)
  
  # Pull out the mean sea ice concentration value within the buffer for every cell
  extracted_mean32 <- terra::extract(mean_summer_sea_ice, domain.pts_buffer32, fun = mean, na.rm = T)
  extracted_mean36 <- terra::extract(mean_summer_sea_ice, domain.pts_buffer36, fun = mean, na.rm = T)
  extracted_mean42 <- terra::extract(mean_summer_sea_ice, domain.pts_buffer42, fun = mean, na.rm = T)
  
  # Classify NAs as 0
  extracted_mean32[is.na(extracted_mean32)] <- 0
  extracted_mean36[is.na(extracted_mean36)] <- 0
  extracted_mean42[is.na(extracted_mean42)] <- 0
  
  # Create a new empty raster matching domain
  sea_ice32 <- domain
  values(sea_ice32) <- NA
  sea_ice36 <- domain
  values(sea_ice36) <- NA
  sea_ice42 <- domain
  values(sea_ice42) <- NA
  
  
  # Find the cell numbers of the domain pts
  cell_ids32 <- cellFromXY(sea_ice32, crds(domain.pts))
  cell_ids36 <- cellFromXY(sea_ice36, crds(domain.pts))
  cell_ids42 <- cellFromXY(sea_ice42, crds(domain.pts))
  
  # Assign extracted sea ice values to these cells
  sea_ice32[cell_ids32] <- extracted_mean32[, 2]
  sea_ice36[cell_ids36] <- extracted_mean36[, 2]
  sea_ice42[cell_ids42] <- extracted_mean42[, 2]
  
  # Save mean for that buffer
  writeRaster(sea_ice32,
              paste0(outpath, "/Mean_Oct_Feb_Sea_Ice_Concentration_32km_HISTORICAL_1995_2014.tif"), gdal=c("COMPRESS=NONE"), overwrite = T)
  writeRaster(sea_ice36,
              paste0(outpath, "/Mean_Oct_Feb_Sea_Ice_Concentration_36km_HISTORICAL_1995_2014.tif"), gdal=c("COMPRESS=NONE"), overwrite = T)
  writeRaster(sea_ice42,
              paste0(outpath, "/Mean_Oct_Feb_Sea_Ice_Concentration_42km_HISTORICAL_1995_2014.tif"), gdal=c("COMPRESS=NONE"), overwrite = T)
  
  
  
  
  
  # FUTURE --------------------------------------------------------------
  
  
  ### Extract the sea-ice concentration for every month
  
  months = seq(1, 12, by=1)
  years_future <- seq(2080, 2100, by = 1) # DO 2080 to not have to remove
  
  
  # Store all monthly means organized by month across years
  
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name # month.name is a built in constant
  
  
  for(y in seq_along(years_future)) {
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the first day of the next month
    
    # Get the rasters for the year, but exclude the file that spans from previous year
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_future[y]),
                                       x = variable_paths, perl = TRUE)]
    year_files <- year_files[!grepl(pattern = paste0((years_future[y] - 1), "12"), x = year_files)]
    
    print(years_future[y])
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
      
      Doy <- get_doy_range(years_future[y],months[m])
      
      # Subset the relevant layers for the month
      r_month <- r[[Doy]]
      
      # Crop out the coast
      r_month <- mask(r_month, coast, inverse = TRUE, touches = F)
      
      # Find South Georgia and South America cells
      IDs <- extract(r_month, SG, cells = TRUE)
      IDs <- IDs[!is.na(IDs[,2]) & IDs[,2] == 100, ]
      # Set those cells to zero in the original raster
      r_month[IDs$cell] <- 0
      
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
  
  # October 2080 - 2099 (for summers 2081 - 2100)
  # November 2080 - 2099 (for summers 2081 - 2100)
  # December 2080 - 2099 (for summers 2081 - 2100)
  # January 2081 - 2100 (for summers 2081 - 2100)
  # February 2081 - 2100 (for summers 2081 - 2100)
  
  
  # CONCENTRATION FIRST -----------------------------------------------------
  
  # OCTOBER
  
  m = 10
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove final layer (October *2100*)
  month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  
  # Calculate mean across all years
  climatological_mean_OCT <- app(month_stack, mean, na.rm = TRUE)
  
  
  # NOVEMBER
  
  m = 11
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove final layer (November *2100*)
  month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  
  # Calculate mean across all years
  climatological_mean_NOV <- app(month_stack, mean, na.rm = TRUE)
  
  
  # DECEMBER
  
  m = 12
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove final layer (December *2100*)
  month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  
  # Calculate mean across all years
  climatological_mean_DEC <- app(month_stack, mean, na.rm = TRUE)
  
  ## JANUARY
  
  m = 1
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove first layer (January 2080)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_JAN <- app(month_stack, mean, na.rm = TRUE)
  
  ## FEBRUARY
  
  m = 2
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove first layer (February 2080)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_FEB <- app(month_stack, mean, na.rm = TRUE)
  
  
  # STEP 3: Calculate mean sea ice from the 5 climatological monthly means
  # Summer = ONDJF (October, November, December, January, February)
  
  summer_stack <- c(climatological_mean_OCT, climatological_mean_NOV, climatological_mean_DEC, climatological_mean_JAN, climatological_mean_FEB)
  mean_summer_sea_ice <- app(summer_stack, mean, na.rm = TRUE)
  
  
  # Save the climatological mean summer sea ice concentration
  writeRaster(mean_summer_sea_ice,
              paste0(outpath, "/Mean_Oct_Feb_Sea_Ice_Concentration_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # SEA ICE CONCENTRATION IN BUFFER NOW  ------------------------------------
  
  # Make all areas that are not sea ice 0, not NA (COASTLINE)
  # So that they're counted in mean sea ice in buffer
  mean_summer_sea_ice <- ifel(is.na(mean_summer_sea_ice), 0, mean_summer_sea_ice)
  
  # Pull out the mean sea ice concentration value within the buffer for every cell
  extracted_mean32 <- terra::extract(mean_summer_sea_ice, domain.pts_buffer32, fun = mean, na.rm = T)
  extracted_mean36 <- terra::extract(mean_summer_sea_ice, domain.pts_buffer36, fun = mean, na.rm = T)
  extracted_mean42 <- terra::extract(mean_summer_sea_ice, domain.pts_buffer42, fun = mean, na.rm = T)
  
  # Classify NAs as 0
  extracted_mean32[is.na(extracted_mean32)] <- 0
  extracted_mean36[is.na(extracted_mean36)] <- 0
  extracted_mean42[is.na(extracted_mean42)] <- 0
  
  # Create a new empty raster matching domain
  sea_ice32 <- domain
  values(sea_ice32) <- NA
  sea_ice36 <- domain
  values(sea_ice36) <- NA
  sea_ice42 <- domain
  values(sea_ice42) <- NA
  
  
  # Find the cell numbers of the domain pts
  cell_ids32 <- cellFromXY(sea_ice32, crds(domain.pts))
  cell_ids36 <- cellFromXY(sea_ice36, crds(domain.pts))
  cell_ids42 <- cellFromXY(sea_ice42, crds(domain.pts))
  
  # Assign extracted sea ice values to these cells
  sea_ice32[cell_ids32] <- extracted_mean32[, 2]
  sea_ice36[cell_ids36] <- extracted_mean36[, 2]
  sea_ice42[cell_ids42] <- extracted_mean42[, 2]
  
  # Save mean for that buffer
  writeRaster(sea_ice32,
              paste0(outpath, "/Mean_Oct_Feb_Sea_Ice_Concentration_32km_FUTURE_2081_2100.tif"), gdal=c("COMPRESS=NONE"), overwrite = T)
  writeRaster(sea_ice36,
              paste0(outpath, "/Mean_Oct_Feb_Sea_Ice_Concentration_36km_FUTURE_2081_2100.tif"), gdal=c("COMPRESS=NONE"), overwrite = T)
  writeRaster(sea_ice42,
              paste0(outpath, "/Mean_Oct_Feb_Sea_Ice_Concentration_42km_FUTURE_2081_2100.tif"), gdal=c("COMPRESS=NONE"), overwrite = T)
  
  
}




# 6. TOTAL ANNUAL PRECIPITATION----------------------------------------------------------
# (kg m-2 s-1 (precipitation flux))
# Daily

if(variable == "total_precip"){
  
  # Find TT all years
  variable_names <- files[grepl(pattern = "pr", files)]
  variable_paths <- file_paths[grepl(pattern = "pr", file_paths)]
  
  
  # HISTORICAL --------------------------------------------------------------
  
  # For every year, calculate the total annual precipitation
  
  years_hist <- seq(1995, 2014, by = 1)
  
  # For every year
  for(y in seq_along(years_hist)) {
    
    r <- terra::rast(variable_paths[grepl(variable_paths, pattern = years_hist[y])])
    
    # Convert daily flux to daily mm by multiplying by seconds
    r <- r*86400
    
    # Sum these to get the yearly precip sum:
    yearly_precip_sum <- app(r, sum, na.rm = TRUE)
    
    # Save the yearly sum
    name <- paste0(outpath, "/Total_Annual_Precipitation_Year_", years_hist[y],".tif")
    
    writeRaster(yearly_precip_sum, name, gdal=c("COMPRESS=NONE"), overwrite = T)
    
    tmp_files <- list.files(tmp_dir, full.names = T, pattern = "^file")
    
    file.remove(tmp_files)
    
  }
  
  
  # MEAN of TOTAL ANNUAL PRECIP FOR ALL HISTORICAL YEARS -----------------------
  
  # List all the yearly mean rasters you just saved
  annual_mean_files <- list.files(outpath, pattern = "^Total_Annual_Precipitation_Year_.*\\.tif$", full.names = TRUE)
  
  # Subset just the relevant ones based on year:
  # Collapse the years into a regex pattern: "1995|1996|...|2014"
  year_pattern <- paste0("(", paste(years_hist, collapse = "|"), ")")
  
  # Subset the files
  annual_mean_files <- annual_mean_files[grepl(year_pattern, annual_mean_files)]
  
  
  # Stack them together
  annual_means <- rast(annual_mean_files)
  
  # Calculate the overall mean
  final_mean <- app(annual_means, mean, na.rm = TRUE)
  
  # Save the final mean raster
  writeRaster(final_mean, paste0(outpath, "/Total_Annual_Precipitation_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  
  # FUTURE YEARS ------------------------------------------------------------
  
  # For every year, calculate the mean
  
  years_future <- seq(2081, 2100, by = 1)
  
  # For every year
  for(y in seq_along(years_future)) {
    
    r <- terra::rast(variable_paths[grepl(variable_paths, pattern = years_future[y])])
    
    # Convert daily flux to daily mm by multiplying by seconds
    r <- r*86400
    
    # Sum these to get the yearly precip sum:
    yearly_precip_sum <- app(r, sum, na.rm = TRUE)
    
    # Save the yearly sum
    name <- paste0(outpath, "/Total_Annual_Precipitation_Year_", years_future[y],".tif")
    
    writeRaster(yearly_precip_sum, name, gdal=c("COMPRESS=NONE"), overwrite = T)
    
    tmp_files <- list.files(tmp_dir, full.names = T, pattern = "^file")
    
    file.remove(tmp_files)
    
  }
  
  
  # ANNUAL MEAN FOR ALL FUTURE YEARS -----------------------------------------------
  
  # List all the yearly mean rasters you just saved
  annual_mean_files <- list.files(outpath, pattern = "^Total_Annual_Precipitation_Year_.*\\.tif$", full.names = TRUE)
  
  # Subset just the relevant ones based on year:
  # Collapse the years into a regex pattern: "2081|2082|...|2100"
  year_pattern <- paste0("(", paste(years_future, collapse = "|"), ")")
  
  # Subset the files
  annual_mean_files <- annual_mean_files[grepl(year_pattern, annual_mean_files)]
  
  # Stack them together
  annual_means <- rast(annual_mean_files)
  
  # Calculate the overall mean
  final_mean <- app(annual_means, mean, na.rm = TRUE)
  
  # Save the final mean raster
  writeRaster(final_mean, paste0(outpath, "/Total_Annual_Precipitation_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
}


# 6b. TOTAL SUMMER PRECIPITATION----------------------------------------------------------
# (kg m-2 s-1 (precipitation flux))
# Daily

if(variable == "total_summer_precip"){
  
  # Find TT all years
  variable_names <- files[grepl(pattern = "pr", files)]
  variable_paths <- file_paths[grepl(pattern = "pr", file_paths)]
  
  
  # MONTHLY VALUES HISTORICAL  -------------------------------------------------------------
  
  months = seq(1, 12, by=1)
  years_hist <- seq(1994, 2014, by = 1) # DO 1994 to not have to remove
  
  library(lubridate)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    return(seq(doy_start, doy_end))
  }
  
  
  annual_rasters <- list()
  monthly_rasters <- list()
  
  for(y in seq_along(years_hist)) {
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_hist[y],months[m])
      
      # Get the rasters for the year then subset for the month of interest
      r <- terra::rast(variable_paths[grepl(variable_paths, pattern = years_hist[y])])
      
      # Subset the relevant layers for the month
      r <- r[[Doy]]
      
      # Convert daily flux to daily mm by multiplying by seconds
      r <- r*86400
      
      monthly_precip_sum <- app(r, sum, na.rm = TRUE)
      
      # Store in list
      monthly_rasters[[paste0("Year_", years_hist[y], "_Month_", months[m])]] <- monthly_precip_sum
      
    }
    
  }
  
  
  # Extract the summer months
  
  summer_total_rasters <- list()
  
  # UPDATE YEARS TO REMOVE 1994
  years_hist <- seq(1995, 2014, by = 1)
  
  
  for(y in years_hist) {
    
    dec_prev_year <- paste0("Year_", y - 1, "_Month_12")
    jan_curr_year <- paste0("Year_", y, "_Month_1")
    feb_curr_year <- paste0("Year_", y, "_Month_2")
    
    # Check which of these actually exist (to avoid errors if edge years are missing)
    existing_months <- c(dec_prev_year, jan_curr_year, feb_curr_year)
    existing_months <- existing_months[existing_months %in% names(monthly_rasters)]
    
    summer_stack <- rast(monthly_rasters[existing_months])
    
    # Take total for the summer
    summer_total <- app(summer_stack, sum, na.rm = TRUE)
    
    # Store
    summer_total_rasters[[paste0("DJF_", y)]] <- summer_total
    
  }
  
  
  # Calculate the mean for summers over all years ---------------------------
  
  # Combine them into one raster
  summer_total_rasters <- rast(summer_total_rasters)
  
  summer_total_ALL <- app(summer_total_rasters, mean, na.rm = TRUE)
  
  # Save the final mean raster
  writeRaster(summer_total_ALL, paste0(outpath, "/Mean_Summer_Total_Precipitation_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  
  # MONTHLY VALUES FUTURE  -----------------------------------------------------
  
  months = seq(1, 12, by=1)
  
  years_future <- seq(2080, 2100, by = 1) # DO 2080 to not have to remove
  
  library(lubridate)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    return(seq(doy_start, doy_end))
  }
  
  
  annual_rasters <- list()
  monthly_rasters <- list()
  
  for(y in seq_along(years_future)) {
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_future[y],months[m])
      
      # Get the rasters for the year then subset for the month of interest
      r <- terra::rast(variable_paths[grepl(variable_paths, pattern = years_future[y])])
      
      # Subset the relevant layers for the month
      r <- r[[Doy]]
      
      # Convert daily flux to daily mm by multiplying by seconds
      r <- r*86400
      
      monthly_precip_sum <- app(r, sum, na.rm = TRUE)
      
      # Store in list
      monthly_rasters[[paste0("Year_", years_future[y], "_Month_", months[m])]] <- monthly_precip_sum
      
    }
    
  }
  
  
  # Extract the summer months
  
  summer_total_rasters <- list()
  
  # UPDATE YEARS TO REMOVE 2080
  years_future <- seq(2081, 2100, by = 1)
  
  
  for(y in years_future) {
    
    dec_prev_year <- paste0("Year_", y - 1, "_Month_12")
    jan_curr_year <- paste0("Year_", y, "_Month_1")
    feb_curr_year <- paste0("Year_", y, "_Month_2")
    
    # Check which of these actually exist (to avoid errors if edge years are missing)
    existing_months <- c(dec_prev_year, jan_curr_year, feb_curr_year)
    existing_months <- existing_months[existing_months %in% names(monthly_rasters)]
    
    summer_stack <- rast(monthly_rasters[existing_months])
    
    # Take total for the summer
    summer_total <- app(summer_stack, sum, na.rm = TRUE)
    
    # Store
    summer_total_rasters[[paste0("DJF_", y)]] <- summer_total
    
  }
  
  
  # Calculate the mean for summers over all years ---------------------------
  
  
  # Combine them into one raster
  summer_total_rasters <- rast(summer_total_rasters)
  
  summer_total_ALL <- app(summer_total_rasters, mean, na.rm = TRUE)
  
  # Save the final mean raster
  writeRaster(summer_total_ALL, paste0(outpath, "/Mean_Summer_Total_Precipitation_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
}




# 6c. MEAN ANNUAL PRECIPITATION -----------------------------------
# (kg m-2 s-1 (precipitation flux))

if(variable == "mean_precip"){
  
  # Find tas all years
  variable_names <- files[grepl(pattern = "pr", files)]
  variable_paths <- file_paths[grepl(pattern = "pr", file_paths)]
  
  
  # HISTORICAL --------------------------------------------------------------
  
  library(lubridate)
  
  years_hist <- seq(1995, 2014, by = 1)
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
    r <- terra::rast(variable_paths[grepl(variable_paths, pattern = years_hist[y])])
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_hist[y], months[m])
      
      # Subset for the month of interest (1 raster per day)
      r_month <- r[[Doy]]
      
      # Convert from Kelvin to Celsius
      r_month <- r_month*86400
      
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
    name <- paste0(outpath, "/Climatological_Monthly_Mean_Precipitation_",
                   month.name[m], "_1995_2014.tif")
    writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
    
    # Store in list for annual calculation
    climatological_monthly_means[[m]] <- climatological_mean
    
    # Clean up temp files
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 3: Calculate mean annual temperature from the 12 climatological monthly means
  annual_stack <- rast(climatological_monthly_means)
  mean_annual_precip <- app(annual_stack, mean, na.rm = TRUE)
  
  # Save the final climatological mean annual temperature
  writeRaster(mean_annual_precip,
              paste0(outpath, "/Mean_Annual_Precipitation_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # FUTURE CLIMATOLOGY (2081-2100) -----------------------------------------
  
  years_future <- seq(2081, 2100, by = 1)
  
  # STEP 1: Calculate monthly means for each year
  # Store all monthly means organized by month across years
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name
  
  for(y in seq_along(years_future)) {
    
    # Get the raster for the year
    r <- terra::rast(variable_paths[grepl(variable_paths, pattern = years_future[y])])
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_future[y], months[m])
      
      # Subset for the month of interest (1 raster per day)
      r_month <- r[[Doy]]
      
      # Convert from Kelvin to Celsius
      r_month <- r_month*86400
      
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
  climatological_monthly_means_future <- list()
  
  for(m in seq_along(months)) {
    
    # Stack all years for this month
    month_stack <- rast(monthly_means_by_month[[m]])
    
    # Calculate mean across all years
    climatological_mean <- app(month_stack, mean, na.rm = TRUE)
    
    # Save the climatological monthly mean
    month_name <- sprintf("%02d", m)
    name <- paste0(outpath, "/Climatological_Monthly_Mean_Precipitation_",
                   month.name[m], "_2081_2100.tif")
    writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
    
    # Store in list for annual calculation
    climatological_monthly_means_future[[m]] <- climatological_mean
    
    # Clean up temp files
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 3: Calculate mean annual temperature from the 12 climatological monthly means
  annual_stack_future <- rast(climatological_monthly_means_future)
  mean_annual_precip_future <- app(annual_stack_future, mean, na.rm = TRUE)
  
  # Save the final climatological mean annual temperature
  writeRaster(mean_annual_precip_future,
              paste0(outpath, "/Mean_Annual_Precipitation_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
}


# 6c. MEAN SUMMER PRECIPITATION -----------------------------------
# (kg m-2 s-1 (precipitation flux))

if(variable == "mean_summer_precip"){
  
  # Find PR all years
  variable_names <- files[grepl(pattern = "pr", files)]
  variable_paths <- file_paths[grepl(pattern = "pr", file_paths)]
  
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
  
  
  # HISTORICAL SUMMER PRECIP (1995-2014) ------------------------------
  # Have to re-calculate to include December 1994
  
  years_hist <- seq(1994, 2014, by = 1)  # Include 1994 to get Dec for 1995
  
  # STEP 1: Calculate monthly means for each year from daily values
  # Store all monthly means organized by month across years
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name
  
  for(y in seq_along(years_hist)) {
    
    # Get the raster for the year
    r <- terra::rast(variable_paths[grepl(variable_paths, pattern = years_hist[y])])
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_hist[y], months[m])
      
      # Subset for the month of interest (1 raster per day)
      r_month <- r[[Doy]]
      
      # Convert daily flux to daily mm by multiplying by seconds
      r_month <- r_month*86400
      
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
  
  # December 1994 - 2013 (for summers 1995 - 2014)
  # January 1995 - 2014 (for summers 1995 - 2014)
  # February 1995 - 2014 (for summers 1995 - 2014)
  
  ## JANUARY
  
  m = 1
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove first layer (January 1994)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_JAN <- app(month_stack, mean, na.rm = TRUE)
  
  ## FEBRUARY
  
  m = 2
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove first layer (February 1994)
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
  
  
  # STEP 3: Calculate mean summer snow cover from the 3 climatological monthly means
  # Summer = DJF (December, January, February)
  
  summer_stack <- c(climatological_mean_DEC, climatological_mean_JAN, climatological_mean_FEB)
  mean_summer_temp <- app(summer_stack, mean, na.rm = TRUE)
  
  
  
  # Save the climatological mean summer temperature
  writeRaster(mean_summer_temp,
              paste0(outpath, "/Mean_Summer_Precipitation_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # FUTURE SUMMER PRECIP (2081-2100) ----------------------------------
  # Have to re-calculate to include December 2080
  
  years_future <- seq(2080, 2100, by = 1)  # Include 2080 to get Dec for 2081
  
  # STEP 1: Calculate monthly means for each year from daily values
  # Store all monthly means organized by month across years
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name
  
  for(y in seq_along(years_future)) {
    
    # Get the raster for the year
    r <- terra::rast(variable_paths[grepl(variable_paths, pattern = years_future[y])])
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_future[y], months[m])
      
      # Subset for the month of interest (1 raster per day)
      r_month <- r[[Doy]]
      
      # Convert daily flux to daily mm by multiplying by seconds
      r_month <- r_month*86400
      
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
  
  # December 2080 - 2099 (for summers 2081 - 2100)
  # January 2081 - 2100 (for summers 2081 - 2100)
  # February 2081 - 2100 (for summers 2081 - 2100)
  
  
  ## JANUARY
  
  m = 1
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # # Remove final layer (January 2100)
  # month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  # Remove first layer (January 2080)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_JAN <- app(month_stack, mean, na.rm = TRUE)
  
  ## FEBRUARY
  
  m = 2
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # # Remove final layer (February 2100)
  # month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  # Remove first layer (February 2080)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_FEB <- app(month_stack, mean, na.rm = TRUE)
  
  ## DECEMBER
  
  m = 12
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove final layer (December *2100*)
  month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  # # Remove other final layer (December *2013*)
  # month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  
  # Calculate mean across all years
  climatological_mean_DEC <- app(month_stack, mean, na.rm = TRUE)
  
  
  # STEP 3: Calculate mean summer snow cover from the 3 climatological monthly means
  # Summer = DJF (December, January, February)
  
  summer_stack <- c(climatological_mean_DEC, climatological_mean_JAN, climatological_mean_FEB)
  mean_summer_temp_future <- app(summer_stack, mean, na.rm = TRUE)
  
  
  # Save the climatological mean summer temperature
  writeRaster(mean_summer_temp_future,
              paste0(outpath, "/Mean_Summer_Precipitation_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
}




# 8. SOLAR RADIATION (Watts per m2?) ---------------------------------------

# These are stored in a MONTHLY format, with daily values per month
# Each monthly raster includes the first day of the next month, which I deal with

# **ALSO** there was one corrupt file which wasn's able to be downloaded: rsus_fp_ANT11_ANT11_ALADIN43_v1_CESM2_r11i1p1f1_historical_day_199705010000-199706010000.nc
# MAY 1997 (31 days)

if(variable == "solar_rad"){
  
  # Find rsus all years
  variable_names_rsus <- files[grepl(pattern = "rsus", files)]
  variable_paths_rsus <- file_paths[grepl(pattern = "rsus", file_paths)]
  
  variable_names_rsds <- files[grepl(pattern = "rsds", files)]
  variable_paths_rsds <- file_paths[grepl(pattern = "rsds", file_paths)]
  
  
  # HISTORICAL --------------------------------------------------------------
  
  # For every year, pull out the monthly rasters, calculate their mean, then summarise for entire year
  
  years_hist <- seq(1995, 2014, by = 1)
  months = seq(1, 12, by=1)
  
  # Get the the days range for each month (what day index is in that month)
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
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month
    # But unlike SICONA - it doesn't appear to actually hold that day in it (so we don't have to remove any layers)
    
    # Get the rasters for the year, but exclude the file that spans from previous year (because of naming)
    year_files_rsus <- variable_paths_rsus[grepl(pattern = paste0("(?<=[_-])", years_hist[y]),
                                                 x = variable_paths_rsus, perl = TRUE)]
    year_files_rsus <- year_files_rsus[!grepl(pattern = paste0((years_hist[y] - 1), "12"), x = year_files_rsus)]
    
    year_files_rsds <- variable_paths_rsds[grepl(pattern = paste0("(?<=[_-])", years_hist[y]),
                                                 x = variable_paths_rsds, perl = TRUE)]
    year_files_rsds <- year_files_rsds[!grepl(pattern = paste0((years_hist[y] - 1), "12"), x = year_files_rsds)]
    
    print(years_hist[y])
    print(length(year_files_rsds))
    print(length(year_files_rsus))
    
    # Load each file then combine
    r_list_rsus <- lapply(year_files_rsus, function(file) {
      r_temp <- terra::rast(file)
    })
    
    r_list_rsds <- lapply(year_files_rsds, function(file) {
      r_temp <- terra::rast(file)
    })
    
    # Combine all into one raster
    r_rsus <- terra::rast(r_list_rsus)
    r_rsds <- terra::rast(r_list_rsds)
    
    r_rsus
    r_rsds
    
    # ---- INSERT NA filler for missing rsus May 1997 (31 days) ----
    # This creates an NA raster for the missing month to preserve the right number of days
    if(years_hist[y] == 1997) {
      may_doys <- get_doy_range(1997, 5)          # days 121-151 (31 days)
      filler <- rast(replicate(length(may_doys),  # 31 NA layers matching grid
                               init(r_rsds[[1]], fun = NA)))
      
      # Split r_rsus at the May insertion point, then reassemble
      r_rsus <- rast(list(
        r_rsus[[1:(may_doys[1] - 1)]],   # days before May
        filler,                            # NA placeholder for May
        r_rsus[[may_doys[1]:nlyr(r_rsus)]] # days after May
      ))
    }
    # ---------------------------------------------------------------
    
    # Checking that rsds and rsus now match in layer numbers
    if(nlyr(r_rsus) != nlyr(r_rsds)) stop(paste("Layer mismatch in year", years_hist[y]))
    
    # Daily net solar radiation over entire year
    SWnet <- r_rsds - r_rsus
    
    for(m in seq_along(months)) {
      
      # ---- SKIP May 1997 entirely (rsus data missing for that month) ----
      if(years_hist[y] == 1997 && m == 5) {
        message("Skipping May 1997 - corrupt rsus file")
        monthly_means_by_month[[5]][[y]] <- NA  # NA preserves the slot
        next
      }
      # -------------------------------------------------------------------
      
      
      Doy <- get_doy_range(years_hist[y],months[m])
      
      # Subset the relevant layers for the month
      r_month <- SWnet[[Doy]]
      
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
    
    # Stack all years for this month - *filter out NAs before stacking (due to May 1997 thing above)
    month_stack <- rast(Filter(function(x) inherits(x, "SpatRaster"), monthly_means_by_month[[m]]))
    
    # Calculate mean across all years
    climatological_mean <- app(month_stack, mean, na.rm = TRUE)
    
    # Save the climatological monthly mean
    month_name <- sprintf("%02d", m)
    name <- paste0(outpath, "/Climatological_Monthly_Mean_Solar_Radiation_",
                   month.name[m], "_1995_2014.tif")
    writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
    
    # Store in list for annual calculation
    climatological_monthly_means[[m]] <- climatological_mean
    
    # Clean up temp files
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 3: Calculate mean annual solar radiation from the 12 climatological monthly means
  annual_stack <- rast(climatological_monthly_means)
  mean_annual_temp <- app(annual_stack, mean, na.rm = TRUE)
  
  # Save the final climatological mean annual temperature
  writeRaster(mean_annual_temp,
              paste0(outpath, "/Mean_Annual_Solar_Radiation_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # FUTURE CLIMATOLOGY (2081 - 2100) -------------------------------------------
  
  # For every year, pull out the monthly rasters, calculate their mean, then summarise for entire year
  
  years_future <- seq(2081, 2100, by = 1)
  months = seq(1, 12, by=1)
  
  # Get the the days range for each month (what day index is in that month)
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
  
  
  for(y in seq_along(years_future)) {
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month
    # But unlike SICONA - it doesn't appear to actually hold that day in it (so we don't have to remove any layers)
    
    # Get the rasters for the year, but exclude the file that spans from previous year (because of naming)
    year_files_rsus <- variable_paths_rsus[grepl(pattern = paste0("(?<=[_-])", years_future[y]),
                                                 x = variable_paths_rsus, perl = TRUE)]
    year_files_rsus <- year_files_rsus[!grepl(pattern = paste0((years_future[y] - 1), "12"), x = year_files_rsus)]
    
    year_files_rsds <- variable_paths_rsds[grepl(pattern = paste0("(?<=[_-])", years_future[y]),
                                                 x = variable_paths_rsds, perl = TRUE)]
    year_files_rsds <- year_files_rsds[!grepl(pattern = paste0((years_future[y] - 1), "12"), x = year_files_rsds)]
    
    print(years_future[y])
    print(length(year_files_rsds))
    print(length(year_files_rsus))
    
    # Load each file then combine
    r_list_rsus <- lapply(year_files_rsus, function(file) {
      r_temp <- terra::rast(file)
    })
    
    r_list_rsds <- lapply(year_files_rsds, function(file) {
      r_temp <- terra::rast(file)
    })
    
    # Combine all into one raster
    r_rsus <- terra::rast(r_list_rsus)
    r_rsds <- terra::rast(r_list_rsds)
    
    r_rsus
    r_rsds
    
    # Daily net solar radiation over entire year
    SWnet <- r_rsds - r_rsus
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_future[y],months[m])
      
      # Subset the relevant layers for the month
      r_month <- SWnet[[Doy]]
      
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
  climatological_monthly_means_future <- list()
  
  for(m in seq_along(months)) {
    
    # Stack all years for this month
    month_stack <- rast(monthly_means_by_month[[m]])
    
    # Calculate mean across all years
    climatological_mean <- app(month_stack, mean, na.rm = TRUE)
    
    # Save the climatological monthly mean
    month_name <- sprintf("%02d", m)
    name <- paste0(outpath, "/Climatological_Monthly_Mean_Solar_Radiation_",
                   month.name[m], "_2081_2100.tif")
    writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
    
    # Store in list for annual calculation
    climatological_monthly_means_future[[m]] <- climatological_mean
    
    # Clean up temp files
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 3: Calculate mean annual solar radiation from the 12 climatological monthly means
  annual_stack_future <- rast(climatological_monthly_means_future)
  mean_annual_temp_future <- app(annual_stack_future, mean, na.rm = TRUE)
  
  # Save the final climatological mean annual temperature
  writeRaster(mean_annual_temp_future,
              paste0(outpath, "/Mean_Annual_Solar_Radiation_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
}



# 9a. MEAN ANNUAL SNOW MELT ---------------------------------------------------
# Surface snow melt flux
# (kg m-2 s-1)
# Mass of ice/snow melted at the surface per second

# These are stored in a MONTHLY format, with daily values per month
# Each monthly raster includes the first day of the next month, which I deal with

if(variable == "mean_melt"){
  
  # Find melt
  variable_names <- files[grepl(pattern = "snm", files)]
  variable_paths <- file_paths[grepl(pattern = "snm", file_paths)]
  
  
  # HISTORICAL --------------------------------------------------------------
  
  # For every year, pull out the monthly rasters, calculate their mean, then summarise for entire year
  
  years_hist <- seq(1995, 2014, by = 1)
  months = seq(1, 12, by=1)
  
  # Get the the days range for each month (what day index is in that month)
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
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month
    # But unlike SICONA - it doesn't appear to actually hold that day in it (so we don't have to remove any layers)
    
    # Get the rasters for the year, but exclude the file that spans from previous year (because of naming)
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_hist[y]),
                                       x = variable_paths, perl = TRUE)]
    year_files <- year_files[!grepl(pattern = paste0((years_hist[y] - 1), "12"), x = year_files)]
    
    print(years_hist[y])
    print(length(year_files))
    
    
    # Load each file then combine
    r_list <- lapply(year_files, function(file) {
      r_temp <- terra::rast(file)
    })
    
    
    # Combine all into one raster
    r <- terra::rast(r_list)
    r
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_hist[y],months[m])
      
      # Subset the relevant layers for the month
      r_month <- r[[Doy]]
      
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
    name <- paste0(outpath, "/Climatological_Monthly_Mean_Melt_",
                   month.name[m], "_1995_2014.tif")
    writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
    
    # Store in list for annual calculation
    climatological_monthly_means[[m]] <- climatological_mean
    
    # Clean up temp files
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 3: Calculate mean annual solar radiation from the 12 climatological monthly means
  annual_stack <- rast(climatological_monthly_means)
  mean_annual_temp <- app(annual_stack, mean, na.rm = TRUE)
  
  # Save the final climatological mean annual temperature
  writeRaster(mean_annual_temp,
              paste0(outpath, "/Mean_Annual_Melt_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # FUTURE CLIMATOLOGY (2081 - 2100) -------------------------------------------
  
  # For every year, pull out the monthly rasters, calculate their mean, then summarise for entire year
  
  years_future <- seq(2081, 2100, by = 1)
  months = seq(1, 12, by=1)
  
  # Get the the days range for each month (what day index is in that month)
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
  
  
  for(y in seq_along(years_future)) {
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month
    # But unlike SICONA - it doesn't appear to actually hold that day in it (so we don't have to remove any layers)
    
    # Get the rasters for the year, but exclude the file that spans from previous year (because of naming)
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_future[y]),
                                       x = variable_paths, perl = TRUE)]
    year_files <- year_files[!grepl(pattern = paste0((years_future[y] - 1), "12"), x = year_files)]
    
    # Load each file then combine
    r_list <- lapply(year_files, function(file) {
      r_temp <- terra::rast(file)
    })
    
    # Combine all into one raster
    r <- terra::rast(r_list)
    r
    
    for(m in seq_along(months)) {
      
      Doy <- get_doy_range(years_future[y],months[m])
      
      # Subset the relevant layers for the month
      r_month <- r[[Doy]]
      
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
  climatological_monthly_means_future <- list()
  
  for(m in seq_along(months)) {
    
    # Stack all years for this month
    month_stack <- rast(monthly_means_by_month[[m]])
    
    # Calculate mean across all years
    climatological_mean <- app(month_stack, mean, na.rm = TRUE)
    
    # Save the climatological monthly mean
    month_name <- sprintf("%02d", m)
    name <- paste0(outpath, "/Climatological_Monthly_Mean_Melt_",
                   month.name[m], "_2081_2100.tif")
    writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
    
    # Store in list for annual calculation
    climatological_monthly_means_future[[m]] <- climatological_mean
    
    # Clean up temp files
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 3: Calculate mean annual melt from the 12 climatological monthly means
  annual_stack_future <- rast(climatological_monthly_means_future)
  mean_annual_temp_future <- app(annual_stack_future, mean, na.rm = TRUE)
  
  # Save the final climatological mean annual temperature
  writeRaster(mean_annual_temp_future,
              paste0(outpath, "/Mean_Annual_Melt_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
}



# 9b. TOTAL ANNUAL SNOW MELT ---------------------------------------------------
# Surface snow melt flux
# (kg m-2 s-1)
# Mass of ice/snow melted at the surface per second

# These are stored in a MONTHLY format, with daily values per month
# Each monthly raster includes the first day of the next month, which I deal with

if(variable == "total_melt"){
  
  # Find melt
  variable_names <- files[grepl(pattern = "snm", files)]
  variable_paths <- file_paths[grepl(pattern = "snm", file_paths)]
  
  
  # HISTORICAL --------------------------------------------------------------
  
  # For every year, pull out the monthly rasters, calculate their total, then summarise for entire year
  
  years_hist <- seq(1995, 2014, by = 1)
  months = seq(1, 12, by=1)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    return(seq(doy_start, doy_end))
  }
  
  
  # STEP 1: Calculate annual total per year
  # Store all annual totals
  
  
  for(y in seq_along(years_hist)) {
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month
    # But unlike SICONA - it doesn't appear to actually hold that day in it (so we don't have to remove any layers)
    
    # Get the rasters for the year, but exclude the file that spans from previous year (because of naming)
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_hist[y]),
                                       x = variable_paths, perl = TRUE)]
    year_files <- year_files[!grepl(pattern = paste0((years_hist[y] - 1), "12"), x = year_files)]
    
    print(years_hist[y])
    print(length(year_files))
    
    # Load each file then combine
    r_list <- lapply(year_files, function(file) {
      r_temp <- terra::rast(file)
    })
    
    
    # Combine all into one raster (for whole year)
    r <- terra::rast(r_list)
    r
    
    # Sum these to get the yearly melt sum:
    yearly_melt_sum <- app(r, sum, na.rm = TRUE)
    
    # Save the yearly sum
    name <- paste0(outpath, "/Total_Annual_Melt_Year_", years_hist[y],".tif")
    
    writeRaster(yearly_melt_sum, name, gdal=c("COMPRESS=NONE"), overwrite = T)
    
    tmp_files <- list.files(tmp_dir, full.names = T, pattern = "^file")
    
    file.remove(tmp_files)
    
  }
  
  # Step 2. MEAN of TOTAL ANNUAL MELT FOR ALL HISTORICAL YEARS ------------------
  
  # List all the yearly mean rasters you just saved
  annual_mean_files <- list.files(outpath, pattern = "^Total_Annual_Melt_Year_.*\\.tif$", full.names = TRUE)
  
  # Subset just the relevant ones based on year:
  # Collapse the years into a regex pattern: "1995|1996|...|2014"
  year_pattern <- paste0("(", paste(years_hist, collapse = "|"), ")")
  
  # Subset the files
  annual_mean_files <- annual_mean_files[grepl(year_pattern, annual_mean_files)]
  
  
  # Stack them together
  annual_means <- rast(annual_mean_files)
  
  # Calculate the overall mean
  final_mean <- app(annual_means, mean, na.rm = TRUE)
  
  # Save the final mean raster
  writeRaster(final_mean, paste0(outpath, "/Total_Annual_Melt_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # FUTURE YEARS ------------------------------------------------------------
  
  years_future <- seq(2081, 2100, by = 1)
  
  months = seq(1, 12, by=1)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    return(seq(doy_start, doy_end))
  }
  
  # STEP 1: Calculate annual total per year
  # Store all annual totals
  
  
  for(y in seq_along(years_future)) {
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month
    # But unlike SICONA - it doesn't appear to actually hold that day in it (so we don't have to remove any layers)
    
    # Get the rasters for the year, but exclude the file that spans from previous year (because of naming)
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_future[y]),
                                       x = variable_paths, perl = TRUE)]
    year_files <- year_files[!grepl(pattern = paste0((years_future[y] - 1), "12"), x = year_files)]
    
    
    print(years_future[y])
    print(length(year_files))
    
    # Load each file then combine
    r_list <- lapply(year_files, function(file) {
      r_temp <- terra::rast(file)
    })
    
    
    # Combine all into one raster (for whole year)
    r <- terra::rast(r_list)
    r
    
    # Sum these to get the yearly melt sum:
    yearly_melt_sum <- app(r, sum, na.rm = TRUE)
    
    # Save the yearly sum
    name <- paste0(outpath, "/Total_Annual_Melt_Year_", years_future[y],".tif")
    
    writeRaster(yearly_melt_sum, name, gdal=c("COMPRESS=NONE"), overwrite = T)
    
    tmp_files <- list.files(tmp_dir, full.names = T, pattern = "^file")
    
    file.remove(tmp_files)
    
  }
  
  # Step 2. MEAN of TOTAL ANNUAL MELT FOR ALL HISTORICAL YEARS ------------------
  
  # List all the yearly mean rasters you just saved
  annual_mean_files <- list.files(outpath, pattern = "^Total_Annual_Melt_Year_.*\\.tif$", full.names = TRUE)
  
  # Subset just the relevant ones based on year:
  # Collapse the years into a regex pattern: "2081|2082|...|2100"
  year_pattern <- paste0("(", paste(years_future, collapse = "|"), ")")
  
  # Subset the files
  annual_mean_files <- annual_mean_files[grepl(year_pattern, annual_mean_files)]
  
  
  # Stack them together
  annual_means <- rast(annual_mean_files)
  
  # Calculate the overall mean
  final_mean <- app(annual_means, mean, na.rm = TRUE)
  
  # Save the final mean raster
  writeRaster(final_mean, paste0(outpath, "/Total_Annual_Melt_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
}



#########################################################################
# 10. MEAN ANNUAL SNOW COVER ---------------------------------------------------
# % "Snow area percentage"

# These are stored in a MONTHLY format, with 6-hourly values per month
# Each monthly raster includes the first 6-hours of the next month, which I deal with

# **ALSO** there was one corrupt file which wasn's able to be downloaded:
# snc_sfx_ANT11_ANT11_ALADIN43_v1_CESM2_r11i1p1f1_historical_6hr_199511010000-199512010000.nc
# NOVEMBER 1995

if(variable == "mean_snow"){
  
  # Find snow cover
  variable_names <- files[grepl(pattern = "snc", files)]
  variable_paths <- file_paths[grepl(pattern = "snc", file_paths)]
  
  
  # HISTORICAL --------------------------------------------------------------
  
  # For every year, pull out the monthly rasters, calculate their mean, then summarise for entire year
  
  years_hist <- seq(1995, 2014, by = 1)
  months = seq(1, 12, by=1)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    # Calculate 6-hourly layer indices
    index_start <- (doy_start - 1) * 4 + 1
    index_end <- doy_end * 4
    
    return(seq(index_start, index_end))
  }
  
  
  
  # STEP 1: Calculate monthly means for each year
  # Store all monthly means organized by month across years
  
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name # month.name is a built in constant
  
  
  for(y in seq_along(years_hist)) {
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month
    # Like SICONA - it appears to hold that day in it (so we have to remove a layer)
    
    # Get the rasters for the year, but exclude the file that spans from previous year (because of naming)
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
    
    # ---- INSERT NA filler for missing rsus November 1995 (31 days) ----
    # This creates an NA raster for the missing month to preserve the right number of days
    if(years_hist[y] == 1995) {
      may_doys <- get_doy_range(1995, 11)          # (30 days)
      filler <- rast(replicate(length(may_doys),  # 30 NA layers matching grid
                               init(r[[1]], fun = NA)))
      
      # Split r at the May insertion point, then reassemble
      r <- rast(list(
        r[[1:(may_doys[1] - 1)]],   # days before May
        filler,                            # NA placeholder for May
        r[[may_doys[1]:nlyr(r)]] # days after May
      ))
    }
    # ---------------------------------------------------------------
    
    
    for(m in seq_along(months)) {
      
      # ---- SKIP November 1995 entirely (data missing for that month) ----
      if(years_hist[y] == 1995 && m == 11) {
        message("Skipping November 1995 - corrupt file")
        monthly_means_by_month[[11]][[y]] <- NA  # NA preserves the slot
        next
      }
      # -------------------------------------------------------------------
      
      
      Doy <- get_doy_range(years_hist[y],months[m])
      
      # Subset the relevant layers for the month
      r_6hr <- r[[Doy]]
      
      #First, need to integrate over the time period to get daily sum
      
      # Number of 4-layer (daily) blocks
      n_days <- nlyr(r_6hr) / 4
      
      # Error check: should be divisible by 4
      if (nlyr(r_6hr) %% 4 != 0) stop("Layer count not divisible by 4")
      
      # Create an index that repeats each group of 4 layers
      index <- rep(1:n_days, each = 4)
      
      # Group every 4 layers and compute daily mean
      r_month <- terra::tapp(r_6hr, index = index, fun = "mean", na.rm = TRUE)
      
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
    
    # Stack all years for this month - *filter out NAs before stacking (due to May 1997 thing above)
    month_stack <- rast(Filter(function(x) inherits(x, "SpatRaster"), monthly_means_by_month[[m]]))
    
    # Calculate mean across all years
    climatological_mean <- app(month_stack, mean, na.rm = TRUE)
    
    # Save the climatological monthly mean
    month_name <- sprintf("%02d", m)
    name <- paste0(outpath, "/Climatological_Monthly_Mean_Snow_Cover_",
                   month.name[m], "_1995_2014.tif")
    writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
    
    # Store in list for annual calculation
    climatological_monthly_means[[m]] <- climatological_mean
    
    # Clean up temp files
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 3: Calculate mean annual snow cover from the 12 climatological monthly means
  annual_stack <- rast(climatological_monthly_means)
  mean_annual_temp <- app(annual_stack, mean, na.rm = TRUE)
  
  # Save the final climatological mean annual temperature
  writeRaster(mean_annual_temp,
              paste0(outpath, "/Mean_Annual_Snow_Cover_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # FUTURE CLIMATOLOGY (2081 - 2100) -------------------------------------------
  
  # For every year, pull out the monthly rasters, calculate their mean, then summarise for entire year
  
  years_future <- seq(2081, 2100, by = 1)
  months = seq(1, 12, by=1)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    # Calculate 8-hourly layer indices
    index_start <- (doy_start - 1) * 4 + 1
    index_end <- doy_end * 4
    
    return(seq(index_start, index_end))
  }
  
  
  
  # STEP 1: Calculate monthly means for each year
  # Store all monthly means organized by month across years
  
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name # month.name is a built in constant
  #
  
  for(y in seq_along(years_future)) {
    
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month
    # Like SICONA - it appears to hold that day in it (so we have to remove a layer)
    
    # Get the rasters for the year, but exclude the file that spans from previous year (because of naming)
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_future[y]),
                                       x = variable_paths, perl = TRUE)]
    year_files <- year_files[!grepl(pattern = paste0((years_future[y] - 1), "12"), x = year_files)]
    
    print(years_future[y])
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
      
      Doy <- get_doy_range(years_future[y],months[m])
      
      # Subset the relevant layers for the month
      r_6hr <- r[[Doy]]
      
      #First, need to integrate over the time period to get daily sum
      
      # Number of 4-layer (daily) blocks
      n_days <- nlyr(r_6hr) / 4
      
      # Error check: should be divisible by 4
      if (nlyr(r_6hr) %% 4 != 0) stop("Layer count not divisible by 4")
      
      # Create an index that repeats each group of 4 layers
      index <- rep(1:n_days, each = 4)
      
      # Group every 4 layers and compute daily mean
      r_month <- terra::tapp(r_6hr, index = index, fun = "mean", na.rm = TRUE)
      
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
  climatological_monthly_means_future <- list()
  
  for(m in seq_along(months)) {
    
    # Stack all years for this month
    month_stack <- rast(monthly_means_by_month[[m]])
    
    # Calculate mean across all years
    climatological_mean <- app(month_stack, mean, na.rm = TRUE)
    
    # Save the climatological monthly mean
    month_name <- sprintf("%02d", m)
    name <- paste0(outpath, "/Climatological_Monthly_Mean_Snow_Cover_",
                   month.name[m], "_2081_2100.tif")
    writeRaster(climatological_mean, name, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
    
    # Store in list for annual calculation
    climatological_monthly_means_future[[m]] <- climatological_mean
    
    # Clean up temp files
    tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
    file.remove(tmp_files)
    
  }
  
  # STEP 3: Calculate mean annual melt from the 12 climatological monthly means
  annual_stack_future <- rast(climatological_monthly_means_future)
  mean_annual_temp_future <- app(annual_stack_future, mean, na.rm = TRUE)
  
  # Save the final climatological mean annual temperature
  writeRaster(mean_annual_temp_future,
              paste0(outpath, "/Mean_Annual_Snow_Cover_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
}


#########################################################################
# 10. MEAN SUMMER SNOW COVER ---------------------------------------------------
# % "Snow area percentage"

# These are stored in a MONTHLY format, with daily values per month
# Each monthly raster includes the first day of the next month, which I deal with

if(variable == "summer_snow"){
  
  # Find snow cover
  variable_names <- files[grepl(pattern = "snc", files)]
  variable_paths <- file_paths[grepl(pattern = "snc", file_paths)]
  
  
  # HISTORICAL --------------------------------------------------------------
  
  # For every year, pull out the monthly rasters, calculate their mean, then summarise for entire year
  
  # Have to re-calculate to include December 1994
  
  years_hist <- seq(1994, 2014, by = 1)  # Include 1994 to get Dec for 1995
  months = seq(1, 12, by=1)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    # Calculate 6-hourly layer indices
    index_start <- (doy_start - 1) * 4 + 1
    index_end <- doy_end * 4
    
    return(seq(index_start, index_end))
  }
  
  
  # STEP 1: Calculate monthly means for each year
  # Store all monthly means organized by month across years
  
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name # month.name is a built in constant
  
  
  for(y in seq_along(years_hist)) {
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month
    # Like SICONA - it appears to hold that day in it (so we have to remove a layer)
    
    # Get the rasters for the year, but exclude the file that spans from previous year (because of naming)
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
      r_6hr <- r[[Doy]]
      
      #First, need to integrate over the time period to get daily sum
      
      # Number of 4-layer (daily) blocks
      n_days <- nlyr(r_6hr) / 4
      
      # Error check: should be divisible by 4
      if (nlyr(r_6hr) %% 4 != 0) stop("Layer count not divisible by 4")
      
      # Create an index that repeats each group of 4 layers
      index <- rep(1:n_days, each = 4)
      
      # Group every 4 layers and compute daily mean
      r_month <- terra::tapp(r_6hr, index = index, fun = "mean", na.rm = TRUE)
      
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
  
  ## December 1994 - 2012 (for summers 1995 - 2013)
  ## January 1995 - 2013 (for summers 1995 - 2013)
  ## February 1995 - 2013 (for summers 1995 - 2013)
  
  # December 1994 - 2013 (for summers 1995 - 2014)
  # January 1995 - 2014 (for summers 1995 - 2014)
  # February 1995 - 2014 (for summers 1995 - 2014)
  
  ## JANUARY
  
  m = 1
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # # Remove final layer (January 2014)
  # month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  # Remove first layer (January 1994)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_JAN <- app(month_stack, mean, na.rm = TRUE)
  
  ## FEBRUARY
  
  m = 2
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # # Remove final layer (February 2014)
  # month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  # Remove first layer (February 1994)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_FEB <- app(month_stack, mean, na.rm = TRUE)
  
  ## DECEMBER
  
  m = 12
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove final layer (December *2014*)
  month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  # # Remove other final layer (December *2013*)
  # month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  
  # Calculate mean across all years
  climatological_mean_DEC <- app(month_stack, mean, na.rm = TRUE)
  
  
  # STEP 3: Calculate mean summer snow cover from the 3 climatological monthly means
  # Summer = DJF (December, January, February)
  
  summer_stack <- c(climatological_mean_DEC, climatological_mean_JAN, climatological_mean_FEB)
  mean_summer_temp <- app(summer_stack, mean, na.rm = TRUE)
  
  
  # Save the climatological mean summer snow cover
  writeRaster(mean_summer_temp,
              paste0(outpath, "/Mean_Summer_Snow_Cover_HISTORICAL_1995_2014.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  # FUTURE --------------------------------------------------------------
  
  # For every year, pull out the monthly rasters, calculate their mean, then summarise for entire year
  
  # Have to re-calculate to include December 2080
  
  years_future <- seq(2080, 2100, by = 1)  # Include 2080 to get Dec for 2081
  months = seq(1, 12, by=1)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(year, month) {
    first_day <- ymd(paste(year, month, "01", sep = "-"))
    last_day <- ymd(paste(year, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    # Calculate 6-hourly layer indices
    index_start <- (doy_start - 1) * 4 + 1
    index_end <- doy_end * 4
    
    return(seq(index_start, index_end))
  }
  
  
  # STEP 1: Calculate monthly means for each year
  # Store all monthly means organized by month across years
  
  monthly_means_by_month <- vector("list", 12)
  names(monthly_means_by_month) <- month.name # month.name is a built in constant
  
  
  for(y in seq_along(years_future)) {
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month
    # Like SICONA - it appears to hold that day in it (so we have to remove a layer)
    
    # Get the rasters for the year, but exclude the file that spans from previous year (because of naming)
    year_files <- variable_paths[grepl(pattern = paste0("(?<=[_-])", years_future[y]),
                                       x = variable_paths, perl = TRUE)]
    year_files <- year_files[!grepl(pattern = paste0((years_future[y] - 1), "12"), x = year_files)]
    
    print(years_future[y])
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
      
      Doy <- get_doy_range(years_future[y],months[m])
      
      # Subset the relevant layers for the month
      r_6hr <- r[[Doy]]
      
      #First, need to integrate over the time period to get daily sum
      
      # Number of 4-layer (daily) blocks
      n_days <- nlyr(r_6hr) / 4
      
      # Error check: should be divisible by 4
      if (nlyr(r_6hr) %% 4 != 0) stop("Layer count not divisible by 4")
      
      # Create an index that repeats each group of 4 layers
      index <- rep(1:n_days, each = 4)
      
      # Group every 4 layers and compute daily mean
      r_month <- terra::tapp(r_6hr, index = index, fun = "mean", na.rm = TRUE)
      
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
  
  # December 2080 - 2099 (for summers 2081 - 2100)
  # January 2081 - 2100 (for summers 2081 - 2100)
  # February 2081 - 2100 (for summers 2081 - 2100)
  
  
  ## JANUARY
  
  m = 1
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # # Remove final layer (January 2100)
  # month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  # Remove first layer (January 2080)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_JAN <- app(month_stack, mean, na.rm = TRUE)
  
  ## FEBRUARY
  
  m = 2
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # # Remove final layer (February 2100)
  # month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  # Remove first layer (February 2080)
  month_stack <- month_stack[[2:nlyr(month_stack)]]
  
  # Calculate mean across all years
  climatological_mean_FEB <- app(month_stack, mean, na.rm = TRUE)
  
  ## DECEMBER
  
  m = 12
  # Stack all years for this month
  month_stack <- rast(monthly_means_by_month[[m]])
  
  # Remove final layer (December *2100*)
  month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  # # Remove other final layer (December *2013*)
  # month_stack <- month_stack[[1:(nlyr(month_stack) - 1)]]
  
  # Calculate mean across all years
  climatological_mean_DEC <- app(month_stack, mean, na.rm = TRUE)
  
  
  # STEP 3: Calculate mean summer snow cover from the 3 climatological monthly means
  # Summer = DJF (December, January, February)
  
  summer_stack <- c(climatological_mean_DEC, climatological_mean_JAN, climatological_mean_FEB)
  mean_summer_temp <- app(summer_stack, mean, na.rm = TRUE)
  
  
  # Save the climatological mean summer snow cover
  writeRaster(mean_summer_temp,
              paste0(outpath, "/Mean_Summer_Snow_Cover_FUTURE_2081_2100.tif"),
              gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  
  
  
}


