

#############################################################
# PolarRES DAILY values for AntAirICE and AWS evaluations -----------------
#############################################################

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

# models     <- list("MAR_MPI_ESM1", "MAR_CESM2", "HCLIM_MPI_ESM1", "HCLIM_CESM2")
models     <- list("HCLIM_CESM2")
years_hist <- seq(1994, 2014, by = 1)

# Build a grid of all model × year combinations
job_grid  <- expand.grid(model = unlist(models), year = years_hist, stringsAsFactors = FALSE)

model <- job_grid$model[job_index]
year  <- job_grid$year[job_index]

# Choose variable
variable <- "tas"

print(paste0("Model: ", model, " | Year: ", year, " | Variable: ", variable))


# -------------------------------------------------------------------------
# Load domain (ice-free mask)
# -------------------------------------------------------------------------
domain <- rast(here("Data/Environmental_predictors/ice_free_upsamp_1km.tif"))

# Set domain values
domain <- ifel(not.na(domain), 1, NA)

tmp_dir <- tempdir()

clean_tmp <- function(tmp_dir) {
  tmp_files <- list.files(tmp_dir, full.names = TRUE, pattern = "^file")
  file.remove(tmp_files)
}

# -------------------------------------------------------------------------
# Step 1: Load and compute daily rasters for this model × year
# -------------------------------------------------------------------------

if (model == "MAR_MPI_ESM1") {
  
  file_paths <- list.files(
    "/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes/MAR_MPI_ESM1_tif",
    full.names = TRUE, recursive = TRUE
  )
  outpath <- here("Data/Environmental_predictors/PolarRes/MAR_MPI_ESM1/Validation")
  
  if(variable == "tas"){
    
    variable_paths <- file_paths[grepl("TT", file_paths)]
    year_paths     <- variable_paths[grepl(as.character(year), variable_paths)]
    
    r      <- terra::rast(year_paths)
    n_days <- nlyr(r) / 3
    if (nlyr(r) %% 3 != 0) stop(paste("Year", year, ": layer count not divisible by 3"))
    
    index   <- rep(1:n_days, each = 3)
    r_daily <- terra::tapp(r, index = index, fun = "mean", na.rm = TRUE)
    clean_tmp(tmp_dir)
    
  }
  
  if(variable == "wind"){
    
    # Find U and V component wind paths
    variable_paths  <- file_paths[grepl("UU", file_paths)]
    variable_paths2 <- file_paths[grepl("VV", file_paths)]
    year_paths      <- variable_paths[grepl(as.character(year), variable_paths)]
    year_paths2     <- variable_paths2[grepl(as.character(year), variable_paths2)]  
    
    UU <- terra::rast(year_paths)
    VV <- terra::rast(year_paths2)
    
    n_days <- nlyr(UU) / 3
    if (nlyr(UU) %% 3 != 0) stop(paste("Year", year, ": layer count not divisible by 3"))
    
    index <- rep(1:n_days, each = 3)
    
    UU_daily <- terra::tapp(UU, index = index, fun = "mean", na.rm = TRUE)
    VV_daily <- terra::tapp(VV, index = index, fun = "mean", na.rm = TRUE)
    
    r_daily <- sqrt((UU_daily^2) + (VV_daily^2))
    clean_tmp(tmp_dir)
    
  }
  
  
} else if (model == "MAR_CESM2") {
  
  file_paths <- list.files(
    "/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes/MAR_CESM2_tif",
    full.names = TRUE, recursive = TRUE
  )
  outpath <- here("Data/Environmental_predictors/PolarRes/MAR_CESM2/Validation")
  
  if(variable == "tas"){
    
    variable_paths <- file_paths[grepl("TT", file_paths)]
    year_paths     <- variable_paths[grepl(as.character(year), variable_paths)]
    
    r      <- terra::rast(year_paths)
    n_days <- nlyr(r) / 3
    if (nlyr(r) %% 3 != 0) stop(paste("Year", year, ": layer count not divisible by 3"))
    
    index   <- rep(1:n_days, each = 3)
    r_daily <- terra::tapp(r, index = index, fun = "mean", na.rm = TRUE)
    clean_tmp(tmp_dir)
    
  }
  
  if(variable == "wind"){
    
    # Find U and V component wind paths
    variable_paths  <- file_paths[grepl("UU", file_paths)]
    variable_paths2 <- file_paths[grepl("VV", file_paths)]
    year_paths      <- variable_paths[grepl(as.character(year), variable_paths)]
    year_paths2     <- variable_paths2[grepl(as.character(year), variable_paths2)]
    
    UU <- terra::rast(year_paths)
    VV <- terra::rast(year_paths2)
    
    n_days <- nlyr(UU) / 3
    if (nlyr(UU) %% 3 != 0) stop(paste("Year", year, ": layer count not divisible by 3"))
    
    index <- rep(1:n_days, each = 3)
    
    UU_daily <- terra::tapp(UU, index = index, fun = "mean", na.rm = TRUE)
    VV_daily <- terra::tapp(VV, index = index, fun = "mean", na.rm = TRUE)
    
    r_daily <- sqrt((UU_daily^2) + (VV_daily^2))
    clean_tmp(tmp_dir)
    
    
  }
  
  
  
} else if (model == "HCLIM_MPI_ESM1") {
  
  file_paths <- list.files(
    "/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes/HCLIM_MPI_ESM1_tif",
    full.names = TRUE, recursive = TRUE
  )
  outpath <- here("Data/Environmental_predictors/PolarRes/HCLIM_MPI_ESM1/Validation")
  
  if(variable == "tas"){
    
    variable_paths <- file_paths[grepl("tas", file_paths)]
    year_paths     <- variable_paths[grepl(
      pattern = paste0("(?<=[_-])", year),
      x       = variable_paths,
      perl    = TRUE
    )]
    
    r_daily <- terra::rast(year_paths)
    r_daily  <- r_daily - 273.15
    clean_tmp(tmp_dir)
    
  }
  
  if(variable == "wind"){
    
    variable_paths <- file_paths[grepl("sfcWind", file_paths)]
    year_paths     <- variable_paths[grepl(
      pattern = paste0("(?<=[_-])", year),
      x       = variable_paths,
      perl    = TRUE
    )]
    
    r_daily <- terra::rast(year_paths)
    clean_tmp(tmp_dir)
    
    
  }
  
  
} else if (model == "HCLIM_CESM2") {
  
  
  library(lubridate)
  
  # Get the the days range for each month (what day index is in that month)
  get_doy_range <- function(yr, month) {
    first_day <- ymd(paste(yr, month, "01", sep = "-"))
    last_day <- ymd(paste(yr, month, days_in_month(first_day), sep = "-"))
    
    doy_start <- yday(first_day)
    doy_end <- yday(last_day)
    
    # Calculate Hourly layer indices
    index_start <- (doy_start - 1) * 24 + 1
    index_end <- doy_end * 24
    
    return(seq(index_start, index_end))
  }
  
  file_paths <- list.files(
    "/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes/HCLIM_CESM2_tif",
    full.names = TRUE, recursive = TRUE
  )
  outpath <- here("Data/Environmental_predictors/PolarRes/HCLIM_CESM2/Validation")
  
  if(variable == "tas"){
    
    variable_paths <- file_paths[grepl("tas", file_paths)]
    year_paths     <- variable_paths[grepl(
      pattern = paste0("(?<=[_-])", year),
      x       = variable_paths,
      perl    = TRUE
    )]
    
    # *NOTE* have to do some fiddly stuff since each monthly raster includes the name of the first day of the next month
    year_paths <- year_paths[!grepl(pattern = paste0((year - 1), "12"), x = year_paths)]
    
    
    # ── DIAGNOSTIC 1: What files were matched for this year? ──────────────────
    message("=== DIAGNOSTICS FOR HCLIM_CESM2 | Year: ", year, " ===")
    message("Number of files matched: ", length(year_paths))
    message("Files matched (in order):")
    for (f in year_paths) message("  ", f)
    
    # ── DIAGNOSTIC 2: Layer count per file, before any trimming ───────────────
    layers_per_file <- sapply(year_paths, function(f) terra::nlyr(terra::rast(f)))
    message("\nLayers per file (before trimming):")
    for (i in seq_along(year_paths)) {
      message("  File ", i, ": ", layers_per_file[i], " layers | ", basename(year_paths[i]))
    }
    message("Total layers across all files (before trimming): ", sum(layers_per_file))
    
    # ── DIAGNOSTIC 3: What does -1 per file actually remove in total? ──────────
    total_removed <- length(year_paths)  # one layer dropped per file
    total_after_trim <- sum(layers_per_file) - total_removed
    message("\nLayers removed by [[1:(nlyr-1)]] logic: ", total_removed, " (one per file)")
    message("Total layers after per-file trimming: ", total_after_trim)
    
    r_list <- lapply(year_paths, function(file) {
      r_temp <- terra::rast(file)
      
      # ── DIAGNOSTIC 4: First and last layer name/time per file ───────────────
      message("\n  File: ", basename(file))
      message("    nlyr         : ", terra::nlyr(r_temp))
      message("    First layer  : ", names(r_temp)[1])
      message("    Last layer   : ", names(r_temp)[terra::nlyr(r_temp)])
      
      # Check if time metadata exists
      t <- try(terra::time(r_temp), silent = TRUE)
      if (!inherits(t, "try-error") && !is.null(t)) {
        message("    First time   : ", t[1])
        message("    Last time    : ", t[length(t)])
      } else {
        message("    Time metadata: not available")
      }
      
      r_trimmed <- r_temp[[1:(terra::nlyr(r_temp) - 1)]]
      message("    Layers after trim: ", terra::nlyr(r_trimmed))
      r_trimmed
    })
    
    r <- terra::rast(r_list)
    
    # ---- INSERT NA filler for missing rsus June 1994 (31 days) ----
    # This creates an NA raster for the missing month to preserve the right number of days
    if(year == 1994) {
      june_doys <- get_doy_range(1994, 6)           # days 152-181 (30 days)
      
      filler_june <- rast(replicate(length(june_doys),  # 30 NA layers matching grid
                                    init(r[[1]], fun = NA)))
      
      # Split r at the June insertion point, then reassemble
      r_1hr <- rast(list(
        r[[1:(june_doys[1] - 1)]],   # days before June
        filler_june,                            # NA placeholder for June
        r[[june_doys[1]:nlyr(r)]] # days after June
      ))
    } else {r_1hr <- r}
    
    # ── DIAGNOSTIC 5: Total after stacking, vs expected ───────────────────────
    expected_hours <- ifelse(lubridate::leap_year(year), 366, 365) * 24
    message("\nTotal layers after stacking trimmed files : ", terra::nlyr(r_1hr))
    message("Expected hours for year ", year, "          : ", expected_hours)
    message("Difference (actual - expected)            : ", terra::nlyr(r_1hr) - expected_hours)
    
    # ── DIAGNOSTIC 6: Does total divide evenly into days? ─────────────────────
    remainder <- terra::nlyr(r_1hr) %% 24
    message("Remainder when dividing stacked layers by 24: ", remainder)
    if (remainder != 0) {
      message("WARNING: Layer count is NOT divisible by 24 — daily means will fail.")
    } else {
      message("Layer count IS divisible by 24 (", terra::nlyr(r_1hr) / 24, " days).")
    }
    
    n_days <- nlyr(r_1hr) / 24
    if (nlyr(r_1hr) %% 24 != 0) stop(paste("Year", year, ": layer count not divisible by 24"))
    
    index   <- rep(1:n_days, each = 24)
    r_daily <- terra::tapp(r_1hr, index = index, fun = "mean", na.rm = TRUE)
    r_daily  <- r_daily - 273.15
    clean_tmp(tmp_dir)
    
  }
  
  
  if(variable == "wind"){
    
    variable_paths <- file_paths[grepl("sfcWind", file_paths)]
    year_paths     <- variable_paths[grepl(
      pattern = paste0("(?<=[_-])", year),
      x       = variable_paths,
      perl    = TRUE
    )]
    year_paths <- year_paths[!grepl(pattern = paste0((year - 1), "12"), x = year_paths)]
    
    # ── DIAGNOSTIC 1: What files were matched for this year? ──────────────────
    message("=== DIAGNOSTICS FOR HCLIM_CESM2 wind | Year: ", year, " ===")
    message("Number of files matched: ", length(year_paths))
    message("Files matched (in order):")
    for (f in year_paths) message("  ", f)
    
    # ── DIAGNOSTIC 2: Layer count per file, before any trimming ───────────────
    layers_per_file <- sapply(year_paths, function(f) terra::nlyr(terra::rast(f)))
    message("\nLayers per file (before trimming):")
    for (i in seq_along(year_paths)) {
      message("  File ", i, ": ", layers_per_file[i], " layers | ", basename(year_paths[i]))
    }
    message("Total layers across all files (before trimming): ", sum(layers_per_file))
    
    # ── DIAGNOSTIC 3: What does -1 per file actually remove in total? ──────────
    total_removed <- length(year_paths)
    total_after_trim <- sum(layers_per_file) - total_removed
    message("\nLayers removed by [[1:(nlyr-1)]] logic: ", total_removed, " (one per file)")
    message("Total layers after per-file trimming: ", total_after_trim)
    
    r_list <- lapply(year_paths, function(file) {
      r_temp <- terra::rast(file)
      
      # ── DIAGNOSTIC 4: First and last layer name/time per file ───────────────
      message("\n  File: ", basename(file))
      message("    nlyr         : ", terra::nlyr(r_temp))
      message("    First layer  : ", names(r_temp)[1])
      message("    Last layer   : ", names(r_temp)[terra::nlyr(r_temp)])
      
      t <- try(terra::time(r_temp), silent = TRUE)
      if (!inherits(t, "try-error") && !is.null(t)) {
        message("    First time   : ", t[1])
        message("    Last time    : ", t[length(t)])
      } else {
        message("    Time metadata: not available")
      }
      
      r_trimmed <- r_temp[[1:(terra::nlyr(r_temp) - 1)]]
      message("    Layers after trim: ", terra::nlyr(r_trimmed))
      r_trimmed
    })
    
    r_1hr <- terra::rast(r_list)
    
    # ── DIAGNOSTIC 5: Total after stacking, vs expected ───────────────────────
    expected_hours <- ifelse(lubridate::leap_year(year), 366, 365) * 24
    message("\nTotal layers after stacking trimmed files : ", terra::nlyr(r_1hr))
    message("Expected hours for year ", year, "          : ", expected_hours)
    message("Difference (actual - expected)            : ", terra::nlyr(r_1hr) - expected_hours)
    
    # ── DIAGNOSTIC 6: Does total divide evenly into days? ─────────────────────
    remainder <- terra::nlyr(r_1hr) %% 24
    message("Remainder when dividing stacked layers by 24: ", remainder)
    if (remainder != 0) {
      message("WARNING: Layer count is NOT divisible by 24 — daily means will fail.")
    } else {
      message("Layer count IS divisible by 24 (", terra::nlyr(r_1hr) / 24, " days).")
    }
    
    n_days <- nlyr(r_1hr) / 24
    if (nlyr(r_1hr) %% 24 != 0) stop(paste("Year", year, ": layer count not divisible by 24"))
    
    index   <- rep(1:n_days, each = 24)
    r_daily <- terra::tapp(r_1hr, index = index, fun = "mean", na.rm = TRUE)
    clean_tmp(tmp_dir)
    
  }
  
  
}
  


# -------------------------------------------------------------------------
# Step 2: Reproject → mask to ice-free domain
# -------------------------------------------------------------------------
r_daily <- terra::project(r_daily, domain, method = "near")
r_daily <- mask(r_daily, domain, maskvalue = NA)

# -------------------------------------------------------------------------
# Step 3: Save — all daily layers for this year in one file
# -------------------------------------------------------------------------
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

if(variable == "tas"){
  
  out_file <- file.path(
    outpath,
    paste0("Daily_Temperature_", model, "_", year, "_ICEFREE.tif")
  )
  
}

if(variable == "wind"){
  
  out_file <- file.path(
    outpath,
    paste0("Daily_Wind_Speed_", model, "_", year, "_ICEFREE.tif")
  )
  
}

writeRaster(r_daily, out_file, overwrite = TRUE)
message("Saved: ", out_file)