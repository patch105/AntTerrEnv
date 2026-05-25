# HPC version
lib_loc <- paste(getwd(),"/r_lib",sep="")

library(terra)
library(glue)
library(gdalraster, lib.loc = lib_loc)
library(vapour, lib.loc = lib_loc)

# HCLIM_MPI_ESM1
u <- "https://ensemblesrt3.dmi.dk/data/prudence/temp/JAT/Charlotte/CESM"
# u <- "https://ensemblesrt3.dmi.dk/data/prudence/temp/JAT/Charlotte/missed"
model_name <- "HCLIM_CESM2_"

dirpath <- "/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes"

# extract the arguments provided in the command line
args <- commandArgs(trailingOnly = TRUE)
# The first argument is now the job index
job_index <- as.integer(args[1])

variables <- list("tas")
# variables <- list("tas", "RSDS", "RSUS", "SNM", "SNC", "sfcWind", "siconca")
# "PR"

variable = variables[[job_index]]


print(paste0("Variable is ", variable))


if(variable == "tas"){
  
  # TAS ---------------------------------------------------------------------
  
  
  df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
  df$url <- glue("{u}/{df$file}")
  
  
  # Subset variable of interest:
  df <- df[grepl("tas", df$file), ]
  
  # system("gdalinfo /vsicurl/https://ensemblesrt3.dmi.dk/data/prudence/temp/JAT/Charlotte/CESM/tas_fp_ANT11_ANT11_ALADIN43_v1_CESM2_r11i1p1f1_ssp370_1hr_209604010000-209605010000.nc")
  
  # Remove files from 1985-1994
  years_to_exclude <- 1985:2092
  pattern_to_exclude <- paste0("_(", paste(years_to_exclude, collapse = "|"), ")")
  df <- df[!grepl(pattern_to_exclude, df$file), ]
  
  for (i in seq_len(nrow(df))) {
    localfile <-  glue("{dirpath}/{df$file[i]}")
    if (!file.exists(localfile)) {
      curl::curl_download(df$url[i],localfile,
                          handle = curl::new_handle(timeout = 600))
    }
    
    # Check if file is valid NetCDF
    header <- readBin(localfile, "raw", n = 4)
    header_str <- rawToChar(header)
    
    if (header_str == "" || !grepl("^CDF", header_str)) {
      message(glue("Skipping i = {i}: Invalid or corrupted file - {df$file[i]}"))
      next
    }
    
    clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
                #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
                "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
    
    ltif <- gsub('\\.nc', '.tif', basename(localfile))
    
    localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
    
    #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
    dsn <- vapour::vapour_sds_names(localfile)[[3]]
    
    warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
    
    
    # Delete the .nc file after processing
    file.remove(localfile)
    
  }
  
  
}





# RSDS --------------------------------------------------------------------

if(variable == "RSDS"){
  
  df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
  df$url <- glue("{u}/{df$file}")
  
  # Subset variable of interest:
  df <- df[grepl("rsds", df$file), ]
  
  for (i in seq_len(nrow(df))) {
    localfile <-  glue("{dirpath}/{df$file[i]}")
    if (!file.exists(localfile)) {
      curl::curl_download(df$url[i],localfile,
                          handle = curl::new_handle(timeout = 600))
    }
    
    # Check if file is valid NetCDF
    header <- readBin(localfile, "raw", n = 4)
    header_str <- rawToChar(header)
    
    if (header_str == "" || !grepl("^CDF", header_str)) {
      message(glue("Skipping i = {i}: Invalid or corrupted file - {df$file[i]}"))
      next
    }
    
    clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
                #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
                "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
    
    ltif <- gsub('\\.nc', '.tif', basename(localfile))
    
    localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
    
    #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
    dsn <- vapour::vapour_sds_names(localfile)[[4]]
    
    warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
    
    
    # Delete the .nc file after processing
    file.remove(localfile)
    
  }
  
}




# RSUS --------------------------------------------------------------------

if(variable == "RSUS"){
  
  df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
  df$url <- glue("{u}/{df$file}")
  
  
  # Subset variable of interest:
  df <- df[grepl("rsus", df$file), ]
  
  # Remove files from 1985-1994
  years_to_exclude <- 1985:1993
  pattern_to_exclude <- paste0("_(", paste(years_to_exclude, collapse = "|"), ")")
  df <- df[!grepl(pattern_to_exclude, df$file), ]
  
  for (i in seq_len(nrow(df))) {
    localfile <-  glue("{dirpath}/{df$file[i]}")
    if (!file.exists(localfile)) {
      curl::curl_download(df$url[i],localfile,
                          handle = curl::new_handle(timeout = 600))
    }
    
    # Check if file is valid NetCDF
    header <- readBin(localfile, "raw", n = 4)
    header_str <- rawToChar(header)
    
    if (header_str == "" || !grepl("^CDF", header_str)) {
      message(glue("Skipping i = {i}: Invalid or corrupted file - {df$file[i]}"))
      next
    }
    
    clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
                #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
                "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
    
    ltif <- gsub('\\.nc', '.tif', basename(localfile))
    
    localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
    
    #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
    dsn <- vapour::vapour_sds_names(localfile)[[4]]
    
    warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
    
    
    # Delete the .nc file after processing
    file.remove(localfile)
    
  }
  
}



# SNM (kg m-2 s-1) --------------------------------------------------------------------

if(variable == "SNM"){
  
  df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
  df$url <- glue("{u}/{df$file}")
  
  
  # Subset variable of interest:
  df <- df[grepl("snm", df$file), ]
  
  
  for (i in seq_len(nrow(df))) {
    localfile <-  glue("{dirpath}/{df$file[i]}")
    if (!file.exists(localfile)) {
      curl::curl_download(df$url[i],localfile,
                          handle = curl::new_handle(timeout = 600))
    }
    
    # Check if file is valid NetCDF
    header <- readBin(localfile, "raw", n = 4)
    header_str <- rawToChar(header)
    
    if (header_str == "" || !grepl("^CDF", header_str)) {
      message(glue("Skipping i = {i}: Invalid or corrupted file - {df$file[i]}"))
      next
    }
    
    clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
                #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
                "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
    
    ltif <- gsub('\\.nc', '.tif', basename(localfile))
    
    localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
    
    #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
    dsn <- vapour::vapour_sds_names(localfile)[[4]]
    
    warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
    
    
    # Delete the .nc file after processing
    file.remove(localfile)
    
  }
  
  
}



# SNC (% "Snow area percentage") --------------------------------------------------------------------

if(variable == "SNC"){
  
  df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
  df$url <- glue("{u}/{df$file}")
  
  
  # Subset variable of interest:
  df <- df[grepl("snc", df$file), ]
  
  # Remove files from 1985-1994
  years_to_exclude <- 1985:1994
  pattern_to_exclude <- paste0("_(", paste(years_to_exclude, collapse = "|"), ")")
  df <- df[!grepl(pattern_to_exclude, df$file), ]
  
  for (i in seq_len(nrow(df))) {
    localfile <-  glue("{dirpath}/{df$file[i]}")
    if (!file.exists(localfile)) {
      curl::curl_download(df$url[i],localfile,
                          handle = curl::new_handle(timeout = 600))
    }
    
    # Check if file is valid NetCDF
    header <- readBin(localfile, "raw", n = 4)
    header_str <- rawToChar(header)
    
    if (header_str == "" || !grepl("^CDF", header_str)) {
      message(glue("Skipping i = {i}: Invalid or corrupted file - {df$file[i]}"))
      next
    }
    
    clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
                #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
                "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
    
    ltif <- gsub('\\.nc', '.tif', basename(localfile))
    
    localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
    
    #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
    dsn <- vapour::vapour_sds_names(localfile)[[3]]
    
    warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
    
    
    # Delete the .nc file after processing
    file.remove(localfile)
    
  }
  
}




# sfcWind (m s-1) --------------------------------------------------------------------

if(variable == "sfcWind"){
  
  df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
  df$url <- glue("{u}/{df$file}")
  
  
  # Subset variable of interest:
  df <- df[grepl("sfcWind", df$file), ]
  
  # # Remove files from 1985-1994
  # years_to_exclude <- 1985:2012
  # 
  # pattern_to_exclude <- paste0("_(", paste(years_to_exclude, collapse = "|"), ")")
  # df <- df[!grepl(pattern_to_exclude, df$file), ]
  
  for (i in seq_len(nrow(df))) {
    localfile <-  glue("{dirpath}/{df$file[i]}")
    if (!file.exists(localfile)) {
      curl::curl_download(df$url[i],localfile,
                          handle = curl::new_handle(timeout = 600))
    }
    
    # Check if file is valid NetCDF
    header <- readBin(localfile, "raw", n = 4)
    header_str <- rawToChar(header)
    
    if (header_str == "" || !grepl("^CDF", header_str)) {
      message(glue("Skipping i = {i}: Invalid or corrupted file - {df$file[i]}"))
      next
    }
    
    clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
                #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
                "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
    
    ltif <- gsub('\\.nc', '.tif', basename(localfile))
    
    localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
    
    #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
    dsn <- vapour::vapour_sds_names(localfile)[[3]]
    
    warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
    
    
    # Delete the .nc file after processing
    file.remove(localfile)
    
  }
  
  
  
}



# PR (kg m-2 s-1 (precipitation flux)) --------------------------------------------------------------------


if(variable == "PR"){

  df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
  df$url <- glue("{u}/{df$file}")


  # Subset variable of interest:
  df <- df[grepl("pr", df$file), ]

  # # Remove files from 2015:2079
  # years_to_exclude <- 2015:2079
  # pattern_to_exclude <- paste0("_(", paste(years_to_exclude, collapse = "|"), ")")
  # df <- df[!grepl(pattern_to_exclude, df$file), ]

  for (i in seq_len(nrow(df))) {
    localfile <-  glue("{dirpath}/{df$file[i]}")
    if (!file.exists(localfile)) {
      curl::curl_download(df$url[i],localfile,
                          handle = curl::new_handle(timeout = 600))
    }

    # Check if file is valid NetCDF
    header <- readBin(localfile, "raw", n = 8)
    header_str <- rawToChar(header[1:3])

    # Check for both NetCDF-3 (CDF) and NetCDF-4/HDF5 (\x89HDF)
    is_netcdf3 <- grepl("^CDF", header_str)
    is_netcdf4 <- header[1] == 0x89 && header[2] == 0x48 && header[3] == 0x44 && header[4] == 0x46

    if (!is_netcdf3 && !is_netcdf4) {
      message(glue("Skipping i = {i}: Invalid or corrupted file - {df$file[i]}"))
      next
    }

    clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
                #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
                "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")

    ltif <- gsub('\\.nc', '.tif', basename(localfile))

    localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")

    #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
    dsn <- vapour::vapour_sds_names(localfile)[[1]]

    warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)


    # Delete the .nc file after processing
    file.remove(localfile)

  }


}



# SICONCA (%) --------------------------------------------------------------------

if(variable == "siconca"){
  
  df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
  df$url <- glue("{u}/{df$file}")
  
  
  # Subset variable of interest:
  df <- df[grepl("siconca", df$file), ]
  
  # Remove files from 1985-1994
  years_to_exclude <- 1985:1993
  pattern_to_exclude <- paste0("_(", paste(years_to_exclude, collapse = "|"), ")")
  df <- df[!grepl(pattern_to_exclude, df$file), ]
  
  for (i in seq_len(nrow(df))) {
    localfile <-  glue("{dirpath}/{df$file[i]}")
    if (!file.exists(localfile)) {
      curl::curl_download(df$url[i],localfile,
                          handle = curl::new_handle(timeout = 600))
    }
    
    # Check if file is valid NetCDF
    header <- readBin(localfile, "raw", n = 4)
    header_str <- rawToChar(header)
    
    if (header_str == "" || !grepl("^CDF", header_str)) {
      message(glue("Skipping i = {i}: Invalid or corrupted file - {df$file[i]}"))
      next
    }
    
    clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
                #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
                "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
    
    ltif <- gsub('\\.nc', '.tif', basename(localfile))
    
    localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
    
    #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
    dsn <- vapour::vapour_sds_names(localfile)[[3]]
    
    warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
    
    
    # Delete the .nc file after processing
    file.remove(localfile)
    
  }
  
}

