# HPC version
lib_loc <- paste(getwd(),"/r_lib",sep="")

library(terra)
library(glue)
library(gdalraster, lib.loc = lib_loc)
library(vapour, lib.loc = lib_loc)

# HCLIM_MPI_ESM1
 u <- "https://ensemblesrt3.dmi.dk/data/prudence/temp/JAT/Charlotte/"
# u <- "https://ensemblesrt3.dmi.dk/data/prudence/temp/JAT/Charlotte/missed"
model_name <- "HCLIM_MPI_ESM1_"

dirpath <- "/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes"



# # TAS ---------------------------------------------------------------------
# 
# 
# df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
# df$url <- glue("{u}/{df$file}")
# 
# 
# # Subset variable of interest:
# df <- df[grepl("tas", df$file), ]
# 
# # system("gdalinfo /vsicurl/https://ensemblesrt3.dmi.dk/data/prudence/temp/JAT/Charlotte/tas_ANT-12_MPI-ESM1-2-LR_ssp370_r1i1p1f1_HCLIMcom-DMI_HCLIM43-ALADIN_v1-r1_day_21000101-21001231.nc")
# # 
# # 
# # system("gdalinfo /vsicurl/https://ensemblesrt3.dmi.dk/data/prudence/temp/JAT/Charlotte/pr_ANT-12_MPI-ESM1-2-LR_historical_r1i1p1f1_HCLIMcom-DMI_HCLIM43-ALADIN_v1-r1_day_19850101-19851231.nc")
# 
# for (i in seq_len(nrow(df))) {
#   localfile <-  glue("{dirpath}/{df$file[i]}")
#   if (!file.exists(localfile)) {
#     curl::curl_download(df$url[i],localfile)
#   }
# 
#   clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
#               #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
#               "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
# 
#   ltif <- gsub('\\.nc', '.tif', basename(localfile))
# 
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
# 
#   #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
#   dsn <- head(vapour::vapour_sds_names(localfile), 1)
# 
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
#   
# 
#   # Delete the .nc file after processing
#   file.remove(localfile)
# 
# }
# 
# 
# 
# RSDS --------------------------------------------------------------------


df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
df$url <- glue("{u}/{df$file}")

# Subset variable of interest:
df <- df[grepl("rsds", df$file), ]

# # Only keep 1994-2014
# years_to_include <- 1994:2014
# pattern_to_include <- paste0("_(", paste(years_to_include, collapse = "|"), ")")
# df <- df[grepl(pattern_to_include, df$file), ]

for (i in seq_len(nrow(df))) {
  localfile <-  glue("{dirpath}/{df$file[i]}")
  if (!file.exists(localfile)) {
    curl::curl_download(df$url[i],localfile)
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

# # RSUS --------------------------------------------------------------------
# 
# 
# df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
# df$url <- glue("{u}/{df$file}")
# 
# 
# # Subset variable of interest:
# df <- df[grepl("rsus", df$file), ]
# 
# 
# for (i in seq_len(nrow(df))) {
#   localfile <-  glue("{dirpath}/{df$file[i]}")
#   if (!file.exists(localfile)) {
#     curl::curl_download(df$url[i],localfile)
#   }
#   
#   clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
#               #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
#               "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
#   
#   ltif <- gsub('\\.nc', '.tif', basename(localfile))
#   
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
#   
#   #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
#   dsn <- vapour::vapour_sds_names(localfile)[[4]]
#   
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
#   
#   
#   # Delete the .nc file after processing
#   file.remove(localfile)
#   
# }

# # SNM (kg m-2 s-1) --------------------------------------------------------------------
# 
# 
# df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
# df$url <- glue("{u}/{df$file}")
# 
# 
# # Subset variable of interest:
# df <- df[grepl("snm", df$file), ]
# 
# 
# for (i in seq_len(nrow(df))) {
#   localfile <-  glue("{dirpath}/{df$file[i]}")
#   if (!file.exists(localfile)) {
#     curl::curl_download(df$url[i],localfile)
#   }
#   
#   clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
#               #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
#               "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
#   
#   ltif <- gsub('\\.nc', '.tif', basename(localfile))
#   
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
#   
#   #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
#   dsn <- vapour::vapour_sds_names(localfile)[[4]]
#   
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
#   
#   
#   # Delete the .nc file after processing
#   file.remove(localfile)
#   
# }


# # SNC (% "Snow area percentage") --------------------------------------------------------------------
# 
# 
# df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
# df$url <- glue("{u}/{df$file}")
# 
# 
# # Subset variable of interest:
# df <- df[grepl("snc", df$file), ]
# 
# 
# for (i in seq_len(nrow(df))) {
#   localfile <-  glue("{dirpath}/{df$file[i]}")
#   if (!file.exists(localfile)) {
#     curl::curl_download(df$url[i],localfile)
#   }
#   
#   clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
#               #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
#               "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
#   
#   ltif <- gsub('\\.nc', '.tif', basename(localfile))
#   
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
#   
#   #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
#   dsn <- vapour::vapour_sds_names(localfile)[[3]]
#   
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
#   
#   
#   # Delete the .nc file after processing
#   file.remove(localfile)
#   
# }


# # sfcWind (m s-1) --------------------------------------------------------------------
# 
# 
# df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
# df$url <- glue("{u}/{df$file}")
# 
# 
# # Subset variable of interest:
# df <- df[grepl("sfcWind", df$file), ]
# 
# 
# for (i in seq_len(nrow(df))) {
#   localfile <-  glue("{dirpath}/{df$file[i]}")
#   if (!file.exists(localfile)) {
#     curl::curl_download(df$url[i],localfile)
#   }
#   
#   clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
#               #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
#               "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
#   
#   ltif <- gsub('\\.nc', '.tif', basename(localfile))
#   
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
#   
#   #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
#   dsn <- vapour::vapour_sds_names(localfile)[[1]]
#   
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
#   
#   
#   # Delete the .nc file after processing
#   file.remove(localfile)
#   
# }


# # PR (kg m-2 s-1 (precipitation flux)) --------------------------------------------------------------------
# 
# 
# df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
# df$url <- glue("{u}/{df$file}")
# 
# 
# # Subset variable of interest:
# df <- df[grepl("pr", df$file), ]
# 
# 
# for (i in seq_len(nrow(df))) {
#   localfile <-  glue("{dirpath}/{df$file[i]}")
#   if (!file.exists(localfile)) {
#     curl::curl_download(df$url[i],localfile)
#   }
#   
#   clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
#               #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
#               "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
#   
#   ltif <- gsub('\\.nc', '.tif', basename(localfile))
#   
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
#   
#   #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
#   dsn <- vapour::vapour_sds_names(localfile)[[1]]
#   
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
#   
#   
#   # Delete the .nc file after processing
#   file.remove(localfile)
#   
# }


# # SICONCA (%) --------------------------------------------------------------------
# 
# 
# df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
# df$url <- glue("{u}/{df$file}")
# 
# 
# # Subset variable of interest:
# df <- df[grepl("siconca", df$file), ]
# 
# 
# for (i in seq_len(nrow(df))) {
#   localfile <-  glue("{dirpath}/{df$file[i]}")
#   if (!file.exists(localfile)) {
#     curl::curl_download(df$url[i],localfile)
#   }
#   
#   clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
#               #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
#               "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
#   
#   ltif <- gsub('\\.nc', '.tif', basename(localfile))
#   
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
#   
#   #Tell it to look at the last band of the file (e.g., "surface_upward_latent_heat_flux)
#   dsn <- vapour::vapour_sds_names(localfile)[[3]]
#   
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
#   
#   
#   # Delete the .nc file after processing
#   file.remove(localfile)
#   
# }