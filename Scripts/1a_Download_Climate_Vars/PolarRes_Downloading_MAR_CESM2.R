# HPC version
lib_loc <- paste(getwd(),"/r_lib",sep="")

# extract the arguments provided in the command line
args <- commandArgs(trailingOnly = TRUE)
# The first argument is now the job index
job_index <- as.integer(args[1])

library(terra)
library(glue)
library(gdalraster, lib.loc = lib_loc)
library(vapour, lib.loc = lib_loc)


# MAR MPI-ESM1
u <- "http://ftp.climato.be/fettweis/MARv3.13/PolarRES/Antarctic/MAR-CESM2/"
model_name <- "MAR_CESM2_"

df <- data.frame(file = vsi_read_dir(glue("/vsicurl/{u}")))
df$url <- glue("{u}/{df$file}")


dirpath <- "/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/PolarRes"


# Based on Job Index, assign a year for the job to download
all_years <- c(1:30, 96:116)


target_year <- all_years[job_index]


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
#   ############ NOTE CHANGE VARIABLE NAME EACH TIME ############
#   ltif <- gsub('\\.nc', '_FRA.tif', basename(localfile))
#   
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
#   
#   ############ NOTE CHANGE VARIABLE NAME EACH TIME ############
#   dsn <- sprintf("vrt://%s?sd_name=FRA&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)
#   
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
#   
#   # Delete the .nc file after processing
#   file.remove(localfile)
#   
# 
# }

# FOR TARGET YEAR:

i <- target_year

localfile <-  glue("{dirpath}/{df$file[i]}")
if (!file.exists(localfile)) {
  curl::curl_download(df$url[i],localfile)
}

clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
            #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
            "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")

############ TEMP ############
ltif <- gsub('\\.nc', '_TT.tif', basename(localfile))

localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")

dsn <- sprintf("vrt://%s?sd_name=TT&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)

warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)


############ SWD ############
ltif <- gsub('\\.nc', '_SWD.tif', basename(localfile))

localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")

dsn <- sprintf("vrt://%s?sd_name=SWD&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)

warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)


############ SWU ############
ltif <- gsub('\\.nc', '_SWU.tif', basename(localfile))

localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")

dsn <- sprintf("vrt://%s?sd_name=SWU&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)

warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)


############ PRH ############
ltif <- gsub('\\.nc', '_PRH.tif', basename(localfile))

localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")

dsn <- sprintf("vrt://%s?sd_name=PRH&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)

warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)


############ UU ############
ltif <- gsub('\\.nc', '_UU.tif', basename(localfile))

localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")

dsn <- sprintf("vrt://%s?sd_name=UU&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)

warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)


############ VV ############
ltif <- gsub('\\.nc', '_VV.tif', basename(localfile))

localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")

dsn <- sprintf("vrt://%s?sd_name=VV&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)

warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)


############ FRA ############
ltif <- gsub('\\.nc', '_FRA.tif', basename(localfile))

localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")

dsn <- sprintf("vrt://%s?sd_name=FRA&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)

warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)


############ SHSN3 ############
ltif <- gsub('\\.nc', '_SHSN3.tif', basename(localfile))

localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")

dsn <- sprintf("vrt://%s?sd_name=SHSN3&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)

warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)


####### FINALLY ############

# Delete the .nc file after processing
file.remove(localfile)

gc()

############################################################################

# SECOND 30 YEARS
# 
# for (i in 96:116) {
#   localfile <-  glue("{dirpath}/{df$file[i]}")
#   if (!file.exists(localfile)) {
#     curl::curl_download(df$url[i],localfile)
#   }
# 
#   clargs <- c("-co", "COMPRESS=ZSTD", "-of", "GTiff",
#               #"-multi", "-wo", "NUM_THREADS=ALL_CPUS",
#               "-co", "INTERLEAVE=BAND", "-co", "TILED=NO", "-overwrite")
# 
#   ############ TEMP ############
#   ltif <- gsub('\\.nc', '_TT.tif', basename(localfile))
# 
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
# 
#   dsn <- sprintf("vrt://%s?sd_name=TT&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)
# 
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
# 
# 
#   ############ SWD ############
#   ltif <- gsub('\\.nc', '_SWD.tif', basename(localfile))
# 
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
# 
#   dsn <- sprintf("vrt://%s?sd_name=SWD&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)
# 
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
# 
# 
#   ############ SWU ############
#   ltif <- gsub('\\.nc', '_SWU.tif', basename(localfile))
# 
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
# 
#   dsn <- sprintf("vrt://%s?sd_name=SWU&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)
# 
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
# 
# 
#   ############ PRH ############
#   ltif <- gsub('\\.nc', '_PRH.tif', basename(localfile))
# 
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
# 
#   dsn <- sprintf("vrt://%s?sd_name=PRH&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)
# 
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
# 
# 
#   ############ UU ############
#   ltif <- gsub('\\.nc', '_UU.tif', basename(localfile))
# 
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
# 
#   dsn <- sprintf("vrt://%s?sd_name=UU&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)
# 
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
# 
# 
#   ############ VV ############
#   ltif <- gsub('\\.nc', '_VV.tif', basename(localfile))
# 
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
# 
#   dsn <- sprintf("vrt://%s?sd_name=VV&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)
# 
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
# 
# 
#   ############ FRA ############
#   ltif <- gsub('\\.nc', '_FRA.tif', basename(localfile))
# 
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
# 
#   dsn <- sprintf("vrt://%s?sd_name=FRA&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)
# 
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
# 
# 
#   ############ SHSN3 ############
#   ltif <- gsub('\\.nc', '_SHSN3.tif', basename(localfile))
# 
#   localfile_tif <- glue("{dirpath}/{model_name}tif/{ltif}")
# 
#   dsn <- sprintf("vrt://%s?sd_name=SHSN3&a_srs=+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=20&a_ullr=-4206600,3693600,3462750,-3129300", localfile)
# 
#   warp(dsn, localfile_tif, "EPSG:3031", cl_arg = clargs)
# 
# 
#   ####### FINALLY ############
# 
#   # Delete the .nc file after processing
#   file.remove(localfile)
# 
#   gc()
# 
# }


