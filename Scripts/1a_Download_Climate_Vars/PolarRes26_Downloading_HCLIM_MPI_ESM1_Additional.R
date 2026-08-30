# HPC version
lib_loc <- paste(getwd(),"/r_lib_new",sep="")

library(terra)
library(glue, lib.loc = lib_loc)
library(gdalraster, lib.loc = lib_loc)
library(vapour, lib.loc = lib_loc)

# HCLIM_MPI_ESM1 -- ADDITIONAL VARIABLES (rsus, snm, snc, siconca)
# --------------------------------------------------------------------------
# Same source/style as the original HCLIM_MPI_ESM1 download script, but:
#   - no reprojection / no .tif conversion -- the .nc files are kept as-is
#   - everything goes into one flat "additional" folder, not split into
#     historical/ssp370 subfolders (that can be sorted out later -- for now
#     this just gets every matching file downloaded)
#
# Per the subset diagnostics work, MPI's files live directly in the base
# Charlotte/ folder (historical) and its midfut/ subfolder (ssp370) -- there
# is no dedicated MPI subfolder the way CESM2 has a "CESM" one.
u <- c(
  "http://prudence.dmi.dk/data/temp/JAT/Charlotte/",
  "http://prudence.dmi.dk/data/temp/JAT/Charlotte/midfut/"
)

model_name <- "HCLIM_MPI_ESM1_"

dirpath <- "/mnt/hpccs01/home/n11222026/patterc2/AntarcticFutureHabitat/Data/PolarRes26"
outdir  <- glue("{dirpath}/{model_name}additional")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Combined listing across both source folders, restricted to files that are
# actually MPI's (CESM2's files live in some of the same folders)
df_list <- list()
for (u_i in u) {
  files <- vsi_read_dir(glue("/vsicurl/{u_i}"))
  df_list[[u_i]] <- data.frame(file = files, url = glue("{u_i}{files}"))
}
df_all <- do.call(rbind, df_list)
df_all <- df_all[grepl("MPI", df_all$file), ]

# NOTE: if the same filename exists under more than one of the source
# folders above, only the first one encountered is kept -- file.exists()
# below skips anything already downloaded, so a later duplicate is just
# left alone rather than re-downloaded or overwritten.


# RSUS --------------------------------------------------------------------

df <- df_all[grepl("rsus", df_all$file), ]

for (i in seq_len(nrow(df))) {
  localfile <- glue("{outdir}/{df$file[i]}")
  if (!file.exists(localfile)) {
    curl::curl_download(df$url[i], localfile)
  }
}


# SNM (kg m-2 s-1) --------------------------------------------------------------------

df <- df_all[grepl("snm", df_all$file), ]

for (i in seq_len(nrow(df))) {
  localfile <- glue("{outdir}/{df$file[i]}")
  if (!file.exists(localfile)) {
    curl::curl_download(df$url[i], localfile)
  }
}


# SNC (% "Snow area percentage") --------------------------------------------------------------------

df <- df_all[grepl("snc", df_all$file), ]

for (i in seq_len(nrow(df))) {
  localfile <- glue("{outdir}/{df$file[i]}")
  if (!file.exists(localfile)) {
    curl::curl_download(df$url[i], localfile)
  }
}


# SICONCA (%) --------------------------------------------------------------------

df <- df_all[grepl("siconca", df_all$file), ]

for (i in seq_len(nrow(df))) {
  localfile <- glue("{outdir}/{df$file[i]}")
  if (!file.exists(localfile)) {
    curl::curl_download(df$url[i], localfile)
  }
}