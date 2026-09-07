
HCLIM_CESM2_hist <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/HCLIM_CESM2/Mean_Annual_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"))
HCLIM_CESM2_fut <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/HCLIM_CESM2/Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif"))

HCLIM_CESM2_diff <- HCLIM_CESM2_fut-HCLIM_CESM2_hist

HCLIM_MPI_ESM1_hist <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/HCLIM_MPI_ESM1/Mean_Annual_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"))
HCLIM_MPI_ESM1_fut <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/HCLIM_MPI_ESM1/Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif"))

HCLIM_MPI_ESM1_diff <- HCLIM_MPI_ESM1_fut-HCLIM_MPI_ESM1_hist

RACMO_CESM2_hist <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_CESM2/Mean_Annual_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"))
RACMO_CESM2_fut <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_CESM2/Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif"))

RACMO_CESM2_diff <- RACMO_CESM2_fut-RACMO_CESM2_hist

RACMO_MPI_ESM1_hist <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_MPI_ESM1/Mean_Annual_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"))
RACMO_MPI_ESM1_fut <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_MPI_ESM1/Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif"))

RACMO_MPI_ESM1_diff <- RACMO_MPI_ESM1_fut-RACMO_MPI_ESM1_hist



CESM2_diff <- mean(HCLIM_CESM2_diff, RACMO_CESM2_diff)

MPI_ESM1_diff <- mean(HCLIM_MPI_ESM1_diff, RACMO_MPI_ESM1_diff)

writeRaster(CESM2_diff, here("Plots/CESM2_diff.tif"))
writeRaster(MPI_ESM1_diff, here("Plots/MPI_ESM1_diff.tif"))
writeRaster(RACMO_CESM2_diff, here("Plots/RACMO_CESM2_diff.tif"))
writeRaster(RACMO_MPI_ESM1_diff, here("Plots/RACMO_MPI_ESM1_diff.tif"))
writeRaster(HCLIM_CESM2_diff, here("Plots/HCLIM_CESM2_diff.tif"))
writeRaster(HCLIM_MPI_ESM1_diff, here("Plots/HCLIM_MPI_ESM1_diff.tif"))
