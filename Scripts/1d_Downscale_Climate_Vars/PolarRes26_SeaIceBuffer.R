
library(terra)
library(here)

template <- rast(here("siconca_ANT-12_CESM2_historical_r11i1p1f1_UU-IMAU_RACMO24P-NN_v1-r1_day_19850101-19851231.nc"))[[1]]

r <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_MPI_ESM1/Climatological_Monthly_Mean_Sea_Ice_Concentration_November_HISTORICAL_1995_2014_CONCENTRATION.tif"))

plot(r)

r0 <- r<=0
plot(r0)


