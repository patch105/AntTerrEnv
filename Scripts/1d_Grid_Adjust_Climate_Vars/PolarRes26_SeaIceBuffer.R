
library(terra)
library(here)

template <- rast(here("siconca_ANT-12_CESM2_historical_r11i1p1f1_UU-IMAU_RACMO24P-NN_v1-r1_day_19850101-19851231.nc"))[[1]]

landmask <- rast(here("Data/ANT11_masks.nc"))$LSM

exfix <- ext(c(144, 210, -28.1, 25))
crsfix <- "GEOGCRS[\"Rotated_pole\",
    BASEGEOGCRS[\"unknown\",
        DATUM[\"unnamed\",
            ELLIPSOID[\"Sphere\",6371229,0,
                LENGTHUNIT[\"metre\",1,
                    ID[\"EPSG\",9001]]]],
        PRIMEM[\"Greenwich\",0,
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]]],
    DERIVINGCONVERSION[\"Pole rotation (netCDF CF convention)\",
        METHOD[\"Pole rotation (netCDF CF convention)\"],
        PARAMETER[\"Grid north pole latitude (netCDF CF convention)\",5,
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]],
        PARAMETER[\"Grid north pole longitude (netCDF CF convention)\",20,
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]],
        PARAMETER[\"North pole grid longitude (netCDF CF convention)\",0,
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]]],
    CS[ellipsoidal,2],
        AXIS[\"latitude\",north,
            ORDER[1],
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]],
        AXIS[\"longitude\",east,
            ORDER[2],
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]]]"

template_fixed <- template

terra::set.crs(template_fixed, crsfix)
terra::set.ext(template_fixed, exfix)

landmask_fixed <- landmask

terra::set.crs(landmask_fixed, crsfix)

# Align boundary with template
landmask_fixed <- crop(landmask_fixed, ext(template_fixed))

terra::set.ext(landmask_fixed, exfix)

plot(landmask_fixed$LSM)

test <- template_fixed - landmask_fixed


r <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_MPI_ESM1/Climatological_Monthly_Mean_Sea_Ice_Concentration_November_HISTORICAL_1995_2014_CONCENTRATION.tif"))

plot(r)

r0 <- r<=0
plot(r0)


