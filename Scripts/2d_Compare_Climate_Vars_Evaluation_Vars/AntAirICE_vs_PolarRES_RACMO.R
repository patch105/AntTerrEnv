# VERSION 2. RACMO VARIABLES (RACMO_MPI_ESM1, RACMO_CESM2, RACMO_ERA5) -----

library(terra)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(here)
library(viridis)

# OUTPATH -----------------------------------------------------------------

outpath <- here("Plots/Evaluation_AntAirICE")

# ============================================================
# HELPER: resample ref onto mod's grid (bilinear) if needed
# ============================================================
align_to_model <- function(ref, mod) {
  if (!compareGeom(ref, mod, stopOnError = FALSE)) {
    ref <- resample(ref, mod, method = "bilinear")
  }
  ref
}

# ============================================================
# HELPER: extract non-NA paired values from two rasters
# ============================================================
extract_pairs <- function(ref, mod, label) {
  ref <- align_to_model(ref, mod)
  vals <- data.frame(
    x = as.vector(values(ref)),
    y = as.vector(values(mod))
  )
  vals <- vals[complete.cases(vals), ]
  vals$model <- label
  vals
}

# ============================================================
# LOAD DATA
# ============================================================
AntAir_annual <- rast(here("Data/Environmental_predictors/Mean_Annual_Temp_ICEFREE.tif"))
AntAir_summer <- rast(here("Data/Environmental_predictors/Mean_Summer_Temp_ICEFREE.tif"))

# --- Annual model rasters ---
racmo_mpi_ann   <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_MPI_ESM1/comparison/Mean_Annual_Temperature_HISTORICAL_2003_2014_ICEFREE.tif"))
racmo_cesm2_ann <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_CESM2/comparison/Mean_Annual_Temperature_HISTORICAL_2003_2014_ICEFREE.tif"))
racmo_era5_ann  <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_ERA5/comparison/Mean_Annual_Temperature_HISTORICAL_2003_2014_ICEFREE.tif"))

# --- Summer model rasters ---
racmo_mpi_sum   <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_MPI_ESM1/comparison/Mean_Summer_Temperature_HISTORICAL_2003_2014_ICEFREE.tif"))
racmo_cesm2_sum <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_CESM2/comparison/Mean_Summer_Temperature_HISTORICAL_2003_2014_ICEFREE.tif"))
racmo_era5_sum  <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_ERA5/comparison/Mean_Summer_Temperature_HISTORICAL_2003_2014_ICEFREE.tif"))

# ============================================================
# PLOT 1: Scatter-density panels  (AntAirICE vs each model)
# ============================================================

model_labels <- c("RACMO-MPI-ESM1", "RACMO-CESM2", "RACMO-ERA5")

ann_pairs <- bind_rows(
  extract_pairs(AntAir_annual, racmo_mpi_ann,   "RACMO-MPI-ESM1"),
  extract_pairs(AntAir_annual, racmo_cesm2_ann, "RACMO-CESM2"),
  extract_pairs(AntAir_annual, racmo_era5_ann,  "RACMO-ERA5")
)
ann_pairs$model <- factor(ann_pairs$model, levels = model_labels)

sum_pairs <- bind_rows(
  extract_pairs(AntAir_summer, racmo_mpi_sum,   "RACMO-MPI-ESM1"),
  extract_pairs(AntAir_summer, racmo_cesm2_sum, "RACMO-CESM2"),
  extract_pairs(AntAir_summer, racmo_era5_sum,  "RACMO-ERA5")
)
sum_pairs$model <- factor(sum_pairs$model, levels = model_labels)

# Shared axis limits (use full combined range)
all_vals <- c(ann_pairs$x, ann_pairs$y, sum_pairs$x, sum_pairs$y)
ax_lim <- range(all_vals, na.rm = TRUE)

# Function to build one faceted scatter-density row
scatter_density_row <- function(df, row_title) {
  ggplot(df, aes(x = x, y = y)) +
    stat_bin2d(bins = 150, aes(fill = after_stat(count))) +
    geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.5) +
    scale_fill_viridis_c(
      option  = "plasma",
      name    = "Number of\ngrid cells",
      trans   = "sqrt",          # sqrt compression like the reference
      limits  = c(1, NA),
      na.value= NA
    ) +
    scale_x_continuous(limits = ax_lim, breaks = seq(-30, 10, 10)) +
    scale_y_continuous(limits = ax_lim, breaks = seq(-30, 10, 10)) +
    coord_fixed() +
    facet_wrap(~ model, nrow = 1, strip.position = "top") +
    labs(
      x    = "Temperature (AntAirICE) [°C]",
      y    = paste0("Temperature (Model) [°C]"),
      title= row_title
    ) +
    theme_classic(base_size = 9) +
    theme(
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold", size = 9),
      legend.position   = "right",
      legend.key.height = unit(1.8, "cm"),
      legend.key.width  = unit(0.35, "cm"),
      panel.spacing     = unit(0.4, "cm"),
      plot.title        = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.title.y      = element_text(size = 8),
      axis.title.x      = element_text(size = 8)
    )
}

p_ann <- scatter_density_row(ann_pairs, "Mean annual temperature")
p_sum <- scatter_density_row(sum_pairs, "Mean summer temperature")

plot1 <- p_ann / p_sum +
  plot_annotation(
    title   = "AntAirICE vs RACMO temperatures",
    theme   = theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))
  )

ggsave(
  file.path(outpath, "AntAirICE_vs_RACMO_scatter_density.png"),
  plot1,
  width  = 11,
  height = 7,
  dpi    = 300,
  bg     = "white"
)
message("Plot 1 saved.")

# ============================================================
# PLOT 2: Bias raster maps  (3 columns x 2 rows)
# ============================================================

# Compute bias rasters (resampling AntAirICE onto each model's grid first)
bias_list <- list(
  list(r = racmo_mpi_ann   - align_to_model(AntAir_annual, racmo_mpi_ann),   label = "RACMO-MPI-ESM1", season = "Annual"),
  list(r = racmo_cesm2_ann - align_to_model(AntAir_annual, racmo_cesm2_ann), label = "RACMO-CESM2",    season = "Annual"),
  list(r = racmo_era5_ann  - align_to_model(AntAir_annual, racmo_era5_ann),  label = "RACMO-ERA5",     season = "Annual"),
  list(r = racmo_mpi_sum   - align_to_model(AntAir_summer, racmo_mpi_sum),   label = "RACMO-MPI-ESM1", season = "Summer"),
  list(r = racmo_cesm2_sum - align_to_model(AntAir_summer, racmo_cesm2_sum), label = "RACMO-CESM2",    season = "Summer"),
  list(r = racmo_era5_sum  - align_to_model(AntAir_summer, racmo_era5_sum),  label = "RACMO-ERA5",     season = "Summer")
)

# Convert each bias raster to a dataframe
bias_df <- bind_rows(lapply(bias_list, function(b) {
  df <- as.data.frame(b$r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "bias"
  df$model  <- b$label
  df$season <- b$season
  df
}))

bias_df$model  <- factor(bias_df$model,  levels = model_labels)
bias_df$season <- factor(bias_df$season, levels = c("Annual", "Summer"))

# Symmetric colour scale
bias_lim <- 10   # fixed ±10°C cap — colours use full red-blue range

plot2 <- ggplot(bias_df, aes(x = x, y = y, fill = bias)) +
  geom_tile() +
  scale_fill_distiller(
    palette  = "RdBu",
    limits   = c(-bias_lim, bias_lim),
    name     = "Bias (°C)\n(Model - AntAirICE)",
    direction= -1,
    na.value = NA,
    oob      = scales::squish
  ) +
  coord_fixed() +
  facet_grid(season ~ model) +
  labs(
    title = "Temperature bias: RACMO - AntAirICE",
    x     = "Longitude",
    y     = "Latitude"
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 9),
    legend.position   = "right",
    legend.key.height = unit(2, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    panel.spacing     = unit(0.3, "cm"),
    plot.title        = element_text(face = "bold", size = 11, hjust = 0.5),
    axis.text         = element_text(size = 7),
    axis.title        = element_text(size = 8),
    panel.border      = element_rect(colour = "grey70", fill = NA, linewidth = 0.3)
  )

ggsave(
  file.path(outpath, "AntAirICE_vs_RACMO_bias_maps.png"),
  plot2,
  width  = 11,
  height = 7,
  dpi    = 300,
  bg     = "white"
)
message("Plot 2 saved.")
