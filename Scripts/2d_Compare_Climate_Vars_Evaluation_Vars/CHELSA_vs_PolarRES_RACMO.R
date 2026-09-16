# ============================================================
# CHELSA vs RACMO — spatial map comparison
# Same structure as the CHELSA vs HCLIM script, updated for RACMO
# (RACMO_MPI_ESM1, RACMO_CESM2, RACMO_ERA5), reading from the
# Regridded/<model>/comparison/ outputs, historical period 1985-2010.
# ============================================================

library(terra)
library(ggplot2)
library(dplyr)
library(patchwork)
library(here)
library(scales)
library(sf)

# ============================================================
# OUTPUT PATH
# ============================================================

outpath <- here("Plots/Evaluation_CHELSA")

if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}

# ============================================================
# SHARED SIZING / HELPERS (defined once, used for both sections)
# ============================================================

PAGE_W   <- 8.27
PAGE_H   <- 11.69
DPI      <- 300
BASE     <- 14

LEGBAR_H <- unit(3.5, "cm")
LEGBAR_W <- unit(0.45, "cm")

rast_to_df <- function(r, val_name = "value") {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- val_name
  df
}

theme_ant <- function(base_size = BASE) {
  theme_void(base_size = base_size) +
    theme(
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      legend.position   = "right",
      legend.title      = element_text(size = base_size, angle = 90,
                                       hjust = 0.5, vjust = 0.5),
      legend.text       = element_text(size = base_size - 1),
      legend.key.height = LEGBAR_H,
      legend.key.width  = LEGBAR_W,
      legend.margin     = margin(0, 0, 0, 6),
      plot.title        = element_text(size = base_size, face = "bold",
                                       hjust = 0, vjust = 1,
                                       margin = margin(b = 4)),
      plot.margin       = margin(5, 5, 5, 5)
    )
}

make_cbar <- function() {
  guide_colorbar(
    title.position  = "right",
    title.hjust     = 0.5,
    barheight       = LEGBAR_H,
    barwidth        = LEGBAR_W,
    ticks.colour    = "black",
    frame.colour    = "black",
    frame.linewidth = 0.4
  )
}

# Diverging colour ramp shared by both the temperature scale and the
# difference scale (white placed at 0 via `values`, computed per-section)
diverging_ramp <- c(
  "#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
  "white",
  "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"
)

make_scale_values <- function(min_val, max_val) {
  white_pos <- (0 - min_val) / (max_val - min_val)
  c(
    0,
    white_pos * 0.25, white_pos * 0.50, white_pos * 0.75, white_pos * 0.95,
    white_pos,
    white_pos + (1 - white_pos) * 0.15,
    white_pos + (1 - white_pos) * 0.40,
    white_pos + (1 - white_pos) * 0.65,
    white_pos + (1 - white_pos) * 0.85,
    1
  )
}

make_temp_panel <- function(r, panel_label, temp_min, temp_max, temp_values,
                            coast, base_size = BASE, xlim = NULL, ylim = NULL) {
  df <- rast_to_df(r, "temp")
  p <- ggplot() +
    geom_sf(data = coast, fill = "grey85", colour = "black", linewidth = 0.15) +
    geom_tile(data = df, aes(x = x, y = y, fill = temp)) +
    scale_fill_gradientn(
      colours = diverging_ramp,
      values  = temp_values,
      limits  = c(temp_min, temp_max),
      breaks  = pretty(c(temp_min, temp_max), n = 8),
      name    = "Temperature (°C)",
      oob     = squish,
      guide   = make_cbar()
    ) +
    labs(title = panel_label) +
    theme_ant(base_size)
  
  if (!is.null(xlim) && !is.null(ylim)) {
    p + coord_sf(xlim = xlim, ylim = ylim, expand = FALSE)
  } else {
    p + coord_sf(expand = FALSE)
  }
}

make_diff_panel <- function(r, panel_label, diff_min, diff_max, diff_breaks,
                            diff_values, coast, base_size = BASE,
                            xlim = NULL, ylim = NULL) {
  df <- rast_to_df(r, "diff")
  p <- ggplot() +
    geom_sf(data = coast, fill = "grey85", colour = "black", linewidth = 0.15) +
    geom_tile(data = df, aes(x = x, y = y, fill = diff)) +
    scale_fill_gradientn(
      colours = diverging_ramp,
      values  = diff_values,
      limits  = c(diff_min, diff_max),
      breaks  = diff_breaks,
      name    = "Difference (°C)\n(RACMO - CHELSA)",
      oob     = squish,
      guide   = make_cbar()
    ) +
    labs(title = panel_label) +
    theme_ant(base_size)
  
  if (!is.null(xlim) && !is.null(ylim)) {
    p + coord_sf(xlim = xlim, ylim = ylim, expand = FALSE)
  } else {
    p + coord_sf(expand = FALSE)
  }
}

# Assemble CHELSA row + model rows (temp | spacer | diff, per model) and
# save. Generalised to any number of models -- the model-row weight scales
# with n_models so a 3-model RACMO page isn't squeezed the way a fixed
# weight (tuned for 2 models) would squeeze it.
assemble_and_save <- function(p_chelsa, model_panels, out_file,
                              width = PAGE_W, height = PAGE_H, dpi = DPI) {
  
  n_models <- length(model_panels) / 2
  
  chelsa_row <- plot_spacer() + p_chelsa + plot_spacer() +
    plot_layout(ncol = 3, widths = c(0.5, 2, 0.5))
  
  model_panels_spaced <- list()
  for (i in seq(1, length(model_panels), by = 2)) {
    model_panels_spaced <- c(
      model_panels_spaced,
      list(model_panels[[i]], plot_spacer(), model_panels[[i + 1]])
    )
  }
  
  model_grid <- wrap_plots(model_panels_spaced, ncol = 3) +
    plot_layout(guides = "keep", widths = c(1, 0.06, 1))
  
  final_plot <- chelsa_row / model_grid +
    plot_layout(heights = c(1.1, 2 * n_models))
  
  ggsave(out_file, final_plot, width = width, height = height,
         dpi = dpi, bg = "white")
  
  message("Saved: ", out_file)
}

# ============================================================
# COASTLINE (shared)
# ============================================================

coast <- st_read(
  here("Data/add_coastline_medium_res_polygon_v7_10.shp"),
  quiet = TRUE
)

# ============================================================
# MODEL LABELS
# ============================================================

model_labels <- c("RACMO-MPI-ESM1", "RACMO-CESM2", "RACMO-ERA5")

# ============================================================================
# 1. HISTORICAL (1981-2010 CHELSA vs 1985-2010 RACMO)
# ============================================================================

CHELSA_hist <- rast(
  here("Data/CHELSA/Mean_Annual_Temperature_HISTORICAL_1981_2010_ICEFREE.tif")
)

racmo_mpi_hist <- rast(
  here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_MPI_ESM1/comparison/Mean_Annual_Temperature_HISTORICAL_1985_2010_ICEFREE.tif")
)

racmo_cesm_hist <- rast(
  here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_CESM2/comparison/Mean_Annual_Temperature_HISTORICAL_1985_2010_ICEFREE.tif")
)

racmo_era5_hist <- rast(
  here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_ERA5/comparison/Mean_Annual_Temperature_HISTORICAL_1985_2010_ICEFREE.tif")
)

# --- Align to model grid + compute bias (RACMO - CHELSA) ---

align_to_model <- function(ref, mod) {
  if (!compareGeom(ref, mod, stopOnError = FALSE)) {
    ref <- resample(ref, mod, method = "bilinear")
  }
  ref
}

RACMO_MPI_hist_diff   <- racmo_mpi_hist  - align_to_model(CHELSA_hist, racmo_mpi_hist)
RACMO_CESM2_hist_diff <- racmo_cesm_hist - align_to_model(CHELSA_hist, racmo_cesm_hist)
RACMO_ERA5_hist_diff  <- racmo_era5_hist - align_to_model(CHELSA_hist, racmo_era5_hist)

# --- Aggregate everything to 10 km for mapping ---

CHELSA_hist_10km          <- aggregate(CHELSA_hist, fact = 10, fun = "mean", na.rm = TRUE)
racmo_mpi_hist_10km       <- aggregate(racmo_mpi_hist, fact = 10, fun = "mean", na.rm = TRUE)
racmo_cesm_hist_10km      <- aggregate(racmo_cesm_hist, fact = 10, fun = "mean", na.rm = TRUE)
racmo_era5_hist_10km      <- aggregate(racmo_era5_hist, fact = 10, fun = "mean", na.rm = TRUE)
RACMO_MPI_hist_diff_10km   <- aggregate(RACMO_MPI_hist_diff, fact = 10, fun = "mean", na.rm = TRUE)
RACMO_CESM2_hist_diff_10km <- aggregate(RACMO_CESM2_hist_diff, fact = 10, fun = "mean", na.rm = TRUE)
RACMO_ERA5_hist_diff_10km  <- aggregate(RACMO_ERA5_hist_diff, fact = 10, fun = "mean", na.rm = TRUE)

# --- Colour scales ---

chelsa_vals <- values(CHELSA_hist_10km, na.rm = TRUE)
temp_min    <- floor(min(chelsa_vals))
temp_max    <- ceiling(max(chelsa_vals))
temp_values <- make_scale_values(temp_min, temp_max)

all_diff_vals <- c(
  values(RACMO_MPI_hist_diff_10km, na.rm = TRUE),
  values(RACMO_CESM2_hist_diff_10km, na.rm = TRUE),
  values(RACMO_ERA5_hist_diff_10km, na.rm = TRUE)
)
diff_min_val <- floor(min(all_diff_vals, na.rm = TRUE))
diff_max_val <- ceiling(max(all_diff_vals, na.rm = TRUE))
diff_breaks  <- pretty(c(diff_min_val, diff_max_val), n = 6)
diff_values  <- make_scale_values(diff_min_val, diff_max_val)

message(
  "Historical bias range: ", round(min(all_diff_vals, na.rm = TRUE), 2),
  " to ", round(max(all_diff_vals, na.rm = TRUE), 2),
  "  ->  scale limits: ", diff_min_val, " to ", diff_max_val, "°C"
)

# --- Build panels ---

p_chelsa_hist <- make_temp_panel(
  CHELSA_hist_10km, "a)   Annual mean CHELSA",
  temp_min, temp_max, temp_values, coast
)

hist_models <- list(
  list(name = "RACMO-MPI-ESM1", rast = racmo_mpi_hist_10km,  diff = RACMO_MPI_hist_diff_10km),
  list(name = "RACMO-CESM2",    rast = racmo_cesm_hist_10km, diff = RACMO_CESM2_hist_diff_10km),
  list(name = "RACMO-ERA5",     rast = racmo_era5_hist_10km, diff = RACMO_ERA5_hist_diff_10km)
)

letter_idx <- 2
hist_model_panels <- list()

for (m in hist_models) {
  p_mod <- make_temp_panel(
    m$rast, paste0(letters[letter_idx], ")   Annual mean ", m$name),
    temp_min, temp_max, temp_values, coast
  )
  letter_idx <- letter_idx + 1
  
  p_dif <- make_diff_panel(
    m$diff, paste0(letters[letter_idx], ")   Difference in annual mean ", m$name),
    diff_min_val, diff_max_val, diff_breaks, diff_values, coast
  )
  letter_idx <- letter_idx + 1
  
  hist_model_panels <- c(hist_model_panels, list(p_mod, p_dif))
}

assemble_and_save(
  p_chelsa_hist, hist_model_panels,
  file.path(outpath, "plot_annual_temp_HISTORICAL_comparison_CHELSA_vs_RACMO.png"),
  height = 13.5
)

# ============================================================================
# 2. FUTURE (2071-2100 CHELSA multi-model mean vs 2071-2100 RACMO)
# ============================================================================

# --- Load CHELSA future members and build the ensemble mean ---

CHELSA_future1 <- rast(here("Data/CHELSA/Validation/Mean_Annual_Temperature_FUTURE_gfdl-esm4_2071_2100_ICEFREE.tif"))
CHELSA_future2 <- rast(here("Data/CHELSA/Validation/Mean_Annual_Temperature_FUTURE_ipsl-cm6a-lr_2071_2100_ICEFREE.tif"))
CHELSA_future3 <- rast(here("Data/CHELSA/Validation/Mean_Annual_Temperature_FUTURE_mpi-esm1-2-hr_2071_2100_ICEFREE.tif"))
CHELSA_future4 <- rast(here("Data/CHELSA/Validation/Mean_Annual_Temperature_FUTURE_mri-esm2-0_2071_2100_ICEFREE.tif"))
CHELSA_future5 <- rast(here("Data/CHELSA/Validation/Mean_Annual_Temperature_FUTURE_ukesm1-0-ll_2071_2100_ICEFREE.tif"))

CHELSA_future <- c(CHELSA_future1, CHELSA_future2, CHELSA_future3,
                   CHELSA_future4, CHELSA_future5)
CHELSA_future <- app(CHELSA_future, mean, na.rm = TRUE)

# --- Load RACMO future ---

racmo_mpi_future <- rast(
  here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_MPI_ESM1/Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif")
)

racmo_cesm_future <- rast(
  here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_CESM2/Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif")
)

# TEMPORARILY REMOVED -- RACMO_ERA5 future not finished yet. Restore by
# uncommenting this block and the RACMO-ERA5 lines below.
# racmo_era5_future <- rast(
#   here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_ERA5/Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif")
# )

# --- Bias ---

RACMO_MPI_future_diff   <- racmo_mpi_future  - align_to_model(CHELSA_future, racmo_mpi_future)
RACMO_CESM2_future_diff <- racmo_cesm_future - align_to_model(CHELSA_future, racmo_cesm_future)
# RACMO_ERA5_future_diff  <- racmo_era5_future - align_to_model(CHELSA_future, racmo_era5_future)

# --- Aggregate to 10 km ---

CHELSA_future_10km           <- aggregate(CHELSA_future, fact = 10, fun = "mean", na.rm = TRUE)
racmo_mpi_future_10km        <- aggregate(racmo_mpi_future, fact = 10, fun = "mean", na.rm = TRUE)
racmo_cesm_future_10km       <- aggregate(racmo_cesm_future, fact = 10, fun = "mean", na.rm = TRUE)
# racmo_era5_future_10km       <- aggregate(racmo_era5_future, fact = 10, fun = "mean", na.rm = TRUE)
RACMO_MPI_future_diff_10km   <- aggregate(RACMO_MPI_future_diff, fact = 10, fun = "mean", na.rm = TRUE)
RACMO_CESM2_future_diff_10km <- aggregate(RACMO_CESM2_future_diff, fact = 10, fun = "mean", na.rm = TRUE)
# RACMO_ERA5_future_diff_10km  <- aggregate(RACMO_ERA5_future_diff, fact = 10, fun = "mean", na.rm = TRUE)

# --- Colour scales ---

chelsa_vals <- values(CHELSA_future_10km, na.rm = TRUE)
temp_min    <- floor(min(chelsa_vals))
temp_max    <- ceiling(max(chelsa_vals))
temp_values <- make_scale_values(temp_min, temp_max)

all_diff_vals <- c(
  values(RACMO_MPI_future_diff_10km, na.rm = TRUE),
  values(RACMO_CESM2_future_diff_10km, na.rm = TRUE)
  # , values(RACMO_ERA5_future_diff_10km, na.rm = TRUE)
)
diff_min_val <- floor(min(all_diff_vals, na.rm = TRUE))
diff_max_val <- ceiling(max(all_diff_vals, na.rm = TRUE))
diff_breaks  <- pretty(c(diff_min_val, diff_max_val), n = 6)
diff_values  <- make_scale_values(diff_min_val, diff_max_val)

message(
  "Future bias range: ", round(min(all_diff_vals, na.rm = TRUE), 2),
  " to ", round(max(all_diff_vals, na.rm = TRUE), 2),
  "  ->  scale limits: ", diff_min_val, " to ", diff_max_val, "°C"
)

# --- Build panels ---

p_chelsa_future <- make_temp_panel(
  CHELSA_future_10km, "a)   Annual mean CHELSA",
  temp_min, temp_max, temp_values, coast
)

future_models <- list(
  list(name = "RACMO-MPI-ESM1", rast = racmo_mpi_future_10km,  diff = RACMO_MPI_future_diff_10km),
  list(name = "RACMO-CESM2",    rast = racmo_cesm_future_10km, diff = RACMO_CESM2_future_diff_10km)
  # , list(name = "RACMO-ERA5", rast = racmo_era5_future_10km, diff = RACMO_ERA5_future_diff_10km)
)

letter_idx <- 2
future_model_panels <- list()

for (m in future_models) {
  p_mod <- make_temp_panel(
    m$rast, paste0(letters[letter_idx], ")   Annual mean ", m$name),
    temp_min, temp_max, temp_values, coast
  )
  letter_idx <- letter_idx + 1
  
  p_dif <- make_diff_panel(
    m$diff, paste0(letters[letter_idx], ")   Difference in annual mean ", m$name),
    diff_min_val, diff_max_val, diff_breaks, diff_values, coast
  )
  letter_idx <- letter_idx + 1
  
  future_model_panels <- c(future_model_panels, list(p_mod, p_dif))
}

assemble_and_save(
  p_chelsa_future, future_model_panels,
  file.path(outpath, "plot_annual_temp_FUTURE_comparison_CHELSA_vs_RACMO.png"),
  height = 11.69
)



# OTHER VERSION CHELSA VS. RACMO ------------------------------------------

# CHELSA vs RACMO
# Historical + Future annual temperature comparison
# Structured to match the AntAirICE vs RACMO evaluation script

library(terra)
library(ggplot2)
library(dplyr)
library(here)
library(viridis)


# OUTPATH -----------------------------------------------------------------

outpath <- here("Plots/Evaluation_CHELSA")

if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}


# ============================================================
# HELPER: resample reference onto model's grid if needed
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


# -------------------------------------------------------------------------
# HISTORICAL (1981-2010 CHELSA vs 1985-2010 RACMO)
# -------------------------------------------------------------------------

CHELSA_hist <- rast(
  here(
    "Data/CHELSA/Mean_Annual_Temperature_HISTORICAL_1981_2010_ICEFREE.tif"
  )
)


# --- RACMO historical ---

racmo_mpi_hist <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_MPI_ESM1/comparison/Mean_Annual_Temperature_HISTORICAL_1985_2010_ICEFREE.tif"))

racmo_cesm_hist <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_CESM2/comparison/Mean_Annual_Temperature_HISTORICAL_1985_2010_ICEFREE.tif"))

racmo_era5_hist <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_ERA5/comparison/Mean_Annual_Temperature_HISTORICAL_1985_2010_ICEFREE.tif"))


# -------------------------------------------------------------------------
# FUTURE
# -------------------------------------------------------------------------

# CHELSA future ensemble / multi-model mean
CHELSA_future1 <- rast(
  here(
    "Data/CHELSA/Validation/",
    "Mean_Annual_Temperature_FUTURE_gfdl-esm4_2071_2100_ICEFREE.tif"
  )
)

CHELSA_future2 <- rast(
  here(
    "Data/CHELSA/Validation/",
    "Mean_Annual_Temperature_FUTURE_ipsl-cm6a-lr_2071_2100_ICEFREE.tif"
  )
)

CHELSA_future3 <- rast(
  here(
    "Data/CHELSA/Validation/",
    "Mean_Annual_Temperature_FUTURE_mpi-esm1-2-hr_2071_2100_ICEFREE.tif"
  )
)

CHELSA_future4 <- rast(
  here(
    "Data/CHELSA/Validation/",
    "Mean_Annual_Temperature_FUTURE_mri-esm2-0_2071_2100_ICEFREE.tif"
  )
)

CHELSA_future5 <- rast(
  here(
    "Data/CHELSA/Validation/",
    "Mean_Annual_Temperature_FUTURE_ukesm1-0-ll_2071_2100_ICEFREE.tif"
  )
)


CHELSA_future <- c(
  CHELSA_future1,
  CHELSA_future2,
  CHELSA_future3,
  CHELSA_future4,
  CHELSA_future5
)

CHELSA_future <- app(
  CHELSA_future,
  mean,
  na.rm = TRUE
)


# --- RACMO future ---

racmo_mpi_future <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_MPI_ESM1/Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif"))

racmo_cesm_future <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_CESM2/Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif"))

# TEMPORARILY REMOVED -- RACMO_ERA5 future not finished yet.
# racmo_era5_future <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded/RACMO_ERA5/Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif"))


# ============================================================
# MODEL LABELS
# ============================================================

model_labels <- c(
  "RACMO-MPI-ESM1",
  "RACMO-CESM2",
  "RACMO-ERA5"
)


# ============================================================
# HISTORICAL: SCATTER-DENSITY
# ============================================================

hist_pairs <- bind_rows(
  extract_pairs(
    CHELSA_hist,
    racmo_mpi_hist,
    "RACMO-MPI-ESM1"
  ),
  
  extract_pairs(
    CHELSA_hist,
    racmo_cesm_hist,
    "RACMO-CESM2"
  ),
  
  extract_pairs(
    CHELSA_hist,
    racmo_era5_hist,
    "RACMO-ERA5"
  )
)

hist_pairs$model <- factor(
  hist_pairs$model,
  levels = model_labels
)


# ============================================================
# FUTURE: SCATTER-DENSITY
# ============================================================

future_pairs <- bind_rows(
  extract_pairs(
    CHELSA_future,
    racmo_mpi_future,
    "RACMO-MPI-ESM1"
  ),
  
  extract_pairs(
    CHELSA_future,
    racmo_cesm_future,
    "RACMO-CESM2"
  )
  
  # TEMPORARILY REMOVED -- RACMO_ERA5 future not finished yet.
  # , extract_pairs(
  #   CHELSA_future,
  #   racmo_era5_future,
  #   "RACMO-ERA5"
  # )
)

future_pairs$model <- factor(
  future_pairs$model,
  levels = model_labels
)


# ============================================================
# SHARED AXIS LIMITS
# ============================================================

all_vals <- c(
  hist_pairs$x,
  hist_pairs$y,
  future_pairs$x,
  future_pairs$y
)

ax_lim <- range(
  all_vals,
  na.rm = TRUE
)


# ============================================================
# SCATTER-DENSITY PLOT FUNCTION
# ============================================================

scatter_density_row <- function(
    df,
    row_title,
    ref_period
) {
  
  ggplot(
    df,
    aes(x = x, y = y)
  ) +
    
    stat_bin2d(
      bins = 150,
      aes(fill = after_stat(count))
    ) +
    
    geom_abline(
      slope = 1,
      intercept = 0,
      colour = "red",
      linewidth = 0.5
    ) +
    
    scale_fill_viridis_c(
      option   = "plasma",
      name     = "Number of\ngrid cells",
      trans    = "sqrt",
      limits   = c(1, NA),
      na.value = NA
    ) +
    
    scale_x_continuous(
      limits = ax_lim,
      breaks = seq(-30, 10, 10)
    ) +
    
    scale_y_continuous(
      limits = ax_lim,
      breaks = seq(-30, 10, 10)
    ) +
    
    coord_fixed() +
    
    facet_wrap(
      ~ model,
      nrow = 1,
      strip.position = "top"
    ) +
    
    labs(
      x = paste0(
        "Temperature (CHELSA) [°C]"
      ),
      y = "Temperature (RACMO) [°C]",
      title = row_title
    ) +
    
    theme_classic(
      base_size = 9
    ) +
    
    theme(
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = 9
      ),
      legend.position = "right",
      legend.key.height = unit(
        1.8,
        "cm"
      ),
      legend.key.width = unit(
        0.35,
        "cm"
      ),
      panel.spacing = unit(
        0.4,
        "cm"
      ),
      plot.title = element_text(
        face = "bold",
        size = 10,
        hjust = 0.5
      ),
      axis.title.y = element_text(
        size = 8
      ),
      axis.title.x = element_text(
        size = 8
      )
    )
}


p_hist <- scatter_density_row(
  hist_pairs,
  "Historical: 1985-2010 RACMO",
  "Historical"
)

p_future <- scatter_density_row(
  future_pairs,
  "Future: 2071-2100 CHELSA",
  "Future"
)


plot_scatter <- p_hist / p_future +
  plot_annotation(
    title = "CHELSA vs RACMO temperatures",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 11,
        hjust = 0.5
      )
    )
  )


# ============================================================
# SAVE SCATTER-DENSITY PLOT
# ============================================================

ggsave(
  file.path(
    outpath,
    "CHELSA_vs_RACMO_scatter_density.png"
  ),
  plot_scatter,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

message("Scatter-density plot saved.")


# ============================================================
# BIAS RASTERS
# ============================================================

# Bias = RACMO - CHELSA

bias_list <- list(
  
  # Historical
  list(
    r = racmo_mpi_hist -
      align_to_model(
        CHELSA_hist,
        racmo_mpi_hist
      ),
    label = "RACMO-MPI-ESM1",
    period = "Historical"
  ),
  
  list(
    r = racmo_cesm_hist -
      align_to_model(
        CHELSA_hist,
        racmo_cesm_hist
      ),
    label = "RACMO-CESM2",
    period = "Historical"
  ),
  
  list(
    r = racmo_era5_hist -
      align_to_model(
        CHELSA_hist,
        racmo_era5_hist
      ),
    label = "RACMO-ERA5",
    period = "Historical"
  ),
  
  # Future
  list(
    r = racmo_mpi_future -
      align_to_model(
        CHELSA_future,
        racmo_mpi_future
      ),
    label = "RACMO-MPI-ESM1",
    period = "Future"
  ),
  
  list(
    r = racmo_cesm_future -
      align_to_model(
        CHELSA_future,
        racmo_cesm_future
      ),
    label = "RACMO-CESM2",
    period = "Future"
  )
  
  # TEMPORARILY REMOVED -- RACMO_ERA5 future not finished yet.
  # , list(
  #   r = racmo_era5_future -
  #     align_to_model(
  #       CHELSA_future,
  #       racmo_era5_future
  #     ),
  #   label = "RACMO-ERA5",
  #   period = "Future"
  # )
)


# ============================================================
# CONVERT BIAS RASTERS TO DATAFRAME
# ============================================================

bias_df <- bind_rows(
  lapply(
    bias_list,
    function(b) {
      
      df <- as.data.frame(
        b$r,
        xy = TRUE,
        na.rm = TRUE
      )
      
      names(df)[3] <- "bias"
      
      df$model  <- b$label
      df$period <- b$period
      
      df
    }
  )
)


bias_df$model <- factor(
  bias_df$model,
  levels = model_labels
)

bias_df$period <- factor(
  bias_df$period,
  levels = c(
    "Historical",
    "Future"
  )
)


# ============================================================
# BIAS COLOUR SCALE
# ============================================================

bias_lim <- 10


# ============================================================
# BIAS MAP
# ============================================================

plot_bias <- ggplot(
  bias_df,
  aes(
    x = x,
    y = y,
    fill = bias
  )
) +
  
  geom_tile() +
  
  scale_fill_distiller(
    palette = "RdBu",
    limits = c(
      -bias_lim,
      bias_lim
    ),
    name = "Bias (°C)\n(RACMO - CHELSA)",
    direction = -1,
    na.value = NA,
    oob = scales::squish
  ) +
  
  coord_fixed() +
  
  facet_grid(
    period ~ model
  ) +
  
  labs(
    title = "Temperature bias: RACMO - CHELSA",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_classic(
    base_size = 9
  ) +
  
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    legend.position = "right",
    legend.key.height = unit(
      2,
      "cm"
    ),
    legend.key.width = unit(
      0.4,
      "cm"
    ),
    panel.spacing = unit(
      0.3,
      "cm"
    ),
    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),
    axis.text = element_text(
      size = 7
    ),
    axis.title = element_text(
      size = 8
    ),
    panel.border = element_rect(
      colour = "grey70",
      fill = NA,
      linewidth = 0.3
    )
  )


# ============================================================
# SAVE BIAS MAP
# ============================================================

ggsave(
  file.path(
    outpath,
    "CHELSA_vs_RACMO_bias_maps.png"
  ),
  plot_bias,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

message("Bias maps saved.")