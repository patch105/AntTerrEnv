# ==============================================================================
# Covariate Summary Maps -- HCLIM_MPI_ESM1 vs HCLIM_CESM2
# Historical / Future / Future-Historical Difference, aggregated to ~10km
# ==============================================================================
#
# Produces one 6-panel figure per variable, laid out like:
#   (a) Hist  -- MPI-ESM1     (b) Hist  -- CESM2
#   (c) Future -- MPI-ESM1    (d) Future -- CESM2
#   (e) Diff  -- MPI-ESM1     (f) Diff  -- CESM2
#
# Colour scale is shared across the two model columns WITHIN a row (so the
# two Historical panels use the same scale, the two Future panels use the
# same scale, and the two Diff panels use the same diverging scale), but
# each row gets its own scale/range -- matching the reference figure.
#
# Reuses the exact filename patterns and period/range-label conventions
# from the Step 2/summary script, but reads the raw per-model rasters
# directly instead of building a Peninsula/Continent tidy data frame, since
# the goal here is spatial maps, not aggregate tables.
# ==============================================================================

library(terra)
library(here)
library(ggplot2)
library(patchwork)

# ---- 1. Configuration -------------------------------------------------------

models <- c("HCLIM_MPI_ESM1", "HCLIM_CESM2")
#models <- c("RACMO_MPI_ESM1", "RACMO_CESM2")

model_labels <- c(HCLIM_MPI_ESM1 = "HCLIM-MPI-ESM1", HCLIM_CESM2 = "HCLIM-CESM2")
# model_labels <- c(RACMO_MPI_ESM1 = "RACMO-MPI-ESM1", RACMO_CESM2 = "RACMO-CESM2")

outpath    <- here("Plots/Covariate_Summary_Maps/HCLIM")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

years_hist   <- seq(1995, 2014, by = 1)
years_future <- seq(2081, 2100, by = 1)

periods <- list(
  HISTORICAL = paste(min(years_hist),   max(years_hist),   sep = "_"),
  FUTURE     = paste(min(years_future), max(years_future), sep = "_")
)

input_base <- here("Data/Environmental_predictors/PolarRes26/Regridded")


# Target aggregation resolution, in the CRS units of the rasters (assumed
# metres -- these are polar-stereographic grids). Adjust if your rasters use
# a different CRS/unit.
target_res_m <- 10000  # 10 km

# Headline variables to plot, with their Step-2 filename patterns (%s %s ->
# period_name, range_label) and a display label/unit for the legend title.
plot_variables <- list(
  AnnualTemp        = list(pattern = "Mean_Annual_Temperature_%s_%s_ICEFREE.tif",       label = "Annual Temp (\u00b0C)"),
  SummerTemp        = list(pattern = "Mean_Summer_Temperature_%s_%s_ICEFREE.tif",       label = "Summer Temp (\u00b0C)"),
  WinterTemp        = list(pattern = "Mean_Winter_Temperature_%s_%s_ICEFREE.tif",       label = "Winter Temp (\u00b0C)"),
  DegreeDaysMinus5  = list(pattern = "Mean_Total_Annual_Degree_Days-5_%s_%s_ICEFREE.tif", label = "Degree Days (-5\u00b0C)"),
  DegreeDays0       = list(pattern = "Mean_Total_Annual_Degree_Days0_%s_%s_ICEFREE.tif",  label = "Degree Days (0\u00b0C)"),
  TotalAnnualPrecip = list(pattern = "Total_Annual_Precipitation_%s_%s_ICEFREE.tif",      label = "Total Annual Precip (mm)"),
  WindSpeed         = list(pattern = "Mean_Annual_Wind_Speed_%s_%s_ICEFREE.tif",          label = "Wind Speed (m/s)"),
  RelHumidity       = list(pattern = "Mean_Annual_Relative_Humidity_%s_%s_ICEFREE.tif",   label = "Relative Humidity (%)")
)


# ---- 2. Helpers --------------------------------------------------------------

# Read one model/period/variable raster, following the same path/pattern
# convention as the Step 2/summary script.
read_var_raster <- function(model, var_pattern, period_name) {
  range_label <- periods[[period_name]]
  path <- file.path(input_base, model, sprintf(var_pattern, period_name, range_label))
  if (!file.exists(path)) {
    warning("Missing: ", path)
    return(NULL)
  }
  rast(path)
}

# Aggregate a raster to ~10km resolution, ignoring NAs in the average.
aggregate_10km <- function(r, target_res = target_res_m) {
  current_res <- res(r)[1]
  fact <- round(target_res / current_res)
  if (fact < 1) fact <- 1
  aggregate(r, fact = fact, fun = mean, na.rm = TRUE)
}

# SpatRaster -> tidy x/y/value data frame (drops NA cells, since the domain
# is already ice-free-masked upstream).
raster_to_df <- function(r, value_name = "value") {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- value_name
  df
}

# One filled map panel. `limits` is the shared min/max for its row;
# `diverging` switches to a zero-centred diverging scale for the diff row.
make_map_panel <- function(df, title, limits, diverging = FALSE) {
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    coord_equal() +
    labs(title = title, fill = NULL, x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0, size = 11),
      legend.key.height = unit(0.8, "cm")
    )
  if (diverging) {
    lim <- max(abs(limits), na.rm = TRUE)
    p + scale_fill_gradient2(low = "#440154", mid = "#21908C", high = "#FDE725",
                             midpoint = 0, limits = c(-lim, lim))
  } else {
    p + scale_fill_viridis_c(limits = limits, na.value = NA)
  }
}

# ---- 3. Build + save one 6-panel figure per variable -------------------------

for (var_name in names(plot_variables)) {
  
  var_pattern <- plot_variables[[var_name]]$pattern
  var_label   <- plot_variables[[var_name]]$label
  message("=== ", var_name, " ===")
  
  # Read + aggregate to 10km for every model x period combo.
  rasters <- list()
  for (model in models) {
    for (period_name in names(periods)) {
      r <- read_var_raster(model, var_pattern, period_name)
      if (is.null(r)) next
      rasters[[paste(model, period_name, sep = "_")]] <- aggregate_10km(r)
    }
  }
  
  if (length(rasters) < length(models) * length(periods)) {
    message("  -- skipping ", var_name, ", missing one or more required rasters")
    next
  }
  
  # Future - Historical difference, per model.
  diffs <- list()
  for (model in models) {
    hist_r   <- rasters[[paste(model, "HISTORICAL", sep = "_")]]
    future_r <- rasters[[paste(model, "FUTURE", sep = "_")]]
    diffs[[model]] <- future_r - hist_r
  }
  
  # ---- Row-wise shared colour scales (pooled across both models) ----
  hist_vals   <- unlist(lapply(models, function(m) values(rasters[[paste(m, "HISTORICAL", sep = "_")]], na.rm = TRUE)))
  future_vals <- unlist(lapply(models, function(m) values(rasters[[paste(m, "FUTURE", sep = "_")]], na.rm = TRUE)))
  diff_vals   <- unlist(lapply(models, function(m) values(diffs[[m]], na.rm = TRUE)))
  
  hist_limits   <- range(hist_vals, na.rm = TRUE)
  future_limits <- range(future_vals, na.rm = TRUE)
  diff_limits   <- range(diff_vals, na.rm = TRUE)
  
  # ---- Build panels ----
  hist_panels <- lapply(models, function(m) {
    df <- raster_to_df(rasters[[paste(m, "HISTORICAL", sep = "_")]])
    make_map_panel(df, model_labels[[m]], hist_limits)
  })
  future_panels <- lapply(models, function(m) {
    df <- raster_to_df(rasters[[paste(m, "FUTURE", sep = "_")]])
    make_map_panel(df, model_labels[[m]], future_limits)
  })
  diff_panels <- lapply(models, function(m) {
    df <- raster_to_df(diffs[[m]])
    make_map_panel(df, model_labels[[m]], diff_limits, diverging = TRUE)
  })
  
  # ---- Assemble: 3 rows (Hist/Future/Diff) x 2 columns (models) ----
  fig <- (hist_panels[[1]]   | hist_panels[[2]]) /
    (future_panels[[1]] | future_panels[[2]]) /
    (diff_panels[[1]]   | diff_panels[[2]]) +
    plot_annotation(
      title = var_label,
      subtitle = paste0("Top: Historical (", periods$HISTORICAL, ")   ",
                        "Middle: Future (", periods$FUTURE, ")   ",
                        "Bottom: Future \u2212 Historical"),
      tag_levels = "a"
    )
  
  out_file <- file.path(outpath, paste0(var_name, "_summary_map.png"))
  ggsave(out_file, fig, width = 10, height = 13, dpi = 300, bg = "white")
  message("  -- saved ", out_file)
}

message("\nDone. Figures written to: ", outpath)
