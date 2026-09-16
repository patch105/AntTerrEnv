# ==============================================================================
# Covariate Summary Maps -- per-RCM GCM comparison + Historical driver comparison
# Historical / Future / Future-Historical Difference
# Rasters are read as-is -- input data is already at ~10km resolution, so no
# aggregation/resampling happens in this script.
# ==============================================================================
#
# For each RCM group (HCLIM, RACMO -- run one, the other, or both in the same
# call via `run_rcms` below) this produces, per variable:
#
#  (1) A 6-panel GCM-storyline figure, laid out like:
#        (a) Hist  -- MPI-ESM1     (b) Hist  -- CESM2
#        (c) Future -- MPI-ESM1    (d) Future -- CESM2
#        (e) Diff  -- MPI-ESM1     (f) Diff  -- CESM2
#      Colour scale is shared across the two model columns WITHIN a row (so
#      both Historical panels share a scale, both Future panels share a
#      scale, and both Diff panels share a diverging scale), but each row
#      gets its own scale/range.
#
#  (2) A 3-panel HISTORICAL-only driver comparison figure:
#        (a) MPI-ESM1   (b) CESM2   (c) ERA5
#      all sharing one colour scale, since ERA5 only ever has a HISTORICAL
#      period (same as in the Step 2/summary scripts) and is the natural
#      "how well does the driver-comparison line up with the evaluation run"
#      check.
#
# The full set of headline variables plotted is pulled straight from the
# Step 2/PolarRes26 summary script's variable -> filename mapping (direct
# annual/seasonal rasters, the two "monthly mean" families TasMin/TasMax,
# and sea ice), so this script and the summary-table script stay in sync.
#
# Reads the raw per-model rasters directly rather than building a tidy data
# frame, since the goal here is spatial maps, not aggregate tables. Rasters
# are plotted at their native resolution (no aggregation/resampling step).
# ==============================================================================

library(terra)
library(here)
library(ggplot2)
library(patchwork)

# ---- 1. Configuration -------------------------------------------------------

# Which RCM group(s) to run in this call. Set to "HCLIM", "RACMO", or
# c("HCLIM", "RACMO") to do both in one run -- each gets its own output
# subfolder so the figures never overwrite one another.
run_rcms <- c("HCLIM", "RACMO")

# GCM-driven models + the matching ERA5 evaluation run for each RCM. Add
# further RCMs here (e.g. MetUM) following the same shape if needed --
# MetUM only has an ERA5 run in this dataset, so it wouldn't get a 6-panel
# GCM-storyline figure, only the historical comparison (with just 1 panel).
rcm_groups <- list(
  HCLIM = list(
    GCM  = c(MPI_ESM1 = "HCLIM_MPI_ESM1", CESM2 = "HCLIM_CESM2"),
    ERA5 = "HCLIM_ERA5"
  ),
  RACMO = list(
    GCM  = c(MPI_ESM1 = "RACMO_MPI_ESM1", CESM2 = "RACMO_CESM2"),
    ERA5 = "RACMO_ERA5"
  )
)

years_hist   <- seq(1995, 2014, by = 1)
years_future <- seq(2081, 2100, by = 1)

periods <- list(
  HISTORICAL = paste(min(years_hist),   max(years_hist),   sep = "_"),
  FUTURE     = paste(min(years_future), max(years_future), sep = "_")
)

input_base    <- here("Data/Environmental_predictors/PolarRes26/Regridded")
outpath_base  <- here("Plots/Covariate_Summary_Maps")

# Sea ice: plain (non-buffered) concentration, averaged across its Oct-Feb
# season -- the only months Step 2 computes for it (matches the summary script).
sea_ice_months <- c("October", "November", "December", "January", "February")

# ---- 2. Headline variables to plot ------------------------------------------
# Pulled from the Step 2/PolarRes26 summary script's variable -> filename
# mapping, so the two scripts stay in sync. Each entry has:
#   type    : "direct" (one annual/seasonal raster per period), "monthly_mean"
#             (12 monthly rasters averaged here), or "sea_ice" (Oct-Feb rasters
#             averaged here)
#   pattern : sprintf pattern -- direct/sea_ice take (period_name, range_label)
#             [sea_ice also takes the month name as its first %s]; monthly_mean
#             takes (month_name, period_name, range_label)
#   label   : display label/unit for the legend title / figure title

direct_variable_defs <- list(
  AnnualTemp        = list(pattern = "Mean_Annual_Temperature_%s_%s_ICEFREE.tif",        label = "Annual Temp (\u00b0C)"),
  SummerTemp        = list(pattern = "Mean_Summer_Temperature_%s_%s_ICEFREE.tif",        label = "Summer Temp (\u00b0C)"),
  WinterTemp        = list(pattern = "Mean_Winter_Temperature_%s_%s_ICEFREE.tif",        label = "Winter Temp (\u00b0C)"),
  DegreeDaysMinus5  = list(pattern = "Mean_Total_Annual_Degree_Days-5_%s_%s_ICEFREE.tif", label = "Degree Days (-5\u00b0C)"),
  DegreeDays0       = list(pattern = "Mean_Total_Annual_Degree_Days0_%s_%s_ICEFREE.tif",  label = "Degree Days (0\u00b0C)"),
  WindSpeed         = list(pattern = "Mean_Annual_Wind_Speed_%s_%s_ICEFREE.tif",          label = "Wind Speed (m/s)"),
  TotalAnnualPrecip = list(pattern = "Total_Annual_Precipitation_%s_%s_ICEFREE.tif",      label = "Total Annual Precip (mm)"),
  TotalSummerPrecip = list(pattern = "Mean_Total_Summer_Precipitation_%s_%s_ICEFREE.tif", label = "Total Summer Precip (mm)"),
  MeanAnnualPrecip  = list(pattern = "Mean_Annual_Precipitation_%s_%s_ICEFREE.tif",       label = "Mean Annual Precip (mm)"),
  MeanSummerPrecip  = list(pattern = "Mean_Summer_Precipitation_%s_%s_ICEFREE.tif",       label = "Mean Summer Precip (mm)"),
  SolarRad          = list(pattern = "Mean_Annual_Solar_Radiation_%s_%s_ICEFREE.tif",     label = "Solar Radiation"),
  MeanMelt          = list(pattern = "Mean_Annual_Melt_%s_%s_ICEFREE.tif",                label = "Mean Annual Melt"),
  TotalMelt         = list(pattern = "Mean_Total_Annual_Melt_%s_%s_ICEFREE.tif",          label = "Total Annual Melt"),
  SnowCover         = list(pattern = "Mean_Annual_Snow_Cover_%s_%s_ICEFREE.tif",          label = "Snow Cover"),
  SummerRelHumidity = list(pattern = "Mean_Summer_Relative_Humidity_%s_%s_ICEFREE.tif",   label = "Summer Rel. Humidity (%)"),
  WinterRelHumidity = list(pattern = "Mean_Winter_Relative_Humidity_%s_%s_ICEFREE.tif",   label = "Winter Rel. Humidity (%)"),
  RelHumidity       = list(pattern = "Mean_Annual_Relative_Humidity_%s_%s_ICEFREE.tif",   label = "Relative Humidity (%)"),
  VPD               = list(pattern = "Mean_Annual_VPD_%s_%s_ICEFREE.tif",                 label = "VPD")
)

# BIO1-BIO19 all follow "BIO<n>_%s_%s_ICEFREE.tif" -- generated rather than
# typed out 19 times.
bio_variable_defs <- setNames(
  lapply(1:19, function(i) list(pattern = sprintf("BIO%d_%%s_%%s_ICEFREE.tif", i),
                                label = paste0("BIO", i))),
  paste0("BIO", 1:19)
)

# Families where Step 2 only saved 12 monthly files -- averaged into a mean
# annual raster here.
monthly_mean_variable_defs <- list(
  TasMin = list(pattern = "Climatological_Monthly_Mean_TasMin_%s_%s_%s_ICEFREE.tif", label = "Min Temp -- TasMin (\u00b0C)"),
  TasMax = list(pattern = "Climatological_Monthly_Mean_TasMax_%s_%s_%s_ICEFREE.tif", label = "Max Temp -- TasMax (\u00b0C)")
)

sea_ice_variable_defs <- list(
  SeaIceConc = list(pattern = "Mean_%s_Sea_Ice_Concentration_%s_%s_ICEFREE.tif", label = "Sea Ice Concentration (%)")
)

# Tag every entry with its type and fold them all into one list to loop over.
plot_variables <- c(
  lapply(direct_variable_defs,        function(v) c(v, type = "direct")),
  lapply(bio_variable_defs,           function(v) c(v, type = "direct")),
  lapply(monthly_mean_variable_defs,  function(v) c(v, type = "monthly_mean")),
  lapply(sea_ice_variable_defs,       function(v) c(v, type = "sea_ice"))
)

# ---- 3. Helpers --------------------------------------------------------------

# Read one model/period/variable raster, dispatching on the variable's type.
# For monthly_mean / sea_ice types this already returns the across-month mean
# (i.e. the annual/seasonal aggregate), matching what the summary script does.
read_variable_raster <- function(model, var_def, period_name) {
  model_dir   <- file.path(input_base, model)
  range_label <- periods[[period_name]]
  
  if (var_def$type == "direct") {
    path <- file.path(model_dir, sprintf(var_def$pattern, period_name, range_label))
    if (!file.exists(path)) { warning("Missing: ", path); return(NULL) }
    return(rast(path))
  }
  
  if (var_def$type == "monthly_mean") {
    paths <- file.path(model_dir, sprintf(var_def$pattern, month.name, period_name, range_label))
    if (any(!file.exists(paths))) { warning("Missing monthly file(s): ", var_def$pattern, " (", model, ", ", period_name, ")"); return(NULL) }
    return(app(rast(paths), mean, na.rm = TRUE))
  }
  
  if (var_def$type == "sea_ice") {
    paths <- file.path(model_dir, sprintf(var_def$pattern, sea_ice_months, period_name, range_label))
    if (any(!file.exists(paths))) { warning("Missing sea-ice monthly file(s) for ", model, ", ", period_name); return(NULL) }
    return(app(rast(paths), mean, na.rm = TRUE))
  }
  
  stop("Unknown variable type: ", var_def$type)
}

# SpatRaster -> tidy x/y/value data frame (drops NA cells, since the domain
# is already ice-free-masked upstream).
raster_to_df <- function(r, value_name = "value") {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- value_name
  df
}

# One filled map panel. `limits` is the shared min/max for its row/figure;
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

# ---- 4. Build + save figures, per RCM group ---------------------------------

for (rcm_name in run_rcms) {
  
  group      <- rcm_groups[[rcm_name]]
  models     <- group$GCM                       # c(MPI_ESM1 = "..._MPI_ESM1", CESM2 = "..._CESM2")
  era5_model <- group$ERA5                      # "..._ERA5"
  
  model_labels <- setNames(gsub("_", "-", models), models)
  era5_label   <- gsub("_", "-", era5_model)
  
  outpath <- file.path(outpath_base, rcm_name)
  dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
  
  message("\n##### RCM group: ", rcm_name, " #####")
  
  era5_dir_exists <- dir.exists(file.path(input_base, era5_model))
  if (!era5_dir_exists) {
    message("  (note: ", era5_model, " directory not found -- historical driver-comparison figures will be skipped)")
  }
  
  for (var_name in names(plot_variables)) {
    
    var_def   <- plot_variables[[var_name]]
    var_label <- var_def$label
    message("=== ", var_name, " ===")
    
    # ---- 4a. 6-panel GCM-storyline figure (Hist / Future / Diff) ----------
    
    rasters <- list()
    for (model in models) {
      for (period_name in names(periods)) {
        r <- read_variable_raster(model, var_def, period_name)
        if (is.null(r)) next
        rasters[[paste(model, period_name, sep = "_")]] <- r
      }
    }
    
    if (length(rasters) < length(models) * length(periods)) {
      message("  -- skipping ", var_name, " GCM-storyline map, missing one or more required rasters")
    } else {
      
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
    
    # ---- 4b. 3-panel HISTORICAL driver comparison (MPI-ESM1, CESM2, ERA5) --
    
    if (!era5_dir_exists) next
    
    hist_rasters <- list()
    for (model in c(models, era5_model)) {
      r <- read_variable_raster(model, var_def, "HISTORICAL")
      if (is.null(r)) next
      hist_rasters[[model]] <- r
    }
    
    if (length(hist_rasters) < length(models) + 1) {
      message("  -- skipping ", var_name, " historical driver-comparison map, missing one or more required rasters")
      next
    }
    
    hist_comp_vals   <- unlist(lapply(hist_rasters, values, na.rm = TRUE))
    hist_comp_limits <- range(hist_comp_vals, na.rm = TRUE)
    
    panel_order  <- c(models, era5_model)
    panel_labels <- c(model_labels, setNames(era5_label, era5_model))
    
    hist_comp_panels <- lapply(panel_order, function(m) {
      df <- raster_to_df(hist_rasters[[m]])
      make_map_panel(df, panel_labels[[m]], hist_comp_limits)
    })
    
    fig_hist_comp <- (hist_comp_panels[[1]] | hist_comp_panels[[2]] | hist_comp_panels[[3]]) +
      plot_annotation(
        title = var_label,
        subtitle = paste0("Historical (", periods$HISTORICAL, ") -- driver comparison"),
        tag_levels = "a"
      )
    
    out_file_hist_comp <- file.path(outpath, paste0(var_name, "_historical_driver_comparison_map.png"))
    ggsave(out_file_hist_comp, fig_hist_comp, width = 14, height = 5, dpi = 300, bg = "white")
    message("  -- saved ", out_file_hist_comp)
  }
  
  message("\nDone with ", rcm_name, ". Figures written to: ", outpath)
}

message("\nAll requested RCM groups complete.")

