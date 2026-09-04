# ==============================================================================
# PolarRes26 -- HISTORICAL / MID / FUTURE CLIMATE SUMMARIES
# ==============================================================================
# Summarises the headline climatology rasters produced by the Step 2
# summarising script, across ALL PolarRes26 models, for Peninsula vs
# Continent -- for HISTORICAL, MID and FUTURE periods, plus MID-HISTORICAL
# and FUTURE-HISTORICAL differences.
#
# One headline raster per variable family is used: whatever annual/seasonal
# aggregate Step 2 already saves for that family (e.g. Mean_Annual_
# Temperature). Three families -- TasMin, TasMax, and sea ice -- don't have
# a saved annual aggregate (Step 2 only saves their monthly/seasonal
# climatology files), so those are computed here as a mean across their
# monthly files instead. Sea ice uses the plain (non-buffered) concentration,
# averaged across its Oct-Feb season.
#
# Ice-free masking is NOT done here -- it's already been applied upstream by
# the Grid_adjust script, which regrids and masks everything onto the common
# grid and writes the "_ICEFREE" versions this script reads. The ice-free
# domain rasters are still loaded, but only as a sanity check that the grid
# actually lines up (no resampling/regridding happens here -- if it doesn't
# line up, this stops with an error rather than silently fixing it).
#
# ERA5-driven models (HCLIM_ERA5, RACMO_ERA5, MetUM_ERA5) only ever have a
# HISTORICAL period -- same as in the Step 2/2b scripts -- so they
# contribute a HISTORICAL row only, with no MID/FUTURE diff.
#
# Model names are RCM_Driver (e.g. HCLIM_CESM2, RACMO_ERA5, MetUM_ERA5).
# RCM and Driver are split out of the model name. "Storyline" (by-driver)
# tables only make sense for the two SSP-driven GCMs -- CESM2 and MPI_ESM1
# -- since ERA5 is an evaluation run, not a climate storyline, so it's
# excluded from those tables (but still included in the by-model tables).
#
# NOTE on diffs: MID_DIFF/FUTURE_DIFF are only computed for variables that
# actually have both a HISTORICAL raster and a MID/FUTURE raster. terra's
# `-` operator on SpatRasters is POSITIONAL (it pairs layers by index, not
# by name), so if a variable is missing from one period but present in the
# other, a naive `target - hist` silently misaligns every layer after that
# point and produces unnamed/garbage columns. See diff_common_layers().
# ==============================================================================

library(dplyr)
library(purrr)
library(terra)
library(here)
library(tidyr)
library(stringr)

# ---- 1. Configuration -----------------------------------------------------

models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5",
            "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5", "MetUM_ERA5")
# models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "RACMO_MPI_ESM1")

years_hist   <- seq(1995, 2014, by = 1)
years_mid    <- seq(2041, 2060, by = 1)
years_future <- seq(2081, 2100, by = 1)

periods <- list(
  HISTORICAL = paste(min(years_hist),   max(years_hist),   sep = "_"),
  MID        = paste(min(years_mid),    max(years_mid),    sep = "_"),
  FUTURE     = paste(min(years_future), max(years_future), sep = "_")
)

input_base <- here("Data/Environmental_predictors/PolarRes26/Regridded")
outpath    <- here("Data/Environmental_predictors/PolarRes26/Scenario_Summaries")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

# Peninsula / Continent boundary
pen_cont_boundary <- vect(here("Data/Peninsula_Continent_Boundary.shp"))

# RCM / Driver derived from each model name (e.g. "HCLIM_CESM2" -> RCM =
# "HCLIM", Driver = "CESM2"; "RACMO_MPI_ESM1" -> RCM = "RACMO",
# Driver = "MPI_ESM1").
model_meta <- tibble(Model = models) %>%
  separate(Model, into = c("RCM", "Driver"), sep = "_", extra = "merge", remove = FALSE)

# Only these two drivers represent an actual SSP storyline -- ERA5 is an
# evaluation run and is excluded from the by-driver ("storyline") tables.
storyline_drivers <- c("CESM2", "MPI_ESM1")

# ---- 2. Headline variable -> filename mapping ------------------------------
# All filenames carry the "_ICEFREE" suffix Grid_adjust writes.

direct_variables <- list(
  AnnualTemp       = "Mean_Annual_Temperature_%s_%s_ICEFREE.tif",
  SummerTemp       = "Mean_Summer_Temperature_%s_%s_ICEFREE.tif",
  WinterTemp       = "Mean_Winter_Temperature_%s_%s_ICEFREE.tif",
  DegreeDaysMinus5 = "Mean_Total_Annual_Degree_Days-5_%s_%s_ICEFREE.tif",
  DegreeDays0      = "Mean_Total_Annual_Degree_Days0_%s_%s_ICEFREE.tif",
  WindSpeed        = "Mean_Annual_Wind_Speed_%s_%s_ICEFREE.tif",
  TotalAnnualPrecip= "Total_Annual_Precipitation_%s_%s_ICEFREE.tif",
  TotalSummerPrecip= "Mean_Total_Summer_Precipitation_%s_%s_ICEFREE.tif",
  MeanAnnualPrecip = "Mean_Annual_Precipitation_%s_%s_ICEFREE.tif",
  MeanSummerPrecip = "Mean_Summer_Precipitation_%s_%s_ICEFREE.tif",
  SolarRad         = "Mean_Annual_Solar_Radiation_%s_%s_ICEFREE.tif",
  MeanMelt         = "Mean_Annual_Melt_%s_%s_ICEFREE.tif",
  TotalMelt        = "Mean_Total_Annual_Melt_%s_%s_ICEFREE.tif",
  SnowCover        = "Mean_Annual_Snow_Cover_%s_%s_ICEFREE.tif",
  SummerRelHumidity= "Mean_Summer_Relative_Humidity_%s_%s_ICEFREE.tif",
  WinterRelHumidity= "Mean_Winter_Relative_Humidity_%s_%s_ICEFREE.tif",
  RelHumidity      = "Mean_Annual_Relative_Humidity_%s_%s_ICEFREE.tif",
  VPD              = "Mean_Annual_VPD_%s_%s_ICEFREE.tif",
  BIO1             = "BIO1_%s_%s_ICEFREE.tif",
  BIO2             = "BIO2_%s_%s_ICEFREE.tif",
  BIO3             = "BIO3_%s_%s_ICEFREE.tif",
  BIO4             = "BIO4_%s_%s_ICEFREE.tif",
  BIO5             = "BIO5_%s_%s_ICEFREE.tif",
  BIO6             = "BIO6_%s_%s_ICEFREE.tif",
  BIO7             = "BIO7_%s_%s_ICEFREE.tif",
  BIO8             = "BIO8_%s_%s_ICEFREE.tif",
  BIO9             = "BIO9_%s_%s_ICEFREE.tif",
  BIO10            = "BIO10_%s_%s_ICEFREE.tif",
  BIO11            = "BIO11_%s_%s_ICEFREE.tif",
  BIO12            = "BIO12_%s_%s_ICEFREE.tif",
  BIO13            = "BIO13_%s_%s_ICEFREE.tif",
  BIO14            = "BIO14_%s_%s_ICEFREE.tif",
  BIO15            = "BIO15_%s_%s_ICEFREE.tif",
  BIO16            = "BIO16_%s_%s_ICEFREE.tif",
  BIO17            = "BIO17_%s_%s_ICEFREE.tif",
  BIO18            = "BIO18_%s_%s_ICEFREE.tif",
  BIO19            = "BIO19_%s_%s_ICEFREE.tif"
)


# Families where Step 2 only saved 12 monthly files -- averaged into a mean
# annual raster here.
monthly_mean_variables <- list(
  TasMin = "Climatological_Monthly_Mean_TasMin_%s_%s_%s_ICEFREE.tif",
  TasMax = "Climatological_Monthly_Mean_TasMax_%s_%s_%s_ICEFREE.tif"
)

# Sea ice: plain (non-buffered) concentration, averaged across its Oct-Feb
# season -- the only months Step 2 computes for it.
sea_ice_pattern <- "Mean_%s_Sea_Ice_Concentration_%s_%s_ICEFREE.tif"
sea_ice_months  <- c("October", "November", "December", "January", "February")

# ---- 3. Helpers -------------------------------------------------------------

read_direct <- function(model_dir, pattern, period_name, range_label) {
  path <- file.path(model_dir, sprintf(pattern, period_name, range_label))
  if (!file.exists(path)) { warning("Missing: ", path); return(NULL) }
  rast(path)
}

read_monthly_mean <- function(model_dir, pattern, period_name, range_label) {
  paths <- file.path(model_dir, sprintf(pattern, month.name, period_name, range_label))
  if (any(!file.exists(paths))) { warning("Missing monthly file(s): ", pattern); return(NULL) }
  app(rast(paths), mean, na.rm = TRUE)
}

read_sea_ice_mean <- function(model_dir, period_name, range_label) {
  paths <- file.path(model_dir, sprintf(sea_ice_pattern, sea_ice_months, period_name, range_label))
  if (any(!file.exists(paths))) { warning("Missing sea-ice monthly file(s) for ", period_name); return(NULL) }
  app(rast(paths), mean, na.rm = TRUE)
}

# Subtract two SpatRasters variable-by-variable, matching on LAYER NAME
# rather than position. terra's `-` on SpatRasters is positional, so if one
# period is missing a variable the other has (e.g. no FUTURE file for some
# covariate), a plain `target - hist` silently pairs up the wrong layers and
# produces an unnamed/garbage result. Here we only diff the variables present
# in BOTH periods and drop the rest, returning NULL if there's no overlap at
# all (so callers can skip that region/period pair entirely).
diff_common_layers <- function(target_r, hist_r) {
  if (is.null(target_r) || is.null(hist_r)) return(NULL)
  common <- intersect(names(target_r), names(hist_r))
  missing_from_target <- setdiff(names(hist_r), names(target_r))
  missing_from_hist    <- setdiff(names(target_r), names(hist_r))
  if (length(missing_from_target) > 0) {
    message("    (diff) skipping vars with no target period raster: ",
            paste(missing_from_target, collapse = ", "))
  }
  if (length(missing_from_hist) > 0) {
    message("    (diff) skipping vars with no HISTORICAL raster: ",
            paste(missing_from_hist, collapse = ", "))
  }
  if (length(common) == 0) return(NULL)
  target_r[[common]] - hist_r[[common]]
}

# Split into Peninsula / Continent SpatRasters (kept as rasters, not data
# frames, so HISTORICAL/MID/FUTURE can still be subtracted cell-by-cell
# before any NA rows get dropped).
split_regions <- function(stack) {
  peninsula <- mask(stack, pen_cont_boundary)
  peninsula <- crop(peninsula, ext(pen_cont_boundary))
  
  continent <- mask(stack, pen_cont_boundary, inverse = TRUE)
  
  list(peninsula = peninsula, continent = continent)
}

# ---- 4. Build a Peninsula/Continent tidy data frame for every model -------

model_dfs <- map(models, function(model) {
  
  model_dir <- file.path(input_base, model)
  if (!dir.exists(model_dir)) {
    message("!! Skipping ", model, " -- directory not found: ", model_dir)
    return(NULL)
  }
  message("=== ", model, " ===")
  
  is_era5 <- endsWith(model, "_ERA5")
  period_names_to_run <- if (is_era5) "HISTORICAL" else names(periods)
  
  period_stacks <- map(period_names_to_run, function(period_name) {
    range_label <- periods[[period_name]]
    message("  -- ", period_name)
    
    layers <- list()
    for (cov in names(direct_variables)) {
      r <- read_direct(model_dir, direct_variables[[cov]], period_name, range_label)
      if (!is.null(r)) layers[[cov]] <- r
    }
    for (cov in names(monthly_mean_variables)) {
      r <- read_monthly_mean(model_dir, monthly_mean_variables[[cov]], period_name, range_label)
      if (!is.null(r)) layers[[cov]] <- r
    }
    r_seaice <- read_sea_ice_mean(model_dir, period_name, range_label)
    if (!is.null(r_seaice)) layers[["SeaIceConc"]] <- r_seaice
    
    if (length(layers) == 0) {
      message("     no files found -- skipping ", period_name)
      return(NULL)
    }
    
    stack <- rast(unname(layers))
    names(stack) <- names(layers)
    
    split_regions(stack)
  })
  names(period_stacks) <- period_names_to_run
  period_stacks <- compact(period_stacks)
  if (length(period_stacks) == 0) return(NULL)
  
  # Differences vs HISTORICAL (only possible where both periods are present)
  # -- matched by variable NAME, not position, so a variable missing from
  # only one of the two periods is dropped from the diff instead of
  # silently misaligning every layer after it.
  diffs <- list()
  if (!is.null(period_stacks$HISTORICAL)) {
    for (target in c("MID", "FUTURE")) {
      if (!is.null(period_stacks[[target]])) {
        message("  -- ", target, "_DIFF vs HISTORICAL")
        pen_diff  <- diff_common_layers(period_stacks[[target]]$peninsula,
                                        period_stacks$HISTORICAL$peninsula)
        cont_diff <- diff_common_layers(period_stacks[[target]]$continent,
                                        period_stacks$HISTORICAL$continent)
        if (!is.null(pen_diff) || !is.null(cont_diff)) {
          diffs[[paste0(target, "_DIFF")]] <- list(
            peninsula = pen_diff,
            continent = cont_diff
          )
        }
      }
    }
  }
  
  all_stacks <- c(period_stacks, diffs)
  
  # Turn every Peninsula/Continent SpatRaster into a tidy data frame, tagged
  # with Model, Region and Period. A NULL region (e.g. a diff with no
  # overlapping variables for that region) is skipped rather than erroring.
  imap(all_stacks, function(regions, period_label) {
    imap(regions, function(r, region_label) {
      if (is.null(r)) return(NULL)
      df <- as.data.frame(r, xy = FALSE, na.rm = TRUE)
      if (nrow(df) == 0) return(NULL)
      df %>% mutate(Model = model, Region = str_to_title(region_label), Period = period_label)
    }) %>% compact() %>% bind_rows()
  }) %>% compact() %>% bind_rows()
})

all_df <- bind_rows(compact(model_dfs)) %>%
  left_join(model_meta, by = "Model")

# ---- 5. Summary tables: mean +/- sd, by Driver (storyline) and by Model ----

summarise_group <- function(df, group_col) {
  df %>%
    group_by(.data[[group_col]]) %>%
    summarise(across(where(is.numeric),
                     list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE)),
                     .names = "{.col}_{.fn}"),
              .groups = "drop") %>%
    pivot_longer(cols = -all_of(group_col), names_to = c("variable", ".value"), names_sep = "_") %>%
    mutate(mean_sd = paste0(round(mean, 2), " +/- ", round(sd, 2))) %>%
    select(all_of(group_col), variable, mean_sd) %>%
    pivot_wider(names_from = all_of(group_col), values_from = mean_sd)
}

region_scopes <- list(All = NULL, Peninsula = "Peninsula", Continent = "Continent")

for (period_label in unique(all_df$Period)) {
  period_df <- all_df %>% filter(Period == period_label)
  
  for (scope_name in names(region_scopes)) {
    scope_region <- region_scopes[[scope_name]]
    scope_df <- if (is.null(scope_region)) period_df else period_df %>% filter(Region == scope_region)
    if (nrow(scope_df) == 0) next
    
    # Storyline (by-driver) table: CESM2/MPI_ESM1 only -- ERA5 is excluded
    # since it's an evaluation run, not a storyline.
    storyline_df <- scope_df %>% filter(Driver %in% storyline_drivers)
    if (nrow(storyline_df) > 0) {
      by_driver <- summarise_group(storyline_df, "Driver")
      write.csv(by_driver, file.path(outpath, sprintf("%s_%s_by_driver.csv", period_label, scope_name)), row.names = FALSE)
    }
    
    # Per-model table: every model, ERA5 included.
    by_model <- summarise_group(scope_df, "Model")
    write.csv(by_model, file.path(outpath, sprintf("%s_%s_by_model.csv", period_label, scope_name)), row.names = FALSE)
    
    message("Wrote ", period_label, " / ", scope_name, " summary tables")
  }
}

message("\nDone. Summary tables written to: ", outpath)

