# ============================================================
# AntAirICE vs HCLIM / RACMO / MetUM temperature evaluation
# Combined Storyline 1 / Storyline 2 / ERA5 Re-analysis figure
# ============================================================
#
# ASSUMPTIONS (check / adjust before running):
#  1. Model rasters all follow:
#       Data/Environmental_predictors/PolarRes26/Regridded/<MODEL>_<DRIVING>/comparison/
#         Mean_<Season>_Temperature_HISTORICAL_2003_2014_ICEFREE.tif
#     e.g. HCLIM_ERA5, RACMO_ERA5, MetUM_ERA5, HCLIM_MPI_ESM1, etc.
#  2. MetUM only exists for the ERA5-driven run (per your layout),
#     so it only appears in row 3. Rows 1-2 (Storyline 1/2) have
#     just HCLIM + RACMO, and those two panels stretch to fill the
#     same overall width as row 3's three panels (they end up a
#     little taller as a result, which is expected/fine).
#  3. No AntAirICE winter observations exist yet, so the winter
#     figure is a placeholder: same row/column template, each panel
#     shows a light-grey "Data pending" tile instead of real data.
#     Swap in `Mean_Winter_Temp_ICEFREE.tif` and re-run once that
#     data exists - the real scatter-density code path is already
#     wired up, it's just skipped for "Winter" below.
#  4. Figure sizing/text sizes are tuned to stay legible when the
#     PNG is pasted into an A4 page (see ggsave() calls at the end).
#  5. The fill (grid-cell count) legend is a single continuous scale
#     computed ONCE from the combined Annual + Summer data, and that
#     same scale (axis limits, bin edges, and colour limits) is reused
#     for every season's figure - so the legend/colours mean exactly
#     the same thing whether you're looking at the Annual, Summer, or
#     Winter (placeholder) PNG.
#
# ============================================================

library(terra)
library(ggplot2)
library(dplyr)
library(purrr)
library(patchwork)
library(here)
library(viridis)
library(Metrics)

# ---------------------------------------------------------------
# OUTPATH
# ---------------------------------------------------------------
outpath <- here("Plots/Evaluation_AntAirICE")

# ---------------------------------------------------------------
# HELPERS (unchanged logic from original script)
# ---------------------------------------------------------------
align_to_model <- function(ref, mod) {
  if (!compareGeom(ref, mod, stopOnError = FALSE)) {
    ref <- resample(ref, mod, method = "bilinear")
  }
  ref
}

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

# ---------------------------------------------------------------
# STYLE (approximating the reference figure: clean sans-serif,
# bold left-aligned titles, muted gridlines/axis text)
# ---------------------------------------------------------------
# Row-title colour lives here so it's defined once and reused by both
# the real panels (build_row_plot) and the placeholder panels
# (build_placeholder_row) - a "darkish grey", not near-black.
row_title_colour <- "grey35"

theme_storyline <- function(base_size = 11, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold", size = 10.5, colour = "grey20"),
      plot.title         = element_text(face = "bold", size = 12.5, hjust = 0,
                                        colour = "grey15", margin = margin(b = 6)),
      axis.title         = element_text(size = 9.5, colour = "grey30"),
      axis.text          = element_text(size = 8.5, colour = "grey45"),
      axis.line          = element_line(colour = "grey60", linewidth = 0.35),
      axis.ticks         = element_line(colour = "grey60", linewidth = 0.35),
      panel.border      = element_rect(colour = "grey60", fill = NA, linewidth = 0.35),
      panel.spacing      = unit(0.45, "cm"),
      legend.title       = element_text(size = 9.5, colour = "grey20"),
      legend.text        = element_text(size = 8.5, colour = "grey40")
    )
}

# ---------------------------------------------------------------
# MODEL CONFIGURATION
#   row_group : which storyline row the panel belongs to
#   col_label : which model column (sub-heading) the panel belongs to
#   folder    : subfolder name under .../Regridded/
# ---------------------------------------------------------------
row_levels <- c("Storyline 1", "Storyline 2", "ERA5 Re-analysis")
# Column order per row comes from the order models appear in
# `model_config` below (HCLIM, RACMO, [MetUM]) - no separate
# level list needed now that rows only build the panels they have.

model_config <- tribble(
  ~model,   ~driving,     ~row_group,
  "HCLIM",  "MPI_ESM1",   "Storyline 1",
  "RACMO",  "MPI_ESM1",   "Storyline 1",
  "HCLIM",  "CESM2",      "Storyline 2",
  "RACMO",  "CESM2",      "Storyline 2",
  "HCLIM",  "ERA5",       "ERA5 Re-analysis",
  "RACMO",  "ERA5",       "ERA5 Re-analysis",
  "MetUM",  "ERA5",       "ERA5 Re-analysis"
) %>%
  mutate(
    col_label = model,
    folder    = paste0(model, "_", driving)
  )

# ---------------------------------------------------------------
# SEASON DEFINITIONS
# ---------------------------------------------------------------
seasons <- c("Annual", "Summer", "Winter")

ref_files <- c(
  Annual = "Mean_Annual_Temp_ICEFREE.tif",
  Summer = "Mean_Summer_Temp_ICEFREE.tif"
  # Winter: not yet available - add e.g. "Mean_Winter_Temp_ICEFREE.tif"
  # here and remove "Winter" from `placeholder_seasons` below once it exists.
)

# ---------------------------------------------------------------
# BUILD PAIRED DATA FOR ONE SEASON
# ---------------------------------------------------------------
build_season_pairs <- function(season) {
  ref <- rast(here("Data/Environmental_predictors", ref_files[[season]]))
  
  pmap_dfr(model_config, function(model, driving, row_group, col_label, folder) {
    mod_path <- here(
      "Data/Environmental_predictors/PolarRes26/Regridded",
      folder, "comparison",
      paste0("Mean_", season, "_Temperature_HISTORICAL_2003_2014_ICEFREE.tif")
    )
    mod <- rast(mod_path)
    df <- extract_pairs(ref, mod, label = paste(model, driving, sep = "_"))
    df$row_group <- row_group
    df$col_label <- col_label
    df
  })
}

# ---------------------------------------------------------------
# SHARED FILL SCALE LIMIT (so the single collected legend is
# valid across every panel, in every row, AND across every season's
# saved figure - see the global computation before the main loop)
# ---------------------------------------------------------------
get_shared_fill_limit <- function(df, breaks) {
  x_bin <- cut(df$x, breaks = breaks, include.lowest = TRUE)
  y_bin <- cut(df$y, breaks = breaks, include.lowest = TRUE)
  grp   <- interaction(df$row_group, df$col_label, drop = TRUE)
  max(table(grp, x_bin, y_bin))
}

# ---------------------------------------------------------------
# BUILD ONE PANEL (a single model vs AntAirICE)
# ---------------------------------------------------------------
build_panel_plot <- function(df_panel, col_label, ax_lim, ax_breaks,
                             common_breaks, fill_limit) {
  df_panel$col_label <- col_label
  
  ggplot(df_panel, aes(x = x, y = y)) +
    stat_bin2d(breaks = list(x = common_breaks, y = common_breaks),
               aes(fill = after_stat(count))) +
    geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.4) +
    scale_fill_viridis_c(
      option   = "plasma",
      name     = "Number of\ngrid cells",
      limits   = c(1, fill_limit),
      oob      = scales::squish,
      na.value = NA
    ) +
    scale_x_continuous(limits = ax_lim, breaks = ax_breaks) +
    scale_y_continuous(limits = ax_lim, breaks = ax_breaks) +
    coord_fixed() +
    facet_wrap(~ col_label, nrow = 1) +
    labs(
      x = "Temperature (AntAir ICE) [\u00B0C]",
      y = "Temperature (Model) [\u00B0C]"
    ) +
    theme_storyline()
}

# ---------------------------------------------------------------
# BUILD ONE ROW (Storyline 1 / Storyline 2 / ERA5 Re-analysis)
# Each row only contains the panels that actually exist for it
# (2 for the storylines, 3 for ERA5). Because every row is wrapped
# to the SAME overall figure width, a 2-panel row automatically
# stretches those panels to fill the full width - no empty "MetUM"
# slot needed. Panels stay square (coord_fixed), so 2-panel rows
# end up a little taller than the 3-panel row, which is fine.
# ---------------------------------------------------------------
build_row_plot <- function(season_df, row_title, ax_lim, ax_breaks,
                           common_breaks, fill_limit) {
  row_config <- filter(model_config, row_group == row_title)
  
  panels <- lapply(row_config$col_label, function(cl) {
    df_panel <- filter(season_df, row_group == row_title, col_label == cl)
    build_panel_plot(df_panel, cl, ax_lim, ax_breaks, common_breaks, fill_limit)
  })
  
  wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title = row_title,
      theme = theme(
        plot.title = element_text(face = "bold", size = 12.5, hjust = 0,
                                  colour = row_title_colour, margin = margin(b = 6),
                                  family = "Helvetica")
      )
    )
}

# ---------------------------------------------------------------
# BUILD FULL SEASON FIGURE (3 rows, shared legend, no overall title)
# ax_lim / ax_breaks / common_breaks / fill_limit are now passed in
# from the GLOBAL computation below, rather than recomputed per
# season, so the continuous fill legend is identical across the
# Annual, Summer, and Winter figures.
# ---------------------------------------------------------------
make_season_plot <- function(season_df, ax_lim, ax_breaks, common_breaks, fill_limit) {
  row_plots <- lapply(row_levels, function(rg) {
    build_row_plot(season_df, rg, ax_lim, ax_breaks, common_breaks, fill_limit)
  })
  
  wrap_plots(row_plots, ncol = 1) +
    plot_layout(guides = "collect") &
    theme(
      legend.position   = "right",
      legend.key.height = unit(1.6, "cm"),
      legend.key.width  = unit(0.35, "cm")
    )
}

# ---------------------------------------------------------------
# PLACEHOLDER FIGURE (used while a season has no observed data,
# e.g. winter AntAirICE not yet available). Keeps the same row
# titles / model sub-headings as the real figures, and the same
# "stretch to fill the row" behaviour for the 2-panel rows, so it
# drops in seamlessly once real data exists - just remove this
# season from `placeholder_seasons` below.
# ---------------------------------------------------------------
placeholder_seasons <- c("Winter")

build_placeholder_panel <- function(col_label) {
  df <- data.frame(col_label = col_label, x = 0.5, y = 0.5)
  
  ggplot(df, aes(x = x, y = y)) +
    geom_text(label = "Data pending", size = 3.6, fontface = "italic",
              colour = "grey55", family = "Helvetica") +
    facet_wrap(~ col_label, nrow = 1) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    coord_fixed() +
    theme_storyline() +
    theme(
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.title       = element_blank(),
      axis.line        = element_blank(),
      panel.background = element_rect(fill = "grey94", colour = NA)
    )
}

build_placeholder_row <- function(row_title) {
  row_config <- filter(model_config, row_group == row_title)
  panels <- lapply(row_config$col_label, build_placeholder_panel)
  
  wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title = row_title,
      theme = theme(
        plot.title = element_text(face = "bold", size = 12.5, hjust = 0,
                                  colour = row_title_colour, margin = margin(b = 6),
                                  family = "Helvetica")
      )
    )
}

make_placeholder_plot <- function() {
  row_plots <- lapply(row_levels, build_placeholder_row)
  wrap_plots(row_plots, ncol = 1)
}

# ---------------------------------------------------------------
# GLOBAL SCALE (computed once, from the combined Annual + Summer
# data, and reused for every season's figure below). This is what
# makes the "Number of grid cells" legend continuous and identical
# across the Annual / Summer / Winter(placeholder) PNGs, instead of
# each figure silently rescaling to its own max count.
#
# The legend is additionally hard-capped at 400: if the raw shared
# max is higher, fill_limit is clamped to 400 and any bin with more
# than 400 grid cells is squished (scales::squish, set in
# build_panel_plot) into the top colour rather than being dropped
# as NA or stretching the scale to a rarely-hit outlier bin.
# ---------------------------------------------------------------
real_seasons <- setdiff(seasons, placeholder_seasons)

season_pairs_cache <- setNames(
  lapply(real_seasons, build_season_pairs),
  real_seasons
)

all_pairs <- bind_rows(season_pairs_cache)

raw_lim   <- range(c(all_pairs$x, all_pairs$y), na.rm = TRUE)
ax_lim    <- c(floor(raw_lim[1] / 5) * 5, ceiling(raw_lim[2] / 5) * 5)
ax_breaks <- seq(ax_lim[1], ax_lim[2], by = 10)

# identical bin edges everywhere -> counts are directly comparable
common_breaks <- seq(ax_lim[1], ax_lim[2], length.out = 151)  # 150 bins

fill_limit <- min(get_shared_fill_limit(all_pairs, common_breaks), 400)

# ---------------------------------------------------------------
# RUN FOR EACH SEASON AND SAVE
# Sizing is tuned for pasting/printing at A4 (210 x 297 mm /
# 8.27 x 11.69 in): a bit narrower and shorter than the page so
# margins remain, with a high dpi so text stays crisp when scaled.
# ---------------------------------------------------------------
a4_width  <- 7.8   # inches
a4_height <- 10.5  # inches
a4_dpi    <- 320

for (season in seasons) {
  message("Processing ", season, " ...")
  
  if (season %in% placeholder_seasons) {
    season_plot <- make_placeholder_plot()
  } else {
    season_pairs <- season_pairs_cache[[season]]
    season_plot  <- make_season_plot(season_pairs, ax_lim, ax_breaks,
                                     common_breaks, fill_limit)
  }
  
  ggsave(
    file.path(outpath, paste0("AntAirICE_vs_Models_", tolower(season), ".png")),
    season_plot,
    width  = a4_width,
    height = a4_height,
    dpi    = a4_dpi,
    bg     = "white"
  )
  message(season, " plot saved.")
}



############################################################################

# ============================================================
# PLOT 2: Bias raster maps  (2 columns x 2 rows)
# ============================================================



# ============================================================
# AntAirICE vs HCLIM / RACMO / MetUM - MEAN BIAS MAPS
# Same Storyline 1 / Storyline 2 / ERA5-Reanalysis row layout as
# AntAirICE_vs_Models_plots.R, showing spatial bias maps
# (Model - AntAirICE) instead of scatter-density panels, styled
# after the CHELSA-vs-HCLIM bias-map script (coastline underlay,
# 10 km aggregation, zero-centred diverging colour scale).
# ============================================================
#
# ASSUMPTIONS (check / adjust before running):
#  1. Model rasters follow the same convention as the scatter-density
#     script:
#       Data/Environmental_predictors/PolarRes26/Regridded/<MODEL>_<DRIVING>/comparison/
#         Mean_<Season>_Temperature_HISTORICAL_2003_2014_ICEFREE.tif
#  2. MetUM only exists for the ERA5-driven run, so it only appears
#     in row 3; rows 1-2 (2 panels) stretch to fill the same overall
#     width as row 3 (3 panels) - exactly as in the scatter-density
#     script, and for the same reason (no empty "MetUM" slot needed).
#  3. No winter AntAirICE observations exist yet, so winter renders
#     as a placeholder ("Data pending" tiles) - swap in the real file
#     in `ref_files` and remove "Winter" from `placeholder_seasons`
#     once it's available.
#  4. Coastline shapefile path follows the CHELSA script:
#       Data/add_coastline_medium_res_polygon_v7_10.shp
#  5. Aggregation to 10 km (fact = 10) assumes native resolution is
#     roughly 1 km, as in the CHELSA comparison script - adjust
#     `agg_fact` below if your model grids are a different resolution.
#
# ============================================================

library(terra)
library(ggplot2)
library(dplyr)
library(purrr)
library(patchwork)
library(here)
library(sf)
library(scales)

# ---------------------------------------------------------------
# OUTPATH
# ---------------------------------------------------------------
outpath <- here("Plots/Evaluation_AntAirICE")
if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

# ---------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------
align_to_model <- function(ref, mod) {
  if (!compareGeom(ref, mod, stopOnError = FALSE)) {
    ref <- resample(ref, mod, method = "bilinear")
  }
  ref
}

rast_to_df <- function(r, val_name = "value") {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- val_name
  df
}

agg_fact <- 10  # aggregation factor - coarser cells "stand out" more visually

# ---------------------------------------------------------------
# DIVERGING COLOUR SCALE - zero always lines up with white, the
# range doesn't have to be symmetric. Same construction as the
# CHELSA vs HCLIM script.
# ---------------------------------------------------------------
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

make_cbar <- function() {
  guide_colorbar(
    barheight       = unit(3.2, "cm"),
    barwidth        = unit(0.4, "cm"),
    ticks.colour    = "black",
    frame.colour    = "black",
    frame.linewidth = 0.4
  )
}

# ---------------------------------------------------------------
# STYLE - no panel border/outline anywhere ("no outline around the
# plots"); bold left-aligned titles to match the scatter-density
# script's look.
# ---------------------------------------------------------------
theme_bias_map <- function(base_size = 11, base_family = "Helvetica") {
  theme_void(base_size = base_size) +
    theme(
      text         = element_text(family = base_family),
      plot.title   = element_text(face = "bold", size = 10.5, hjust = 0,
                                  colour = "grey20", margin = margin(b = 4)),
      legend.title = element_text(size = 9.5, colour = "grey20"),
      legend.text  = element_text(size = 8.5, colour = "grey40"),
      plot.margin  = margin(4, 4, 4, 4)
    )
}

row_title_theme <- theme(
  plot.title = element_text(face = "bold", size = 12.5, hjust = 0,
                            colour = "grey15", margin = margin(b = 6),
                            family = "Helvetica")
)

# ---------------------------------------------------------------
# COASTLINE - lighter grey than the CHELSA script's grey85, no
# outline of its own so it sits quietly behind the bias tiles.
# ---------------------------------------------------------------
coast <- st_read(
  here("Data/add_coastline_medium_res_polygon_v7_10.shp"),
  quiet = TRUE
)
coast_fill <- "grey95"

# ---------------------------------------------------------------
# MODEL CONFIGURATION (identical to the scatter-density script)
# ---------------------------------------------------------------
row_levels <- c("Storyline 1", "Storyline 2", "ERA5 Reanalysis")

model_config <- tribble(
  ~model,   ~driving,     ~row_group,
  "HCLIM",  "MPI_ESM1",   "Storyline 1",
  "RACMO",  "MPI_ESM1",   "Storyline 1",
  "HCLIM",  "CESM2",      "Storyline 2",
  "RACMO",  "CESM2",      "Storyline 2",
  "HCLIM",  "ERA5",       "ERA5 Reanalysis",
  "RACMO",  "ERA5",       "ERA5 Reanalysis",
  "MetUM",  "ERA5",       "ERA5 Reanalysis"
) %>%
  mutate(
    col_label = model,
    folder    = paste0(model, "_", driving)
  )

# ---------------------------------------------------------------
# SEASON DEFINITIONS
# ---------------------------------------------------------------
seasons <- c("Annual", "Summer", "Winter")
placeholder_seasons <- c("Winter")

ref_files <- c(
  Annual = "Mean_Annual_Temp_ICEFREE.tif",
  Summer = "Mean_Summer_Temp_ICEFREE.tif"
  # Winter: not yet available - add here + remove from
  # placeholder_seasons above once it exists.
)

# ---------------------------------------------------------------
# BUILD BIAS RASTERS FOR ONE SEASON (Model - AntAirICE, aggregated
# to 10 km)
# ---------------------------------------------------------------
build_season_bias <- function(season) {
  ref <- rast(here("Data/Environmental_predictors", ref_files[[season]]))
  
  pmap(model_config, function(model, driving, row_group, col_label, folder) {
    mod_path <- here(
      "Data/Environmental_predictors/PolarRes26/Regridded",
      folder, "comparison",
      paste0("Mean_", season, "_Temperature_HISTORICAL_2003_2014_ICEFREE.tif")
    )
    mod       <- rast(mod_path)
    diff      <- mod - align_to_model(ref, mod)
    diff_10km <- aggregate(diff, fact = agg_fact, fun = "mean", na.rm = TRUE)
    
    list(row_group = row_group, col_label = col_label, r = diff_10km)
  })
}

# ---------------------------------------------------------------
# SHARED DIVERGING SCALE FOR ONE SEASON - every panel in that
# season's figure uses this single scale, so the one collected
# legend is valid everywhere.
# ---------------------------------------------------------------
get_season_scale <- function(bias_list) {
  all_vals <- unlist(lapply(bias_list, function(b) values(b$r, na.rm = TRUE)))
  diff_min <- floor(min(all_vals, na.rm = TRUE))
  diff_max <- ceiling(max(all_vals, na.rm = TRUE))
  list(
    min    = diff_min,
    max    = diff_max,
    breaks = pretty(c(diff_min, diff_max), n = 6),
    values = make_scale_values(diff_min, diff_max)
  )
}

# ---------------------------------------------------------------
# BUILD ONE BIAS PANEL (a single model's bias map vs AntAirICE)
# ---------------------------------------------------------------
build_bias_panel <- function(r, col_label, scale_info) {
  df <- rast_to_df(r, "bias")
  
  ggplot() +
    geom_sf(data = coast, fill = coast_fill, colour = "grey35", linewidth = 0.15) +
    geom_tile(data = df, aes(x = x, y = y, fill = bias)) +
    scale_fill_gradientn(
      colours = diverging_ramp,
      values  = scale_info$values,
      limits  = c(scale_info$min, scale_info$max),
      breaks  = scale_info$breaks,
      name    = "Bias (\u00B0C)",
      oob     = squish,
      guide   = make_cbar()
    ) +
    coord_sf(expand = FALSE) +
    labs(title = col_label) +
    theme_bias_map()
}

# ---------------------------------------------------------------
# BUILD ONE ROW (mirrors build_row_plot() in the scatter-density
# script: only the panels that exist for this row are built, and
# because every row is stacked into the same overall figure width,
# a 2-panel row automatically stretches to fill it - no empty
# "MetUM" slot needed for rows 1-2).
# ---------------------------------------------------------------
build_bias_row <- function(row_title, bias_list, scale_info) {
  row_items <- Filter(function(b) b$row_group == row_title, bias_list)
  
  panels <- lapply(row_items, function(b) {
    build_bias_panel(b$r, b$col_label, scale_info)
  })
  
  wrap_plots(panels, nrow = 1) +
    plot_annotation(title = row_title, theme = row_title_theme)
}

# ---------------------------------------------------------------
# BUILD FULL SEASON BIAS FIGURE (3 rows, one shared legend, no
# overall title)
# ---------------------------------------------------------------
make_season_bias_plot <- function(season) {
  bias_list  <- build_season_bias(season)
  scale_info <- get_season_scale(bias_list)
  
  row_plots <- lapply(row_levels, function(rg) {
    build_bias_row(rg, bias_list, scale_info)
  })
  
  wrap_plots(row_plots, ncol = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
}

# ---------------------------------------------------------------
# PLACEHOLDER FIGURE (winter - no observations yet)
# ---------------------------------------------------------------
build_placeholder_panel <- function(col_label) {
  df <- data.frame(x = 0.5, y = 0.5, label = "Data pending")
  ggplot(df, aes(x = x, y = y, label = label)) +
    geom_text(size = 3.6, fontface = "italic", colour = "grey55",
              family = "Helvetica") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    coord_fixed() +
    labs(title = col_label) +
    theme_bias_map() +
    theme(panel.background = element_rect(fill = "grey94", colour = NA))
}

build_placeholder_row <- function(row_title) {
  row_config <- filter(model_config, row_group == row_title)
  panels <- lapply(row_config$col_label, build_placeholder_panel)
  
  wrap_plots(panels, nrow = 1) +
    plot_annotation(title = row_title, theme = row_title_theme)
}

make_placeholder_plot <- function() {
  row_plots <- lapply(row_levels, build_placeholder_row)
  wrap_plots(row_plots, ncol = 1)
}

# ---------------------------------------------------------------
# RUN FOR EACH SEASON AND SAVE (A4-friendly sizing, matching the
# scatter-density script)
# ---------------------------------------------------------------
a4_width  <- 7.8
a4_height <- 10.5
a4_dpi    <- 320

for (season in seasons) {
  message("Processing bias maps: ", season, " ...")
  
  season_plot <- if (season %in% placeholder_seasons) {
    make_placeholder_plot()
  } else {
    make_season_bias_plot(season)
  }
  
  ggsave(
    file.path(outpath, paste0("AntAirICE_vs_Models_bias_", tolower(season), ".png")),
    season_plot,
    width  = a4_width,
    height = a4_height,
    dpi    = a4_dpi,
    bg     = "white"
  )
  message(season, " bias plot saved.")
}


##############################################################################
################ PART 3 - Summary Statistics ################################
##############################################################################

# ============================================================
# MODEL PERFORMANCE METRICS
# Mean Bias, RMSE, MAE and R2
# ============================================================

# ---------------------------------------------------------------
# CALCULATE METRICS FOR ONE MODEL / ONE SEASON
# ---------------------------------------------------------------
calculate_metrics <- function(ref, mod) {
  
  # Align AntAirICE to the model grid, exactly as in the
  # scatter-density plots
  ref <- align_to_model(ref, mod)
  
  # Extract paired values
  vals <- data.frame(
    actual    = as.vector(values(ref)),
    predicted = as.vector(values(mod))
  )
  
  # Keep only grid cells where BOTH are non-NA
  vals <- vals[complete.cases(vals), ]
  
  # Mean bias: model - AntAirICE
  mean_bias <- mean(vals$predicted - vals$actual)
  
  # RMSE
  rmse <- Metrics::rmse(
    actual    = vals$actual,
    predicted = vals$predicted
  )
  
  # Mean absolute error
  mae <- Metrics::mae(
    actual    = vals$actual,
    predicted = vals$predicted
  )
  
  # R2 from a linear model
  lm_fit <- lm(predicted ~ actual, data = vals)
  r2 <- summary(lm_fit)$r.squared
  
  data.frame(
    n_grid_cells = nrow(vals),
    mean_bias    = mean_bias,
    RMSE         = rmse,
    MAE          = mae,
    R2           = r2
  )
}


# ---------------------------------------------------------------
# CALCULATE AND SAVE METRICS
# ---------------------------------------------------------------

# Only seasons with an AntAirICE reference raster
metric_seasons <- names(ref_files)

for (season in metric_seasons) {
  
  message("Calculating metrics for ", season, " ...")
  
  # AntAirICE reference
  ref <- rast(
    here(
      "Data/Environmental_predictors",
      ref_files[[season]]
    )
  )
  
  # Calculate metrics for every model configuration
  metrics_df <- pmap_dfr(
    model_config,
    function(model, driving, row_group, col_label, folder) {
      
      mod_path <- here(
        "Data/Environmental_predictors/PolarRes26/Regridded",
        folder,
        "comparison",
        paste0(
          "Mean_",
          season,
          "_Temperature_HISTORICAL_2003_2014_ICEFREE.tif"
        )
      )
      
      mod <- rast(mod_path)
      
      metrics <- calculate_metrics(ref, mod)
      
      metrics$model     <- model
      metrics$driving   <- driving
      metrics$row_group <- row_group
      metrics$season    <- season
      
      metrics
    }
  )
  
  # Put identifying columns first
  metrics_df <- metrics_df %>%
    select(
      season,
      row_group,
      model,
      driving,
      n_grid_cells,
      mean_bias,
      RMSE,
      MAE,
      R2
    )
  
  # Save one CSV per season
  write.csv(
    metrics_df,
    file.path(
      outpath,
      paste0(
        "AntAirICE_model_metrics_",
        tolower(season),
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  message(season, " metrics saved.")
}


##############################################################################
################ PART 3B - MEAN BIAS PLOT ################################
##############################################################################

# ===============================================================
# BIAS HEATMAP
#
# Rows:
#   Storyline 1
#   Storyline 2
#   ERA5
#
# Column groups:
#   HCLIM   ALL DJF JJA
#   RACMO   ALL DJF JJA
#   MetUM   ALL DJF JJA
# ===============================================================


# ---------------------------------------------------------------
# LOAD METRICS OUTPUT
# ---------------------------------------------------------------

metric_files <- c(
  Annual = file.path(outpath, "AntAirICE_model_metrics_annual.csv"),
  Summer = file.path(outpath, "AntAirICE_model_metrics_summer.csv"),
  Winter = file.path(outpath, "AntAirICE_model_metrics_winter.csv")
)

# Load files that exist
bias_df <- bind_rows(
  lapply(names(metric_files), function(season) {
    
    if (file.exists(metric_files[[season]])) {
      df <- read.csv(metric_files[[season]])
      df$season <- season
      df
      
    } else {
      data.frame()
    }
  })
)


# ---------------------------------------------------------------
# ADD EMPTY WINTER DATA IF WINTER FILE DOES NOT EXIST YET
# ---------------------------------------------------------------

if (!"Winter" %in% bias_df$season) {
  
  winter_empty <- tribble(
    ~model,   ~driving,
    "HCLIM",  "MPI_ESM1",
    "RACMO",  "MPI_ESM1",
    "HCLIM",  "CESM2",
    "RACMO",  "CESM2",
    "HCLIM",  "ERA5",
    "RACMO",  "ERA5",
    "MetUM",  "ERA5"
  )
  
  winter_empty$season <- "Winter"
  winter_empty$mean_bias <- NA_real_
  
  bias_df <- bind_rows(
    bias_df,
    winter_empty
  )
}


# ---------------------------------------------------------------
# DEFINE MODEL / STORYLINE STRUCTURE
# ---------------------------------------------------------------

model_structure <- tribble(
  ~row_group,          ~model,  ~driving,    ~model_group,
  "Storyline 1",       "HCLIM", "MPI_ESM1",  "HCLIM",
  "Storyline 1",       "RACMO", "MPI_ESM1",  "RACMO",
  "Storyline 2",       "HCLIM", "CESM2",     "HCLIM",
  "Storyline 2",       "RACMO", "CESM2",     "RACMO",
  "ERA5",              "HCLIM", "ERA5",      "HCLIM",
  "ERA5",              "RACMO", "ERA5",      "RACMO",
  "ERA5",              "MetUM", "ERA5",      "MetUM"
)


# ---------------------------------------------------------------
# ADD MODEL / STORYLINE INFORMATION
# ---------------------------------------------------------------

bias_df <- bias_df %>%
  mutate(
    model_key = paste(model, driving, sep = "_")
  ) %>%
  left_join(
    model_structure,
    by = c("model", "driving")
  ) %>%
  mutate(
    season_label = recode(
      season,
      "Annual" = "ALL",
      "Summer" = "DJF",
      "Winter" = "JJA"
    )
  )


# ---------------------------------------------------------------
# CREATE COMPLETE 3 x 3 x 3 GRID
#
# This ensures:
#   - all three seasons are present
#   - MetUM cells are empty for Storyline 1 and Storyline 2
#   - JJA cells are empty until winter data exist
# ---------------------------------------------------------------

plot_grid <- expand_grid(
  row_group = c(
    "Storyline 1",
    "Storyline 2",
    "ERA5"
  ),
  model_group = c(
    "HCLIM",
    "RACMO",
    "MetUM"
  ),
  season_label = c(
    "ALL",
    "DJF",
    "JJA"
  )
) %>%
  left_join(
    model_structure,
    by = c("row_group", "model_group")
  ) %>%
  left_join(
    bias_df %>%
      select(
        model,
        driving,
        season_label,
        mean_bias
      ),
    by = c(
      "model",
      "driving",
      "season_label"
    )
  )


# ---------------------------------------------------------------
# SET ORDER
# ---------------------------------------------------------------

plot_grid$row_group <- factor(
  plot_grid$row_group,
  levels = rev(c(
    "Storyline 1",
    "Storyline 2",
    "ERA5"
  ))
)

plot_grid$model_group <- factor(
  plot_grid$model_group,
  levels = c(
    "HCLIM",
    "RACMO",
    "MetUM"
  )
)

plot_grid$season_label <- factor(
  plot_grid$season_label,
  levels = c(
    "ALL",
    "DJF",
    "JJA"
  )
)


# ---------------------------------------------------------------
# SYMMETRIC BIAS SCALE
# ---------------------------------------------------------------

max_bias <- max(
  abs(plot_grid$mean_bias),
  na.rm = TRUE
)

# Round limit up to nearest 0.5
bias_limit <- ceiling(max_bias * 2) / 2


# ---------------------------------------------------------------
# BUILD HEATMAP
# ---------------------------------------------------------------

bias_plot <- ggplot(
  plot_grid,
  aes(
    x = season_label,
    y = row_group
  )
) +
  
  geom_tile(
    aes(fill = mean_bias),
    colour = "white",
    linewidth = 0.8,
    na.rm = FALSE
  ) +
  
  geom_text(
    aes(
      label = ifelse(
        is.na(mean_bias),
        "",
        sprintf("%.2f", mean_bias)
      ),
      colour = ifelse(
        is.na(mean_bias),
        "black",
        ifelse(
          abs(mean_bias) > bias_limit * 0.45,
          "white",
          "black"
        )
      )
    ),
    size = 4
  ) +
  
  scale_colour_identity() +
  
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = c(
      -bias_limit,
      bias_limit
    ),
    name = "Mean bias (\u00B0C)",
    na.value = "grey95"
  ) +
  
  facet_grid(
    . ~ model_group
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  coord_fixed() +
  
  theme_classic() +
  
  theme(
    # Row labels
    axis.text.y = element_text(
      size = 11.5,
      face = "bold",
      colour = "black"
    ),
    
    # ALL / DJF / JJA
    axis.text.x = element_text(
      size = 10.5,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      colour = "grey60",
      linewidth = 0.35
    ),
    
    axis.line = element_blank(),
    
    # HCLIM / RACMO / MetUM headings
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 10.5,
      colour = "black"
    ),
    
    # Remove gaps between the three column groups
    panel.spacing.x = unit(0.8, "lines"),
    
    legend.title = element_text(
      size = 10.5
    ),
    
    legend.text = element_text(
      size = 9.5
    ),
    
    legend.key.height = unit(1, "cm"),
    
    legend.key.width = unit(0.35, "cm"),
    
    plot.margin = margin(
      5, 5, 5, 5
    )
  )


# ---------------------------------------------------------------
# SAVE
# ---------------------------------------------------------------

ggsave(
  file.path(
    outpath,
    "AntAirICE_model_bias_heatmap_grouped.png"
  ),
  bias_plot,
  width = 7.5,
  height = 3.8,
  dpi = 300
)





