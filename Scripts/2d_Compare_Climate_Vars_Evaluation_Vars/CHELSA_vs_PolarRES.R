# ============================================================
# CHELSA vs HCLIM / RACMO — FUTURE ONLY spatial bias comparison
# PENINSULA-ONLY VERSION
#
# Builds on the per-panel-legend version, with these additions:
#   1. Loads Data/Peninsula_Continent_Boundary.shp and uses it to
#      crop + mask every raster (CHELSA ensemble mean, each model's
#      temperature, each model's bias) down to just the Peninsula,
#      rather than plotting the full continent and zooming.
#      -> Cropping BEFORE computing the shared colour scales means
#         the temperature/difference colour ramps are stretched over
#         the Peninsula's actual value range, instead of wasting most
#         of the ramp on continent-wide extremes that don't occur here.
#   2. The full-continent coastline layer is swapped for the Peninsula
#      boundary polygon itself as the basemap outline in every panel.
#   3. `agg_fact` is reduced (10 -> 5 km) since the Peninsula covers a
#      much smaller area than the full continent, so a coarser 10 km
#      aggregation would throw away most of the detail. Adjust back up
#      if the plot renders slowly.
#   4. Saved to a different output filename so it doesn't overwrite
#      the full-continent plot.
#
# ASSUMPTION TO CHECK: the Peninsula boundary shapefile is assumed to
# live at Data/Peninsula_Continent_Boundary.shp (same folder as the
# existing coastline shapefile). Adjust `peninsula_path` below if it's
# stored elsewhere.
# ============================================================

library(terra)
library(ggplot2)
library(dplyr)
library(tibble)
library(purrr)
library(patchwork)
library(here)
library(scales)
library(sf)

# ============================================================
# OUTPUT PATH
# ============================================================

outpath <- here("Plots/Evaluation_CHELSA")
if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

# ============================================================
# SIZING
# ============================================================

PAGE_W <- 8.27
PAGE_H <- 11.69
DPI    <- 320
BASE   <- 8.5
FONT   <- "Helvetica"

LEGBAR_H <- unit(1.9, "cm")
LEGBAR_W <- unit(0.3,  "cm")

agg_fact <- 5  # aggregate to ~5 km — Peninsula is much smaller than full continent

# ============================================================
# HELPERS
# ============================================================

rast_to_df <- function(r, val_name = "value") {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- val_name
  df
}

align_to_model <- function(ref, mod) {
  if (!compareGeom(ref, mod, stopOnError = FALSE)) {
    ref <- resample(ref, mod, method = "bilinear")
  }
  ref
}

# Crop + mask a raster down to the Peninsula boundary polygon
crop_mask_peninsula <- function(r, boundary_vect) {
  r <- crop(r, boundary_vect)
  r <- mask(r, boundary_vect)
  r
}

theme_ant <- function(base_size = BASE, base_family = FONT) {
  theme_void(base_size = base_size, base_family = base_family) +
    theme(
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.background   = element_rect(fill = "white", colour = NA),
      legend.position    = "right",
      legend.box.spacing = unit(1, "pt"),
      legend.spacing     = unit(1, "pt"),
      legend.title       = element_text(size = base_size, angle = 90,
                                        hjust = 0.5, vjust = 0.5),
      legend.text        = element_text(size = base_size - 1, margin = margin(l = 1)),
      legend.key.height  = LEGBAR_H,
      legend.key.width   = LEGBAR_W,
      legend.margin      = margin(0, 0, 0, 0),
      plot.title         = element_text(size = base_size + 1, face = "bold",
                                        hjust = 0.4, vjust = 1, colour = "grey15",
                                        margin = margin(b = 1)),
      plot.margin        = margin(1, 1, 1, 1)
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
    frame.linewidth = 0.35
  )
}

row_title_theme <- theme(
  plot.title = element_text(face = "bold", size = BASE + 2, hjust = 0.4,
                            colour = "grey15", margin = margin(b = 2),
                            family = FONT)
)

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

# ============================================================
# PENINSULA BOUNDARY — used both as the basemap outline and as the
# crop/mask extent for every raster
# ============================================================

peninsula_path     <- here("Data/Peninsula_Continent_Boundary.shp")
peninsula_boundary <- st_read(peninsula_path, quiet = TRUE)
peninsula_vect     <- vect(peninsula_boundary)
peninsula_bbox      <- st_bbox(peninsula_boundary)

coast_fill    <- "grey95"
coast_colour  <- "grey35"

# ============================================================
# PANEL BUILDERS — basemap swapped to peninsula_boundary, view locked
# to its bounding box
# ============================================================

make_temp_panel <- function(r, panel_label, temp_min, temp_max, temp_values) {
  df <- rast_to_df(r, "temp")
  ggplot() +
    geom_sf(data = peninsula_boundary, fill = coast_fill, colour = coast_colour, linewidth = 0.15) +
    geom_tile(data = df, aes(x = x, y = y, fill = temp)) +
    scale_fill_gradientn(
      colours = diverging_ramp,
      values  = temp_values,
      limits  = c(temp_min, temp_max),
      breaks  = pretty(c(temp_min, temp_max), n = 6),
      name    = "Temperature (\u00B0C)",
      oob     = squish,
      guide   = make_cbar()
    ) +
    coord_sf(xlim = c(peninsula_bbox[["xmin"]], peninsula_bbox[["xmax"]]),
             ylim = c(peninsula_bbox[["ymin"]], peninsula_bbox[["ymax"]]),
             expand = FALSE) +
    labs(title = panel_label) +
    theme_ant()
}

make_diff_panel <- function(r, diff_min, diff_max, diff_breaks, diff_values) {
  df <- rast_to_df(r, "diff")
  ggplot() +
    geom_sf(data = peninsula_boundary, fill = coast_fill, colour = coast_colour, linewidth = 0.15) +
    geom_tile(data = df, aes(x = x, y = y, fill = diff)) +
    scale_fill_gradientn(
      colours = diverging_ramp,
      values  = diff_values,
      limits  = c(diff_min, diff_max),
      breaks  = diff_breaks,
      name    = "Difference (\u00B0C)",
      oob     = squish,
      guide   = make_cbar()
    ) +
    coord_sf(xlim = c(peninsula_bbox[["xmin"]], peninsula_bbox[["xmax"]]),
             ylim = c(peninsula_bbox[["ymin"]], peninsula_bbox[["ymax"]]),
             expand = FALSE) +
    labs(title = "Difference") +
    theme_ant()
}

# ============================================================
# MODEL CONFIGURATION
# ============================================================

model_config <- tribble(
  ~model,   ~driving,   ~row_group,
  "HCLIM",  "MPI_ESM1", "Storyline 1",
  "RACMO",  "MPI_ESM1", "Storyline 1",
  "HCLIM",  "CESM2",    "Storyline 2",
  "RACMO",  "CESM2",    "Storyline 2"
) %>%
  mutate(
    col_label = model,
    folder    = paste0(model, "_", driving)
  )

storyline_levels <- c("Storyline 1", "Storyline 2")

# ============================================================
# 1. CHELSA FUTURE ENSEMBLE MEAN (2071-2100, 5 GCMs) — cropped/masked
#    to the Peninsula immediately after building the ensemble mean
# ============================================================

CHELSA_future1 <- rast(here("Data/CHELSA/comparison/Mean_Annual_Temperature_FUTURE_gfdl-esm4_2071_2100_ICEFREE.tif"))
CHELSA_future2 <- rast(here("Data/CHELSA/comparison/Mean_Annual_Temperature_FUTURE_ipsl-cm6a-lr_2071_2100_ICEFREE.tif"))
CHELSA_future3 <- rast(here("Data/CHELSA/comparison/Mean_Annual_Temperature_FUTURE_mpi-esm1-2-hr_2071_2100_ICEFREE.tif"))
CHELSA_future4 <- rast(here("Data/CHELSA/comparison/Mean_Annual_Temperature_FUTURE_mri-esm2-0_2071_2100_ICEFREE.tif"))
CHELSA_future5 <- rast(here("Data/CHELSA/comparison/Mean_Annual_Temperature_FUTURE_ukesm1-0-ll_2071_2100_ICEFREE.tif"))

CHELSA_future <- c(CHELSA_future1, CHELSA_future2, CHELSA_future3,
                   CHELSA_future4, CHELSA_future5)
CHELSA_future <- app(CHELSA_future, mean, na.rm = TRUE)
CHELSA_future <- crop_mask_peninsula(CHELSA_future, peninsula_vect)

# ============================================================
# 2. LOAD EACH FUTURE MODEL, CROP/MASK TO PENINSULA, COMPUTE BIAS
# ============================================================

load_future_model <- function(model, driving, row_group, col_label, folder) {
  mod_path <- here(
    "Data/Environmental_predictors/PolarRes26/Regridded",
    folder, "comparison",
    "Mean_Annual_Temperature_FUTURE_2071_2100_ICEFREE.tif"
  )
  mod  <- rast(mod_path)
  mod  <- crop_mask_peninsula(mod, peninsula_vect)
  diff <- mod - align_to_model(CHELSA_future, mod)
  diff <- crop_mask_peninsula(diff, peninsula_vect)
  
  list(
    row_group = row_group,
    col_label = col_label,
    temp_10km = aggregate(mod,  fact = agg_fact, fun = "mean", na.rm = TRUE),
    diff_10km = aggregate(diff, fact = agg_fact, fun = "mean", na.rm = TRUE)
  )
}

model_data <- pmap(model_config, load_future_model)
names(model_data) <- paste(model_config$row_group, model_config$col_label, sep = " - ")

CHELSA_future_10km <- aggregate(CHELSA_future, fact = agg_fact, fun = "mean", na.rm = TRUE)

# ============================================================
# 3. SHARED COLOUR SCALES — now computed over Peninsula-only values,
#    so the ramp isn't stretched across continent-wide extremes
# ============================================================

all_temp_vals <- c(
  values(CHELSA_future_10km, na.rm = TRUE),
  unlist(lapply(model_data, function(m) values(m$temp_10km, na.rm = TRUE)))
)
temp_min    <- floor(min(all_temp_vals))
temp_max    <- ceiling(max(all_temp_vals))
temp_values <- make_scale_values(temp_min, temp_max)

all_diff_vals <- unlist(lapply(model_data, function(m) values(m$diff_10km, na.rm = TRUE)))
diff_min    <- floor(min(all_diff_vals, na.rm = TRUE))
diff_max    <- ceiling(max(all_diff_vals, na.rm = TRUE))
diff_breaks <- pretty(c(diff_min, diff_max), n = 6)
diff_values <- make_scale_values(diff_min, diff_max)

message(
  "Peninsula future bias range: ", round(min(all_diff_vals), 2), " to ",
  round(max(all_diff_vals), 2), " \u00B0C  ->  scale limits: ",
  diff_min, " to ", diff_max, " \u00B0C"
)

# ============================================================
# 4. BUILD PANELS
# ============================================================

p_chelsa <- make_temp_panel(CHELSA_future_10km, "CHELSA", temp_min, temp_max, temp_values)

model_panels <- lapply(model_data, function(m) {
  list(
    temp = make_temp_panel(m$temp_10km, m$col_label, temp_min, temp_max, temp_values),
    diff = make_diff_panel(m$diff_10km, diff_min, diff_max, diff_breaks, diff_values)
  )
})

# ============================================================
# 5. ASSEMBLE — per-panel legends (no guides = "collect"/"keep"
#    anywhere), matching the previous per-panel-legend version
# ============================================================

build_model_row <- function(panels) {
  panels$temp + plot_spacer() + panels$diff +
    plot_layout(ncol = 3, widths = c(1, 0.015, 1))
}

build_storyline_block <- function(storyline_label) {
  rows <- model_config %>% filter(row_group == storyline_label) %>% pull(col_label)
  keys <- paste(storyline_label, rows, sep = " - ")
  
  row_plots <- lapply(model_panels[keys], build_model_row)
  
  wrap_plots(row_plots, ncol = 1) +
    plot_annotation(title = storyline_label, theme = row_title_theme)
}

chelsa_block <- plot_spacer() + p_chelsa + plot_spacer() +
  plot_layout(ncol = 3, widths = c(0.06, 1, 0.06))

storyline_blocks <- lapply(storyline_levels, build_storyline_block)

final_plot <- (chelsa_block / storyline_blocks[[1]] / storyline_blocks[[2]]) +
  plot_layout(heights = c(1, 2, 2))

# ============================================================
# 6. SAVE — different filename so the full-continent version is
#    preserved
# ============================================================

out_file <- file.path(outpath, "plot_annual_temp_FUTURE_comparison_CHELSA_vs_HCLIM_RACMO_PENINSULA.png")

ggsave(out_file, final_plot, width = PAGE_W, height = PAGE_H,
       dpi = DPI, bg = "white")

message("Saved: ", out_file)