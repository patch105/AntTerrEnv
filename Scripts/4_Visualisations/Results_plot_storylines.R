library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)   # for facetted_pos_scales (per-panel x-axis control)
library(here)
library(svglite)


# ------------------------------------------------------------------
# 0. SETTINGS
# ------------------------------------------------------------------

# Final variable set + display order (left -> right, top -> bottom)
subset_vars <- c(
  "AnnualTemp", "SummerTemp", "WinterTemp", "DegreeDays0",
  "TotalAnnualPrecip", "WindSpeed", "SeaIceConc", "RelHumidity"
)

# Nice facet labels (edit to match your variable naming / units)
var_labels <- c(
  AnnualTemp        = "Annual temp (\u00B0C)",
  SummerTemp        = "Summer temp (\u00B0C)",
  WinterTemp        = "Winter temp (\u00B0C)",
  DegreeDays0       = "Degree days <0",
  TotalAnnualPrecip = "Total annual precip (mm)",
  WindSpeed         = "Wind speed (m/s)",
  SeaIceConc        = "Mean sea ice conc. (%)",
  RelHumidity       = "Relative humidity (%)"
)

driver_colours <- c(CESM2 = "#2a78d6", MPI_ESM1 = "#eb6834")

# Paler tints of the driver colours, used for the individual-RCM ticks so
# each tick reads as "belonging to" its driver without a separate legend
# entry (matches the softer orange/blue used in the annotated mock-up)
pale_driver_colours <- c(CESM2 = "#a9c9ee", MPI_ESM1 = "#f3b79a")

# Neutral grey for the driver-mean-to-driver-mean bar and the full-width
# baseline. Both share one linewidth (axis_linewidth below) so the
# "shaded" (data range) and "unshaded" (rest of the axis) portions of the
# horizontal line read as one continuous track, just different colours.
rcm_bar_colour   <- "grey60"
baseline_colour  <- "grey88"
axis_linewidth   <- 3
tick_linewidth   <- 1.1

# Point + text styling
point_size    <- 4.5
point_stroke  <- 1.1
title_colour  <- "grey45"   # panel/scope titles - light grey, slightly
axis_num_colour <- "grey60" # darker than the axis numbers
label_family <- "sans"      # swap for a specific installed font if you want
# an exact match, e.g. "Helvetica Neue" / "Roboto"

# Variables whose x-axis should be forced to a shared range so the
# panels can be visually compared to one another
shared_axis_group <- c("AnnualTemp", "SummerTemp", "WinterTemp")

# Axis-minimum "floors" -- applied after the automatic 0-floor logic
# below. The axis will extend further than this if the data needs it
# (e.g. an RCM below -0.1), but will never be *less* extreme than this
# even if all the data happens to sit above it.
axis_min_floor <- c(WindSpeed = -0.1)

# ------------------------------------------------------------------
# 1. DATA PREP FUNCTION (per scope)
# ------------------------------------------------------------------

prepare_storyline_data <- function(model_df, subset_vars) {
  
  present_vars <- intersect(subset_vars, unique(model_df$variable))
  missing_vars <- setdiff(subset_vars, present_vars)
  
  model_means <- model_df %>%
    filter(variable %in% present_vars) %>%
    pivot_longer(
      cols = c(HCLIM_CESM2, HCLIM_MPI_ESM1, RACMO_CESM2, RACMO_MPI_ESM1),
      names_to = "Model",
      values_to = "Value"
    ) %>%
    mutate(
      Value  = as.numeric(sub("\\s*\\+/-.*$", "", Value)),
      RCM    = case_when(grepl("^HCLIM", Model) ~ "HCLIM",
                         grepl("^RACMO", Model) ~ "RACMO"),
      Driver = case_when(grepl("CESM2$", Model)    ~ "CESM2",
                         grepl("MPI_ESM1$", Model) ~ "MPI_ESM1")
    ) %>%
    select(Model, Driver, RCM, variable, Value) %>%
    rename(Variable = variable)
  
  # Add placeholder rows for any requested variable not yet in the data
  # (e.g. SeaIceConc) so the facet panel still reserves its slot.
  if (length(missing_vars) > 0) {
    placeholder <- expand.grid(
      Model = NA_character_, Driver = c("CESM2", "MPI_ESM1"),
      RCM = NA_character_, Variable = missing_vars,
      stringsAsFactors = FALSE
    ) %>% mutate(Value = NA_real_)
    model_means <- bind_rows(model_means, placeholder)
  }
  
  model_means <- model_means %>%
    mutate(Variable = factor(Variable, levels = subset_vars))
  
  storyline_means <- model_means %>%
    group_by(Driver, Variable) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(Variable = factor(Variable, levels = subset_vars))
  
  segment_df <- storyline_means %>%
    pivot_wider(names_from = Driver, values_from = Value) %>%
    mutate(Variable = factor(Variable, levels = subset_vars))
  
  list(
    model_means     = model_means,
    storyline_means = storyline_means,
    segment_df      = segment_df,
    missing_vars    = missing_vars
  )
}

# ------------------------------------------------------------------
# 2. PER-VARIABLE X-AXIS SCALES (0 fixed at left, single max break,
#    shared range for the temperature trio, manual overrides supported)
# ------------------------------------------------------------------

build_x_scales <- function(model_means, subset_vars, shared_axis_group,
                           axis_min_floor = c()) {
  
  axis_range <- model_means %>%
    group_by(Variable) %>%
    summarise(min_val = min(Value, na.rm = TRUE),
              max_val = max(Value, na.rm = TRUE))
  
  axis_min <- tibble::deframe(axis_range %>% select(Variable, min_val))
  axis_max <- tibble::deframe(axis_range %>% select(Variable, max_val))
  
  # replace any +/-Inf (all-NA placeholder variables) with nominal values
  axis_min[!is.finite(axis_min)] <- 0
  axis_max[!is.finite(axis_max)] <- 1
  
  # 0 is always kept as the *lower bound* unless the data actually goes
  # negative (e.g. WindSpeed), in which case drop the lower limit to fit it
  axis_min <- pmin(axis_min, 0)
  
  # apply any floors (e.g. WindSpeed -> at least -0.1) -- takes whichever
  # is more extreme so a lower RCM value still gets fully captured
  for (v in names(axis_min_floor)) {
    if (v %in% names(axis_min)) {
      axis_min[[v]] <- min(axis_min[[v]], axis_min_floor[[v]])
    }
  }
  
  shared_max <- max(axis_max[shared_axis_group], na.rm = TRUE)
  axis_max[shared_axis_group] <- shared_max
  
  range_df <- tibble::tibble(
    Variable = factor(subset_vars, levels = subset_vars),
    xmin     = axis_min[subset_vars],
    xmax     = axis_max[subset_vars]
  )
  
  # Build the scale list IN THE SAME ORDER as subset_vars/facet levels
  scales <- lapply(subset_vars, function(v) {
    lower <- axis_min[[v]]
    upper <- axis_max[[v]]
    scale_x_continuous(
      limits = c(lower, upper),
      breaks = c(round(lower, 1), round(upper, 1)),
      expand = expansion(mult = c(0.03, 0.08))
    )
  })
  
  list(scales = scales, range_df = range_df)
}

# ------------------------------------------------------------------
# 3. PLOTTING FUNCTION
# ------------------------------------------------------------------

make_storyline_plot <- function(model_df, subset_vars, var_labels,
                                driver_colours, shared_axis_group,
                                axis_min_floor = c(), scope_title = NULL) {
  
  d  <- prepare_storyline_data(model_df, subset_vars)
  xs <- build_x_scales(d$model_means, subset_vars, shared_axis_group,
                       axis_min_floor)
  
  p <- ggplot() +
    
    # pale full-width baseline ("unshaded" track), one per panel, spanning
    # that panel's axis -- same linewidth as the shaded bar below so the
    # two read as one continuous line
    geom_segment(data = xs$range_df,
                 aes(x = xmin, xend = xmax, y = 1, yend = 1),
                 colour = baseline_colour, linewidth = axis_linewidth,
                 lineend = "round") +
    
    # grey bar ("shaded" portion) joining the two driver (storyline) means
    geom_segment(data = d$segment_df,
                 aes(x = CESM2, xend = MPI_ESM1, y = 1, yend = 1),
                 colour = rcm_bar_colour, linewidth = axis_linewidth,
                 lineend = "round") +
    
    # individual RCM marks as short vertical ticks, coloured by driver
    # (pale blue for CESM2-driven runs, pale orange for MPI_ESM1-driven)
    geom_segment(data = d$model_means,
                 aes(x = Value, xend = Value, y = 0.9, yend = 1.1,
                     colour = Driver),
                 linewidth = tick_linewidth, na.rm = TRUE, show.legend = FALSE) +
    
    # driver (storyline) means as larger dots with a white outline
    geom_point(data = d$storyline_means,
               aes(x = Value, y = 1, fill = Driver),
               shape = 21, colour = "white", stroke = point_stroke,
               size = point_size, na.rm = TRUE) +
    
    facet_wrap2(~ Variable, scales = "free_x", ncol = 4, drop = FALSE,
                labeller = as_labeller(var_labels)) +
    facetted_pos_scales(x = xs$scales) +
    
    scale_colour_manual(values = pale_driver_colours, guide = "none") +
    scale_fill_manual(values = driver_colours, name = NULL) +
    guides(fill = guide_legend(override.aes = list(size = point_size))) +
    
    scale_y_continuous(limits = c(0.5, 1.5), breaks = NULL) +
    labs(x = NULL, y = NULL, title = scope_title) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid           = element_blank(),
      strip.text           = element_text(face = "plain", hjust = 0,
                                          colour = title_colour,
                                          family = label_family,
                                          size = rel(1.05),
                                          margin = margin(b = 3)),
      strip.placement       = "outside",
      axis.text.x           = element_text(colour = axis_num_colour,
                                           family = label_family,
                                           size = rel(0.85)),
      axis.ticks.x          = element_blank(),
      axis.line.x           = element_blank(),
      panel.spacing.x       = unit(2.2, "lines"),
      panel.spacing.y       = unit(2, "lines"),
      legend.position        = "top",
      legend.justification   = "left",
      legend.text             = element_text(family = label_family, size = rel(0.9)),
      plot.title               = element_text(face = "plain", colour = title_colour,
                                              family = label_family),
      plot.title.position      = "plot"
    )
  
  # Placeholder label for any variable not yet in the data (e.g. sea ice)
  if (length(d$missing_vars) > 0) {
    placeholder_range <- xs$range_df %>% filter(Variable %in% d$missing_vars)
    placeholder_df <- placeholder_range %>%
      mutate(x = (xmin + xmax) / 2, y = 1)
    p <- p + geom_text(data = placeholder_df,
                       aes(x = x, y = y, label = "Data pending"),
                       colour = "grey60", size = 3, fontface = "italic")
  }
  
  p
}

# ------------------------------------------------------------------
# 4. RUN FOR THE THREE SCOPES
# ------------------------------------------------------------------

all_model_df       <- read.csv(here("Data/Environmental_predictors/PolarRes26/Scenario_summaries/FUTURE_DIFF_All_by_model.csv"))
peninsula_model_df  <- read.csv(here("Data/Environmental_predictors/PolarRes26/Scenario_summaries/FUTURE_DIFF_PENINSULA_by_model.csv"))
continent_model_df  <- read.csv(here("Data/Environmental_predictors/PolarRes26/Scenario_summaries/FUTURE_DIFF_CONTINENT_by_model.csv"))

# Each scope gets its own axis ranges (computed from its own data only).
# Within a given scope's plot, all 8 panels are still the same physical
# width regardless of their individual data ranges -- that's inherent to
# facet_wrap's grid layout, not something the axis ranges affect.
plot_all <- make_storyline_plot(all_model_df, subset_vars, var_labels,
                                driver_colours, shared_axis_group,
                                axis_min_floor, scope_title = "Ice-free land")

plot_peninsula <- make_storyline_plot(peninsula_model_df, subset_vars, var_labels,
                                      driver_colours, shared_axis_group,
                                      axis_min_floor, scope_title = "Peninsula")

plot_continent <- make_storyline_plot(continent_model_df, subset_vars, var_labels,
                                      driver_colours, shared_axis_group,
                                      axis_min_floor, scope_title = "Continent")

# ------------------------------------------------------------------
# 5. EXPORT AS TRUE VECTOR SVG (identical size => panels stay aligned
#    when the three files are combined in Inkscape)
# ------------------------------------------------------------------

panel_width  <- 10   # inches
panel_height <- 5    # inches

ggsave(here("Plots/storyline_all.svg"),       plot_all,       device = "svg",
       width = panel_width, height = panel_height)
ggsave(here("Plots/storyline_peninsula.svg"), plot_peninsula, device = "svg",
       width = panel_width, height = panel_height)
ggsave(here("Plots/storyline_continent.svg"), plot_continent, device = "svg",
       width = panel_width, height = panel_height)