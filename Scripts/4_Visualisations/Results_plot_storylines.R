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

# Variables whose x-axis should be forced to a shared range so the
# panels can be visually compared to one another
shared_axis_group <- c("AnnualTemp", "SummerTemp", "WinterTemp")

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
    mutate(
      Variable = factor(Variable, levels = subset_vars),
      Mid      = (CESM2 + MPI_ESM1) / 2
    )
  
  list(
    model_means     = model_means,
    storyline_means = storyline_means,
    segment_df      = segment_df,
    missing_vars    = missing_vars
  )
}

# ------------------------------------------------------------------
# 2. PER-VARIABLE X-AXIS SCALES (0 fixed at left, single max break,
#    shared range for the temperature trio)
# ------------------------------------------------------------------

build_x_scales <- function(model_means, subset_vars, shared_axis_group) {
  
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
  
  shared_max <- max(axis_max[shared_axis_group], na.rm = TRUE)
  axis_max[shared_axis_group] <- shared_max
  
  # Build the scale list IN THE SAME ORDER as subset_vars/facet levels
  lapply(subset_vars, function(v) {
    lower <- axis_min[[v]]
    upper <- axis_max[[v]]
    scale_x_continuous(
      limits = c(lower, upper),
      breaks = c(round(lower, 1), round(upper, 1)),
      expand = expansion(mult = c(0.03, 0.08))
    )
  })
}

# ------------------------------------------------------------------
# 3. PLOTTING FUNCTION
# ------------------------------------------------------------------

make_storyline_plot <- function(model_df, subset_vars, var_labels,
                                driver_colours, shared_axis_group,
                                scope_title = NULL) {
  
  d <- prepare_storyline_data(model_df, subset_vars)
  x_scales <- build_x_scales(d$model_means, subset_vars, shared_axis_group)
  
  p <- ggplot() +
    
    # individual RCM marks, coloured by driver
    geom_point(data = d$model_means,
               aes(x = Value, y = 1, fill = Driver, colour = Driver),
               shape = 21, size = 1, stroke = 0.3) +
    
    # storyline (driver) means
    geom_point(data = d$storyline_means,
               aes(x = Value, y = 1, colour = Driver, fill = Driver),
               size = 1.8) +
    
    facet_wrap2(~ Variable, scales = "free_x", ncol = 4, drop = FALSE,
                labeller = as_labeller(var_labels)) +
    facetted_pos_scales(x = x_scales) +
    
    scale_colour_manual(values = driver_colours, guide = "none") +
    scale_fill_manual(values = driver_colours, guide = "none") +
    
    scale_y_continuous(limits = c(0.5, 1.5), breaks = NULL) +
    labs(x = NULL, y = NULL, title = scope_title) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor    = element_blank(),
      strip.text          = element_text(face = "bold"),
      plot.title           = element_text(face = "bold")
    )
  
  # Placeholder label for any variable not yet in the data (e.g. sea ice)
  if (length(d$missing_vars) > 0) {
    placeholder_df <- data.frame(
      Variable = factor(d$missing_vars, levels = subset_vars),
      x = 0.5, y = 1
    )
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

plot_all <- make_storyline_plot(all_model_df, subset_vars, var_labels,
                                driver_colours, shared_axis_group,
                                scope_title = "Ice-free land")

plot_peninsula <- make_storyline_plot(peninsula_model_df, subset_vars, var_labels,
                                      driver_colours, shared_axis_group,
                                      scope_title = "Peninsula")

plot_continent <- make_storyline_plot(continent_model_df, subset_vars, var_labels,
                                      driver_colours, shared_axis_group,
                                      scope_title = "Continent")

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

