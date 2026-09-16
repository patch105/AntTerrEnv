

# HISTORICAL TO FUTURE CLIMATE SUMMARIES ----------------------------------

# This code summarises changes in climate variables from historical (1995-2014) to future (2081-2100) conditions for our climate scenarios

library(dplyr)
library(purrr)
library(terra)
library(here)
library(RColorBrewer)
library(ggplot2)
library(sf)
library(ggpubr)


# Set Climate Scenarios ----------------------------------------------------

climate_scenarios <- list("MAR_MPI_ESM", "MAR_CESM2", "HCLIM_MPI_ESM1", "HCLIM_CESM2")


# Load Peninsula to Continent boundary ------------------------------------

# Load the boundary
pen_cont_boundary <- vect(here("Data/Environmental_predictors/Peninsula_Continent_Boundary.shp"))


# Read historical and future variables ------------------------------------

climate_dfs_list <-  imap(climate_scenarios, function(climate_source, name) {
  
  climate_dfs_list <- list()
  
  annual_temp_hist <- rast(here("Data/Environmental_predictors", climate_source, "Mean_Annual_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"))
  annual_temp_future <- rast(here("Data/Environmental_predictors", climate_source,"Mean_Annual_Temperature_FUTURE_2081_2100_ICEFREE.tif"))
  
  summer_temp_hist <- rast(here("Data/Environmental_predictors", climate_source, "Mean_Summer_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"))
  summer_temp_future <- rast(here("Data/Environmental_predictors", climate_source,"Mean_Summer_Temperature_FUTURE_2081_2100_ICEFREE.tif"))
  
  winter_temp_hist <- rast(here("Data/Environmental_predictors", climate_source, "Mean_Winter_Temperature_HISTORICAL_1995_2014_ICEFREE.tif"))
  winter_temp_future <- rast(here("Data/Environmental_predictors", climate_source,"Mean_Winter_Temperature_FUTURE_2081_2100_ICEFREE.tif"))
  
  nov_sea_ice_100km_hist <- rast(here("Data/Environmental_predictors", climate_source ,"Mean_November_Sea_Ice_Concentration_100km_HISTORICAL_1995_2014_ICEFREE.tif")) 
  nov_sea_ice_100km_future <- rast(here("Data/Environmental_predictors", climate_source, "Mean_November_Sea_Ice_Concentration_100km_FUTURE_2081_2100_ICEFREE.tif")) 
  
  oct_feb_sea_ice_42km_hist <- rast(here("Data/Environmental_predictors", climate_source ,"Mean_Oct_Feb_Sea_Ice_Concentration_42km_HISTORICAL_1995_2014_ICEFREE.tif")) 
  oct_feb_sea_ice_42km_future <- rast(here("Data/Environmental_predictors", climate_source, "Mean_Oct_Feb_Sea_Ice_Concentration_42km_FUTURE_2081_2100_ICEFREE.tif")) 
  
  # Degree Days
  degree_days_hist <- rast(here("Data/Environmental_predictors", climate_source, 
                                "Mean_Annual_Total_Degree_Days-5_HISTORICAL_1995_2014_ICEFREE.tif"))
  degree_days_future <- rast(here("Data/Environmental_predictors", climate_source, 
                                  "Mean_Annual_Total_Degree_Days-5_FUTURE_2081_2100_ICEFREE.tif"))
  
  
  # Solar Radiation
  solar_rad_hist <- rast(here("Data/Environmental_predictors", climate_source, 
                              "Mean_Annual_Solar_Radiation_HISTORICAL_1995_2014_ICEFREE.tif"))
  solar_rad_future <- rast(here("Data/Environmental_predictors", climate_source, 
                                "Mean_Annual_Solar_Radiation_FUTURE_2081_2100_ICEFREE.tif"))
  
  
  # Summer Precipitation
  summer_precip_hist <- rast(here("Data/Environmental_predictors", climate_source, 
                                  "Mean_Summer_Total_Precipitation_HISTORICAL_1995_2014_ICEFREE.tif"))
  summer_precip_future <- rast(here("Data/Environmental_predictors", climate_source, 
                                    "Mean_Summer_Total_Precipitation_FUTURE_2081_2100_ICEFREE.tif"))
  
  
  # Annual Precipitation
  annual_precip_hist <- rast(here("Data/Environmental_predictors", climate_source, 
                                  "Total_Annual_Precipitation_HISTORICAL_1995_2014_ICEFREE.tif"))
  annual_precip_future <- rast(here("Data/Environmental_predictors", climate_source, 
                                    "Total_Annual_Precipitation_FUTURE_2081_2100_ICEFREE.tif"))
  
  # # Snow Area
  # snow_area_hist <- rast(here("Data/Environmental_predictors", climate_source, 
  #                             "Mean_Snow_Area_Percentage_HISTORICAL_1995_2014_ICEFREE.tif"))
  # snow_area_future <- rast(here("Data/Environmental_predictors", climate_source, 
  #                               "Mean_Snow_Area_Percentage_FUTURE_2081_2100_ICEFREE.tif"))
  
  # Wind Speed
  wind_speed_hist <- rast(here("Data/Environmental_predictors", climate_source, 
                               "Mean_Annual_Wind_Speed_HISTORICAL_1995_2014_ICEFREE.tif"))
  wind_speed_future <- rast(here("Data/Environmental_predictors", climate_source, 
                                 "Mean_Annual_Wind_Speed_FUTURE_2081_2100_ICEFREE.tif"))
  
  
  # Stack historical --------------------------------------------------------
  
  historical <- c(annual_temp_hist, summer_temp_hist, winter_temp_hist, nov_sea_ice_100km_hist, oct_feb_sea_ice_42km_hist, degree_days_hist, solar_rad_hist, summer_precip_hist, annual_precip_hist, wind_speed_hist)
  
  future <- c(annual_temp_future, summer_temp_future, winter_temp_future, nov_sea_ice_100km_future, oct_feb_sea_ice_42km_future, degree_days_future, solar_rad_future, summer_precip_future, annual_precip_future, wind_speed_future)
  
  cov_names <- c("AnnualTemp", "SummerTemp", "WinterTemp", "NovSeaIce100km", "OctFebSeaIce42km", "DegreeDays", "SolarRad", "SummerPrecip", "AnnualPrecip", "WindSpeed")
  
  names(historical) <- cov_names
  names(future) <- cov_names
  
  # Crop covariates to Peninsula / Continent --------------------------------
  
  peninsula_historical <- mask(historical, pen_cont_boundary)
  peninsula_historical <- crop(peninsula_historical, ext(pen_cont_boundary))
  peninsula_historical_df <- as.data.frame(peninsula_historical, xy = F, na.rm = T)
  
  peninsula_future <- mask(future, pen_cont_boundary)
  peninsula_future <- crop(peninsula_future, ext(pen_cont_boundary))
  peninsula_future_df <- as.data.frame(peninsula_future, xy = F, na.rm = T)
  
  peninsula_diff <- peninsula_future - peninsula_historical
  peninsula_diff_df <- as.data.frame(peninsula_diff, xy = F, na.rm = T)
  
  continent_historical <- mask(historical, pen_cont_boundary, inverse = T)
  continent_historical <- continent_historical %>% 
    setNames(cov_names)  
  continent_historical_df <- as.data.frame(continent_historical, xy = F, na.rm = T)
  
  continent_future <- mask(future, pen_cont_boundary, inverse = T)
  continent_future <- continent_future %>% 
    setNames(cov_names)  
  continent_future_df <- as.data.frame(continent_future, xy = F, na.rm = T)
  
  continent_diff <- continent_future - continent_historical
  continent_diff_df <- as.data.frame(continent_diff, xy = F, na.rm = T)
  
  
  result_list <- list(peninsula_historical_df = peninsula_historical_df,
                      peninsula_future_df = peninsula_future_df,
                      peninsula_diff_df = peninsula_diff_df,
                      continent_historical_df = continent_historical_df,
                      continent_future_df = continent_future_df,
                      continent_diff_df = continent_diff_df)
  
  climate_dfs_list[[name]] <- result_list
  
  
} )


##########################################
# MAKE A TABLE TO SUMMARISE THE DIFFERENCES -------------------------------
##########################################

outpath <- here("Outputs/RESULTS/CLIMATE_SCENARIO_SUMMARIES")

library(tidyr)

# 1. Summarise entire continent climate changes (DIFF) -------------------------

# Peninsula:

peninsula_diff_df1 <- climate_dfs_list[[1]]$peninsula_diff_df %>% 
  mutate(Scenario = "1") %>% 
  mutate(RCM = "1_MAR")

peninsula_diff_df2 <- climate_dfs_list[[2]]$peninsula_diff_df %>% 
  mutate(Scenario = "2") %>% 
  mutate(RCM = "2_MAR")

peninsula_diff_df3 <- climate_dfs_list[[3]]$peninsula_diff_df %>% 
  mutate(Scenario = "1") %>% 
  mutate(RCM = "1_HCLIM")

peninsula_diff_df4 <- climate_dfs_list[[4]]$peninsula_diff_df %>% 
  mutate(Scenario = "2") %>% 
  mutate(RCM = "2_HCLIM")

peninsula_diff_df <- rbind(peninsula_diff_df1, peninsula_diff_df2, peninsula_diff_df3, peninsula_diff_df4)


# Now continent:

continent_diff_df1 <- climate_dfs_list[[1]]$continent_diff_df %>% 
  mutate(Scenario = "1") %>% 
  mutate(RCM = "1_MAR")

continent_diff_df2 <- climate_dfs_list[[2]]$continent_diff_df %>% 
  mutate(Scenario = "2") %>% 
  mutate(RCM = "2_MAR")

continent_diff_df3 <- climate_dfs_list[[3]]$continent_diff_df %>% 
  mutate(Scenario = "1") %>% 
  mutate(RCM = "1_HCLIM")

continent_diff_df4 <- climate_dfs_list[[4]]$continent_diff_df %>% 
  mutate(Scenario = "2") %>% 
  mutate(RCM = "2_HCLIM")

continent_diff_df <- rbind(continent_diff_df1, continent_diff_df2, continent_diff_df3, continent_diff_df4)


# Combined version:

all_diff_df <- rbind(continent_diff_df, peninsula_diff_df)
names(all_diff_df) <- c("AnnualTemp", "SummerTemp", "WinterTemp", "NovSeaIce100km", "OctFebSeaIce42km", "DegreeDays", "SolarRad", "SummerPrecip", "AnnualPrecip", "WindSpeed", "Scenario", "RCM")


diff_summary_stats_Storylines <- all_diff_df %>%
  group_by(Scenario) %>%
  summarise(across(where(is.numeric), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE)),
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = -Scenario,
               names_to = c("variable", ".value"),
               names_sep = "_") %>%
  mutate(mean_sd = paste0(round(mean, 2), " +/- ", round(sd, 2))) %>%
  select(Scenario, variable, mean_sd) %>%
  pivot_wider(names_from = Scenario, values_from = mean_sd)

print(diff_summary_stats_Storylines)

diff_summary_stats_Scenarios <- all_diff_df %>%
  group_by(RCM) %>%
  summarise(across(where(is.numeric), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE)),
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = -RCM,
               names_to = c("variable", ".value"),
               names_sep = "_") %>%
  mutate(mean_sd = paste0(round(mean, 2), " +/- ", round(sd, 2))) %>%
  select(RCM, variable, mean_sd) %>%
  pivot_wider(names_from = RCM, values_from = mean_sd)

print(diff_summary_stats_Scenarios)


# 2. MEAN HISTORICAL  -------------------------

# Peninsula:

peninsula_hist_df1 <- climate_dfs_list[[1]]$peninsula_historical_df %>% 
  mutate(Scenario = "1") %>% 
  mutate(RCM = "1_MAR")

peninsula_hist_df2 <- climate_dfs_list[[2]]$peninsula_historical_df %>% 
  mutate(Scenario = "2") %>% 
  mutate(RCM = "2_MAR")

peninsula_hist_df3 <- climate_dfs_list[[3]]$peninsula_historical_df %>% 
  mutate(Scenario = "1") %>% 
  mutate(RCM = "1_HCLIM")

peninsula_hist_df4 <- climate_dfs_list[[4]]$peninsula_historical_df %>% 
  mutate(Scenario = "2") %>% 
  mutate(RCM = "2_HCLIM")

peninsula_hist_df <- rbind(peninsula_hist_df1, peninsula_hist_df2, peninsula_hist_df3, peninsula_hist_df4)


# Now continent:

continent_hist_df1 <- climate_dfs_list[[1]]$continent_historical_df %>% 
  mutate(Scenario = "1") %>% 
  mutate(RCM = "1_MAR")

continent_hist_df2 <- climate_dfs_list[[2]]$continent_historical_df %>% 
  mutate(Scenario = "2") %>% 
  mutate(RCM = "2_MAR")

continent_hist_df3 <- climate_dfs_list[[3]]$continent_historical_df %>% 
  mutate(Scenario = "1") %>% 
  mutate(RCM = "1_HCLIM")

continent_hist_df4 <- climate_dfs_list[[4]]$continent_historical_df %>% 
  mutate(Scenario = "2") %>% 
  mutate(RCM = "2_HCLIM")

continent_hist_df <- rbind(continent_hist_df1, continent_hist_df2, continent_hist_df3, continent_hist_df4)


# Combined version:

all_hist_df <- rbind(continent_hist_df, peninsula_hist_df)
names(all_hist_df) <- c("AnnualTemp", "SummerTemp", "WinterTemp", "NovSeaIce100km", "OctFebSeaIce42km", "DegreeDays", "SolarRad", "SummerPrecip", "AnnualPrecip", "WindSpeed", "Scenario", "RCM")


hist_summary_stats_Storylines <- all_hist_df %>%
  group_by(Scenario) %>%
  summarise(across(where(is.numeric), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE)),
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = -Scenario,
               names_to = c("variable", ".value"),
               names_sep = "_") %>%
  mutate(mean_sd = paste0(round(mean, 2), " +/- ", round(sd, 2))) %>%
  select(Scenario, variable, mean_sd) %>%
  pivot_wider(names_from = Scenario, values_from = mean_sd)

print(hist_summary_stats_Storylines)

hist_summary_stats_Scenarios <- all_hist_df %>%
  group_by(RCM) %>%
  summarise(across(where(is.numeric), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE)),
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = -RCM,
               names_to = c("variable", ".value"),
               names_sep = "_") %>%
  mutate(mean_sd = paste0(round(mean, 2), " +/- ", round(sd, 2))) %>%
  select(RCM, variable, mean_sd) %>%
  pivot_wider(names_from = RCM, values_from = mean_sd)

print(hist_summary_stats_Scenarios)


# 3. MEAN FUTURE  -------------------------

# Peninsula:

peninsula_fut_df1 <- climate_dfs_list[[1]]$peninsula_future_df %>% 
  mutate(Scenario = "1") %>% 
  mutate(RCM = "1_MAR")

peninsula_fut_df2 <- climate_dfs_list[[2]]$peninsula_future_df %>% 
  mutate(Scenario = "2") %>% 
  mutate(RCM = "2_MAR")

peninsula_fut_df3 <- climate_dfs_list[[3]]$peninsula_future_df %>% 
  mutate(Scenario = "1") %>% 
  mutate(RCM = "1_HCLIM")

peninsula_fut_df4 <- climate_dfs_list[[4]]$peninsula_future_df %>% 
  mutate(Scenario = "2") %>% 
  mutate(RCM = "2_HCLIM")

peninsula_fut_df <- rbind(peninsula_fut_df1, peninsula_fut_df2, peninsula_fut_df3, peninsula_fut_df4)


# Now continent:

continent_fut_df1 <- climate_dfs_list[[1]]$continent_future_df %>% 
  mutate(Scenario = "1") %>% 
  mutate(RCM = "1_MAR")

continent_fut_df2 <- climate_dfs_list[[2]]$continent_future_df %>% 
  mutate(Scenario = "2") %>% 
  mutate(RCM = "2_MAR")

continent_fut_df3 <- climate_dfs_list[[3]]$continent_future_df %>% 
  mutate(Scenario = "1") %>% 
  mutate(RCM = "1_HCLIM")

continent_fut_df4 <- climate_dfs_list[[4]]$continent_future_df %>% 
  mutate(Scenario = "2") %>% 
  mutate(RCM = "2_HCLIM")

continent_fut_df <- rbind(continent_fut_df1, continent_fut_df2, continent_fut_df3, continent_fut_df4)


# Combined version:

all_fut_df <- rbind(continent_fut_df, peninsula_fut_df)
names(all_fut_df) <- c("AnnualTemp", "SummerTemp", "WinterTemp", "NovSeaIce100km", "OctFebSeaIce42km", "DegreeDays", "SolarRad", "SummerPrecip", "AnnualPrecip", "WindSpeed", "Scenario", "RCM")


fut_summary_stats_Storylines <- all_fut_df %>%
  group_by(Scenario) %>%
  summarise(across(where(is.numeric), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE)),
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = -Scenario,
               names_to = c("variable", ".value"),
               names_sep = "_") %>%
  mutate(mean_sd = paste0(round(mean, 2), " +/- ", round(sd, 2))) %>%
  select(Scenario, variable, mean_sd) %>%
  pivot_wider(names_from = Scenario, values_from = mean_sd)

print(fut_summary_stats_Storylines)

fut_summary_stats_Scenarios <- all_fut_df %>%
  group_by(RCM) %>%
  summarise(across(where(is.numeric), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE)),
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = -RCM,
               names_to = c("variable", ".value"),
               names_sep = "_") %>%
  mutate(mean_sd = paste0(round(mean, 2), " +/- ", round(sd, 2))) %>%
  select(RCM, variable, mean_sd) %>%
  pivot_wider(names_from = RCM, values_from = mean_sd)

print(fut_summary_stats_Scenarios)


# 4. SAVE ALL SUMMARY TABLES -------------------------

CLIMATE_SCENARIO_SUMMARIES <- list(
  diff_summary_stats_Storylines  = diff_summary_stats_Storylines,
  diff_summary_stats_Scenarios   = diff_summary_stats_Scenarios,
  hist_summary_stats_Storylines  = hist_summary_stats_Storylines,
  hist_summary_stats_Scenarios   = hist_summary_stats_Scenarios,
  fut_summary_stats_Storylines   = fut_summary_stats_Storylines,
  fut_summary_stats_Scenarios    = fut_summary_stats_Scenarios
)


# 4. SAVE ALL SUMMARY TABLES -------------------------

write.csv(diff_summary_stats_Storylines, file.path(outpath, "diff_summary_stats_Storylines.csv"), row.names = FALSE)
write.csv(diff_summary_stats_Scenarios,  file.path(outpath, "diff_summary_stats_Scenarios.csv"),  row.names = FALSE)
write.csv(hist_summary_stats_Storylines, file.path(outpath, "hist_summary_stats_Storylines.csv"), row.names = FALSE)
write.csv(hist_summary_stats_Scenarios,  file.path(outpath, "hist_summary_stats_Scenarios.csv"),  row.names = FALSE)
write.csv(fut_summary_stats_Storylines,  file.path(outpath, "fut_summary_stats_Storylines.csv"),  row.names = FALSE)
write.csv(fut_summary_stats_Scenarios,   file.path(outpath, "fut_summary_stats_Scenarios.csv"),   row.names = FALSE)


# 5. SAVE THE SUMMARY DATAFRAMES ------------------------------------------

# 5. REPEAT BY REGION (PENINSULA & CONTINENT SEPARATELY) -------------------------

# Helper function to summarise by Storyline and Scenario for a given region df
summarise_region <- function(df, region_name) {
  
  by_storyline <- df %>%
    group_by(Scenario) %>%
    summarise(across(where(is.numeric),
                     list(mean = ~mean(., na.rm = TRUE),
                          sd   = ~sd(.,   na.rm = TRUE)),
                     .names = "{.col}_{.fn}")) %>%
    pivot_longer(cols = -Scenario,
                 names_to = c("variable", ".value"),
                 names_sep = "_") %>%
    mutate(mean_sd = paste0(round(mean, 2), " +/- ", round(sd, 2)),
           Region  = region_name) %>%
    select(Region, Scenario, variable, mean_sd) %>%
    pivot_wider(names_from = Scenario, values_from = mean_sd)
  
  by_scenario <- df %>%
    group_by(RCM) %>%
    summarise(across(where(is.numeric),
                     list(mean = ~mean(., na.rm = TRUE),
                          sd   = ~sd(.,   na.rm = TRUE)),
                     .names = "{.col}_{.fn}")) %>%
    pivot_longer(cols = -RCM,
                 names_to = c("variable", ".value"),
                 names_sep = "_") %>%
    mutate(mean_sd = paste0(round(mean, 2), " +/- ", round(sd, 2)),
           Region  = region_name) %>%
    select(Region, RCM, variable, mean_sd) %>%
    pivot_wider(names_from = RCM, values_from = mean_sd)
  
  list(Storylines = by_storyline, Scenarios = by_scenario)
}

# --- DIFF ---

peninsula_diff_df_named <- peninsula_diff_df
names(peninsula_diff_df_named) <- c("AnnualTemp", "SummerTemp", "WinterTemp", "NovSeaIce100km", "OctFebSeaIce42km", "DegreeDays", "SolarRad", "SummerPrecip", "AnnualPrecip", "WindSpeed", "Scenario", "RCM")

continent_diff_df_named <- continent_diff_df
names(continent_diff_df_named) <- c("AnnualTemp", "SummerTemp", "WinterTemp", "NovSeaIce100km", "OctFebSeaIce42km", "DegreeDays", "SolarRad", "SummerPrecip", "AnnualPrecip", "WindSpeed", "Scenario", "RCM")

diff_peninsula <- summarise_region(peninsula_diff_df_named, "Peninsula")
diff_continent <- summarise_region(continent_diff_df_named, "Continent")

diff_region_Storylines <- rbind(diff_peninsula$Storylines, diff_continent$Storylines)
diff_region_Scenarios  <- rbind(diff_peninsula$Scenarios,  diff_continent$Scenarios)

print(diff_region_Storylines)
print(diff_region_Scenarios)

write.csv(diff_region_Storylines, file.path(outpath, "diff_region_Storylines.csv"), row.names = FALSE)
write.csv(diff_region_Scenarios,  file.path(outpath, "diff_region_Scenarios.csv"),  row.names = FALSE)


# --- HISTORICAL ---

peninsula_hist_df_named <- peninsula_hist_df
names(peninsula_hist_df_named) <- c("AnnualTemp", "SummerTemp", "WinterTemp", "NovSeaIce100km", "OctFebSeaIce42km", "DegreeDays", "SolarRad", "SummerPrecip", "AnnualPrecip", "WindSpeed", "Scenario", "RCM")

continent_hist_df_named <- continent_hist_df
names(continent_hist_df_named) <- c("AnnualTemp", "SummerTemp", "WinterTemp", "NovSeaIce100km", "OctFebSeaIce42km", "DegreeDays", "SolarRad", "SummerPrecip", "AnnualPrecip", "WindSpeed", "Scenario", "RCM")

hist_peninsula <- summarise_region(peninsula_hist_df_named, "Peninsula")
hist_continent <- summarise_region(continent_hist_df_named, "Continent")

hist_region_Storylines <- rbind(hist_peninsula$Storylines, hist_continent$Storylines)
hist_region_Scenarios  <- rbind(hist_peninsula$Scenarios,  hist_continent$Scenarios)

print(hist_region_Storylines)
print(hist_region_Scenarios)

write.csv(hist_region_Storylines, file.path(outpath, "hist_region_Storylines.csv"), row.names = FALSE)
write.csv(hist_region_Scenarios,  file.path(outpath, "hist_region_Scenarios.csv"),  row.names = FALSE)


# --- FUTURE ---

peninsula_fut_df_named <- peninsula_fut_df
names(peninsula_fut_df_named) <- c("AnnualTemp", "SummerTemp", "WinterTemp", "NovSeaIce100km", "OctFebSeaIce42km", "DegreeDays", "SolarRad", "SummerPrecip", "AnnualPrecip", "WindSpeed", "Scenario", "RCM")

continent_fut_df_named <- continent_fut_df
names(continent_fut_df_named) <- c("AnnualTemp", "SummerTemp", "WinterTemp", "NovSeaIce100km", "OctFebSeaIce42km", "DegreeDays", "SolarRad", "SummerPrecip", "AnnualPrecip", "WindSpeed", "Scenario", "RCM")

fut_peninsula <- summarise_region(peninsula_fut_df_named, "Peninsula")
fut_continent <- summarise_region(continent_fut_df_named, "Continent")

fut_region_Storylines <- rbind(fut_peninsula$Storylines, fut_continent$Storylines)
fut_region_Scenarios  <- rbind(fut_peninsula$Scenarios,  fut_continent$Scenarios)

print(fut_region_Storylines)
print(fut_region_Scenarios)

write.csv(fut_region_Storylines, file.path(outpath, "fut_region_Storylines.csv"), row.names = FALSE)
write.csv(fut_region_Scenarios,  file.path(outpath, "fut_region_Scenarios.csv"),  row.names = FALSE)



# 6. PLOT: STANDARDISED FUTURE VALUES BY REGION & STORYLINE -------------------------

library(ggplot2)
library(tidyr)
library(dplyr)

# Add Region column to peninsula and continent future dfs, then combine
peninsula_fut_plot_df <- peninsula_fut_df_named %>% mutate(Region = "Peninsula")
continent_fut_plot_df <- continent_fut_df_named %>% mutate(Region = "Continent")

all_fut_region_df <- rbind(peninsula_fut_plot_df, continent_fut_plot_df)

# Pivot to long format
all_fut_long <- all_fut_region_df %>%
  pivot_longer(cols = c(AnnualTemp, SummerTemp, WinterTemp, NovSeaIce100km,
                        OctFebSeaIce42km, DegreeDays, SolarRad,
                        SummerPrecip, AnnualPrecip, WindSpeed),
               names_to  = "Variable",
               values_to = "Value")

# Standardise each variable to [-1, 1] using global min/max across both regions
# so both panels share the same scale
all_fut_long <- all_fut_long %>%
  group_by(Variable) %>%
  mutate(Value_std = (Value - min(Value, na.rm = TRUE)) /
           (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE)) * 2 - 1) %>%
  ungroup()

# Variable order and labels to match reference plot
var_levels <- c("AnnualTemp", "SummerTemp", "WinterTemp",
                "AnnualPrecip", "SummerPrecip",
                "OctFebSeaIce42km", "NovSeaIce100km",
                "DegreeDays", "SolarRad", "WindSpeed")

var_labels <- c("Annual temperature", "Summer temperature", "Winter temperature",
                "Annual precipitation", "Summer precipitation",
                "Oct-Feb sea ice (42 km)", "Nov sea ice (100 km)",
                "Degree days -5C", "Solar radiation", "Wind speed")

all_fut_long <- all_fut_long %>%
  filter(Variable %in% var_levels) %>%
  mutate(
    Variable = factor(Variable, levels = rev(var_levels), labels = rev(var_labels)),
    Scenario = factor(Scenario, levels = c("2", "1")),
    Region   = factor(Region, levels = c("Peninsula", "Continent"))
  )

# Storyline colours matching reference (purple = 1, blue = 2)
storyline_colours <- c("1" = "#7B2D8B", "2" = "#4575B4")

p_fut_std <- ggplot(all_fut_long, aes(x = Value_std, y = Variable, fill = Scenario)) +
  geom_boxplot(
    aes(group = interaction(Variable, Scenario)),
    position     = position_dodge(width = 0.6),
    width        = 0.5,
    alpha        = 0.7,
    outlier.size = 0.8
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(
    values = storyline_colours,
    name   = "Climate Storyline",
    labels = c("1", "2")
  ) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  facet_wrap(~ Region, ncol = 2) +
  labs(
    x = "Standardised future value (-1 to 1)",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_line(colour = "grey92"),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(face = "bold", hjust = 0),
    strip.background   = element_blank(),
    legend.position    = "bottom",
    legend.title       = element_text(size = 10),
    axis.text.y        = element_text(size = 9),
    axis.text.x        = element_text(size = 9)
  )


ggsave(
  file.path(outpath, "fut_standardised_by_region_storyline.png"),
  plot   = p_fut_std,
  width  = 10,
  height = 6,
  dpi    = 300
)

