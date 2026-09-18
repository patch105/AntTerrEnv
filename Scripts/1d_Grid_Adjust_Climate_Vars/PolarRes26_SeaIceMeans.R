# ==============================================================================
# PolarRes26 -- STEP 4: ANNUAL / DJF / JJA MEANS OF SEA-ICE CONCENTRATION & BUFFER
# ==============================================================================
# Takes every monthly sea-ice output Script 3 wrote --
#   <stem>_CONCENTRATION.tif
#   <stem>_BUFFER_<radius>km.tif
# -- and, for each (model, period, product) combination, averages the 12
# monthly files into three climatological products on the same grid:
#   - ..._Annual_Mean_..._<product>.tif  -- mean of all 12 months
#   - ..._DJF_Mean_...<product>.tif      -- mean of Dec, Jan, Feb
#   - ..._JJA_Mean_...<product>.tif      -- mean of Jun, Jul, Aug
#
# Input filenames (as written by Script 3) look like:
#   Climatological_Monthly_Mean_Sea_Ice_Concentration_October_HISTORICAL_1995_2014_BUFFER_10km.tif
#   Climatological_Monthly_Mean_Sea_Ice_Concentration_January_FUTURE_2081_2100_CONCENTRATION.tif
#
# A "group" is one (model, period, year1, year2, product) combination --
# i.e. one physical quantity, one period, on one grid. Each group must have
# EXACTLY one file per calendar month (12 total); groups with a missing or
# duplicated month are reported and skipped rather than silently averaged
# over an incomplete year.
#
# Deliberately does NOT hard-code which periods/products exist -- it parses
# whatever monthly files Script 3 actually produced. One job = one group
# (i.e. one output triple of Annual/DJF/JJA for one model+period+product),
# selected via a single job_index (same pattern as Scripts 2 and 3) -- run
# with no argument (or an out-of-range one) to print the full table of
# every group that will be processed, and how many jobs that is.
# ==============================================================================

# ---- 0. Setup ------------------------------------------------------------------

lib_loc <- paste(getwd(), "/r_lib_new", sep = "")

library(dplyr, lib.loc = lib_loc)
library(tidyr, lib.loc = lib_loc)
library(purrr, lib.loc = lib_loc)
library(terra)
library(here)

# ---- 1. Configuration -----------------------------------------------------------

models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5", "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5", "MetUM_ERA5")

# Script 3's output directory is this script's input directory. Outputs
# from this script are written alongside them, in the same per-model folder.
regrid_base <- here("Data/Environmental_predictors/PolarRes26/Regridded")

month_levels <- c("January", "February", "March", "April", "May", "June",
                  "July", "August", "September", "October", "November", "December")
djf_months <- c("December", "January", "February")
jja_months <- c("June", "July", "August")

# Matches exactly the filenames Script 3 writes for sea-ice products:
#   Climatological_Monthly_Mean_Sea_Ice_Concentration_<Month>_<PERIOD>_<y1>_<y2>_<PRODUCT>.tif
sic_pattern <- paste0(
  "^Climatological_Monthly_Mean_Sea_Ice_Concentration_",
  "(", paste(month_levels, collapse = "|"), ")_",
  "(HISTORICAL|MID|FUTURE)_",
  "([0-9]{4})_([0-9]{4})_",
  "(CONCENTRATION|BUFFER_[0-9]+km)\\.tif$"
)

# If TRUE, skip a group entirely (no reading/averaging at all) when every
# output it would produce already exists on disk. Set to FALSE to always
# reprocess and overwrite.
SKIP_EXISTING <- TRUE

# ---- 2. Build the job table: every complete 12-month group, across every model --

file_table <- map_dfr(models, function(model) {
  model_dir <- file.path(regrid_base, model)
  if (!dir.exists(model_dir)) {
    message("Model directory not found, skipping: ", model_dir)
    return(tibble())
  }
  files <- list.files(model_dir, full.names = TRUE)
  matched <- grepl(sic_pattern, basename(files))
  if (!any(matched)) return(tibble())
  
  m <- regmatches(basename(files[matched]), regexec(sic_pattern, basename(files[matched])))
  tibble(
    model      = model,
    input_path = files[matched],
    month      = vapply(m, `[`, character(1), 2),
    period     = vapply(m, `[`, character(1), 3),
    year1      = vapply(m, `[`, character(1), 4),
    year2      = vapply(m, `[`, character(1), 5),
    product    = vapply(m, `[`, character(1), 6)
  )
})

if (nrow(file_table) == 0) {
  stop("No Script-3 sea-ice CONCENTRATION/BUFFER files found under ", regrid_base,
       " -- has Script 3 been run yet?")
}

grouped <- file_table %>%
  group_by(model, period, year1, year2, product) %>%
  summarise(
    n_months = n(),
    n_distinct_months = n_distinct(month),
    months = list(month),
    paths  = list(input_path),
    .groups = "drop"
  )

# A group is only usable if it has exactly 12 files and no repeated month.
grouped <- grouped %>%
  mutate(complete = n_months == 12 & n_distinct_months == 12)

incomplete <- grouped %>% filter(!complete)
if (nrow(incomplete) > 0) {
  message(nrow(incomplete), " group(s) skipped for having != 12 distinct monthly files:")
  for (i in seq_len(nrow(incomplete))) {
    g <- incomplete[i, ]
    message("  ", g$model, " / ", g$period, "_", g$year1, "_", g$year2, " / ", g$product,
            " -- found ", g$n_months, " file(s), ", g$n_distinct_months, " distinct month(s): ",
            paste(sort(g$months[[1]]), collapse = ", "))
  }
}

job_table <- grouped %>% filter(complete) %>% select(-complete)

if (nrow(job_table) == 0) {
  stop("No complete 12-month groups found -- nothing to average.")
}

message(nrow(job_table), " complete group(s) found across ",
        length(unique(job_table$model)), " model(s).")

# ---- 3. job_index selects ONE group (same pattern as Scripts 2 and 3) ----------

args <- commandArgs(trailingOnly = TRUE)
job_index <- suppressWarnings(as.integer(args[1]))

if (length(args) == 0 || is.na(job_index) || job_index < 1 || job_index > nrow(job_table)) {
  print(job_table %>% select(-months, -paths), n = Inf)
  message(nrow(job_table), " job(s) total. Pass a job_index between 1 and ",
          nrow(job_table), " to run a single job.")
  quit(save = "no", status = 0)
}

this_job <- job_table[job_index, ]
model    <- this_job$model
period   <- this_job$period
year1    <- this_job$year1
year2    <- this_job$year2
product  <- this_job$product

month_paths <- setNames(this_job$paths[[1]], this_job$months[[1]])

message("Job ", job_index, "/", nrow(job_table), " -> ", model, " / ", period, "_",
        year1, "_", year2, " / ", product)

# ---- 4. Output paths -------------------------------------------------------------

out_dir <- file.path(regrid_base, model)

make_stat_path <- function(stat_label) {
  file.path(out_dir, paste0(
    "Mean_", stat_label, "_Sea_Ice_Concentration_",
    period, "_", year1, "_", year2, "_", product, ".tif"
  ))
}

# DJF = austral summer, JJA = austral winter.
annual_path <- make_stat_path("Annual")
djf_path    <- make_stat_path("Summer")
jja_path    <- make_stat_path("Winter")

expected_outputs <- c(annual_path, djf_path, jja_path)

# ---- 4b. Skip-existing check ------------------------------------------------------

if (SKIP_EXISTING && all(file.exists(expected_outputs))) {
  message("  all expected output(s) already exist -- skipping: ",
          paste(basename(expected_outputs), collapse = ", "))
  quit(save = "no", status = 0)
}

# ---- 5. Load the 12 monthly rasters and average ----------------------------------
# All 12 files in a group are Script 3 outputs for the same model/period/
# product, so they already share one grid -- stack them directly rather
# than re-checking/re-aligning geometry here.

monthly_stack <- rast(unname(month_paths[month_levels]))
names(monthly_stack) <- month_levels

if (!SKIP_EXISTING || !file.exists(annual_path)) {
  r_annual <- mean(monthly_stack, na.rm = TRUE)
  writeRaster(r_annual, annual_path, overwrite = TRUE)
  message("  wrote ", annual_path)
}

if (!SKIP_EXISTING || !file.exists(djf_path)) {
  r_djf <- mean(monthly_stack[[djf_months]], na.rm = TRUE)
  writeRaster(r_djf, djf_path, overwrite = TRUE)
  message("  wrote ", djf_path)
}

if (!SKIP_EXISTING || !file.exists(jja_path)) {
  r_jja <- mean(monthly_stack[[jja_months]], na.rm = TRUE)
  writeRaster(r_jja, jja_path, overwrite = TRUE)
  message("  wrote ", jja_path)
}

message("Done: ", model, " / ", period, "_", year1, "_", year2, " / ", product)