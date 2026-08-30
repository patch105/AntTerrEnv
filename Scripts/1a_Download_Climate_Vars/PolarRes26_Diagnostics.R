# ==============================================================================
# PolarRes26 -- STEP 1: DIAGNOSTICS
# ==============================================================================
# Purpose
# -------
# Before writing the new climatology-summarising code for the PolarRes26
# netCDF downloads, this script characterises exactly how each
# model / variable / scenario is structured on disk, so Script 2 (the
# summarising script) can be written to match reality instead of assuming
# it still works the same way it did for the old PolarRes .tif workflow.
#
# It does NOT compute any climatologies. It only inspects files and reports.
#
# For every model in `models`, it walks every .nc file found under it and,
# per file, records:
#   - which VARIABLE folder it lives in (matched on the folder name itself,
#     not the filename -- since variables are now organised into their own
#     folders rather than just named in the filename)
#   - which SCENARIO folder it lives in (historical / sspNNN / evaluation --
#     the ERA5-driven models use "evaluation" instead of historical/sspNNN)
#   - which FREQUENCY folder it lives in (expected to always be "day")
#   - the number of time layers in the file, and the first/last timestamp
#     actually stored in the netCDF (read via terra + ncdf4, not guessed
#     from the filename)
#   - the nominal start/end implied by the filename's trailing date range
#   - whether the file's *actual* last date sits after its *nominal* end
#     date -- i.e. whether the old "file spills into the first day(s) of
#     the next chunk" idiosyncrasy is still happening
#   - the time step (days between layers) and whether it's uniform, so we
#     can confirm every variable really is daily now
#   - the calendar attribute on the time dimension (standard/noleap/360_day
#     etc. -- this matters for anything doing day-of-year or degree-day math)
#   - the units attribute stored in the netCDF itself
#
# It then rolls all of this up, per (model, variable, scenario), into:
#   - a coverage check against a full daily sequence (missing dates,
#     duplicated dates -- the same kind of problem as the old missing
#     rsus months)
#   - whether units/calendars are consistent across all files for that
#     variable, and whether they match what the OLD scripts assumed
#     (e.g. tas in Kelvin, pr as a kg m-2 s-1 flux, siconca as %)
#
# Outputs (written to Data/Diagnostics/PolarRes26/):
#   - PolarRes26_diagnostics_file_level.csv       (one row per .nc file)
#   - PolarRes26_diagnostics_variable_summary.csv (one row per model/var/scenario)
#   - PolarRes26_diagnostics_date_gaps.csv         (missing/duplicated dates)
#   - PolarRes26_diagnostics_units_by_variable.csv (units across ALL models)
# Plus a console summary flagging anything worth a human look.
# ==============================================================================

# ---- 0. Setup ----------------------------------------------------------------

# HPC
lib_loc <- paste(getwd(), "/r_lib_new", sep = "")
# .libPaths(lib_loc)  # uncomment if r_lib isn't already on the search path


library(terra)
library(ncdf4, lib.loc = lib_loc)
library(here)
library(lubridate)
library(stringr)
library(dplyr, lib.loc = lib_loc)
library(purrr, lib.loc = lib_loc)
library(tidyr, lib.loc = lib_loc)
library(readr)
library(tibble)

# ncdf4 is optional but strongly recommended: it lets us read the "units"
# and "calendar" attributes directly from the file rather than relying on
# terra to have parsed/exposed them. If it's missing, install it into your
# r_lib, e.g.: install.packages("ncdf4", lib = lib_loc)
has_ncdf4 <- requireNamespace("ncdf4", quietly = TRUE)



# ---- 1. Configuration ---------------------------------------------------------

base_dir <- "/mnt/hpccs01/home/patterc2/n11222026/AntarcticFutureHabitat/Data/PolarRes26"

models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5", "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5", "MetUM_ERA5")

# Expected variables per the new download script (includes the new 'hurs').
# If a model's actual folder names differ (e.g. RACMO uses "huss" instead of
# "hurs", or "uas"/"vas" instead of "sfcWind"), those files will show up
# under the "did not match any expected variable folder" warning below --
# check that output and update this vector if needed.
expected_variables <- c("tas", "tasmax", "tasmin", "sfcWind", "siconca",
                        "pr", "snc", "snm", "rsus", "rsds", "hurs")

# What the OLD (PolarRes/.tif) summarising scripts assumed each variable's
# units were, so we can flag if that assumption no longer holds.
expected_units <- c(
  tas      = "K",
  tasmax   = "K",
  tasmin   = "K",
  sfcWind  = "m s-1",
  siconca  = "%",
  pr       = "kg m-2 s-1",
  snc      = "%",
  snm      = "kg m-2 s-1",
  rsus     = "W m-2",
  rsds     = "W m-2",
  hurs     = "%"
)

years_hist   <- seq(1995, 2014, by = 1)
years_mid    <- seq(2041, 2060, by = 1)
years_future <- seq(2081, 2100, by = 1)

diag_outpath <- here("Data/PolarRes26/Diagnostics")
dir.create(diag_outpath, recursive = TRUE, showWarnings = FALSE)

# Allow running one model at a time on HPC (job array style, matching the
# old summarising scripts), or all four in one go if no argument is given.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && !is.na(suppressWarnings(as.integer(args[1])))) {
  models_to_run <- models[as.integer(args[1])]
} else {
  models_to_run <- models
}
message("Running diagnostics for: ", paste(models_to_run, collapse = ", "))

# ---- 2. Helper functions -------------------------------------------------------

# 2a. Find every .nc file for one model and work out, from the PATH itself,
#     which variable/scenario/frequency folder it lives under.
discover_files <- function(model_dir) {
  all_files <- list.files(model_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
  if (length(all_files) == 0) return(tibble())
  
  components_list <- str_split(all_files, "/")
  
  variable <- map_chr(components_list, function(cmp) {
    hit <- intersect(cmp, expected_variables)
    if (length(hit) == 1) hit else NA_character_
  })
  
  frequency <- map_chr(components_list, function(cmp) {
    idx <- which(cmp %in% c("day", "mon", "6hr", "3hr", "1hr", "yr"))
    if (length(idx) >= 1) cmp[idx[1]] else NA_character_
  })
  
  tibble(
    path      = all_files,
    filename  = basename(all_files),
    variable  = variable,
    frequency = frequency,
    # ERA5-driven models (HCLIM_ERA5, RACMO_ERA5, MetUM_ERA5) use an
    # "evaluation" folder in place of historical/sspNNN.
    scenario  = str_extract(all_files, "historical|ssp[0-9]+|evaluation")
  )
}

# 2b. Pull the nominal start/end out of the filename's trailing
#     "..._<start>-<end>.nc" token. Handles YYYY, YYYYMM and YYYYMMDD forms.
#     Vectorised over a character vector of filenames.
parse_nominal_range <- function(filenames) {
  m <- str_match(filenames, "_(\\d{4,8})-(\\d{4,8})\\.nc$")
  list(start_raw = m[, 2], end_raw = m[, 3])
}

to_date_vec <- function(x, end = FALSE) {
  # x: character vector of YYYY / YYYYMM / YYYYMMDD tokens (may contain NA)
  out <- as.Date(rep(NA_character_, length(x)))
  n <- nchar(x)
  
  idx4 <- which(n == 4)
  if (length(idx4)) {
    out[idx4] <- ymd(paste0(x[idx4], if (end) "1231" else "0101"))
  }
  idx6 <- which(n == 6)
  if (length(idx6)) {
    d <- ymd(paste0(x[idx6], "01"))
    out[idx6] <- if (end) (ceiling_date(d, "month") - days(1)) else d
  }
  idx8 <- which(n == 8)
  if (length(idx8)) {
    out[idx8] <- ymd(x[idx8])
  }
  out
}

# 2c. Read the metadata we need from one netCDF file, WITHOUT loading the
#     pixel data itself (terra::rast() is lazy, so this stays fast even for
#     large Antarctic-domain files).
read_file_meta <- function(path) {
  out <- list(
    nlyr = NA_integer_, first_time = as.POSIXct(NA), last_time = as.POSIXct(NA),
    step_days_median = NA_real_, step_uniform = NA, calendar = NA_character_,
    units = NA_character_, main_var = NA_character_, error = NA_character_
  )
  
  main_var <- NA_character_
  
  # Use ncdf4 (if available) to positively identify the main data variable,
  # its units, and the calendar -- this avoids terra having to guess which
  # variable/band to use in files that contain more than one variable
  # (e.g. bounds variables, a rotated_pole grid-mapping variable, etc.)
  if (has_ncdf4) {
    tryCatch({
      nc <- ncdf4::nc_open(path)
      on.exit(ncdf4::nc_close(nc), add = TRUE)
      candidates <- names(nc$var)
      candidates <- candidates[!grepl(
        "_bnds$|bounds|^lat$|^lon$|rotated|crs|height|^time$",
        candidates, ignore.case = TRUE
      )]
      if (length(candidates) >= 1) {
        main_var <- candidates[1]
        u <- ncdf4::ncatt_get(nc, main_var, "units")
        if (u$hasatt) out$units <- u$value
      }
      if (!is.null(nc$dim$time)) {
        cal <- ncdf4::ncatt_get(nc, "time", "calendar")
        if (cal$hasatt) out$calendar <- cal$value
      }
    }, error = function(e) NULL)
  }
  out$main_var <- main_var
  
  tryCatch({
    r <- if (!is.na(main_var)) terra::rast(path, subds = main_var) else terra::rast(path)
    out$nlyr <- terra::nlyr(r)
    
    tm <- terra::time(r)
    if (!is.null(tm) && length(tm) == out$nlyr && !all(is.na(tm))) {
      tm_sorted <- sort(tm)
      out$first_time <- min(tm, na.rm = TRUE)
      out$last_time  <- max(tm, na.rm = TRUE)
      if (length(tm_sorted) > 1) {
        steps <- as.numeric(diff(tm_sorted), units = "days")
        out$step_days_median <- median(steps)
        out$step_uniform <- (length(unique(round(steps, 4))) == 1)
      }
    }
    
    if (is.na(out$units)) {
      u <- try(terra::units(r), silent = TRUE)
      if (!inherits(u, "try-error") && length(u) > 0) out$units <- unique(u)[1]
    }
  }, error = function(e) out$error <<- conditionMessage(e))
  
  out
}

# ---- 3. Main loop: build the file-level table ----------------------------------

file_level <- list()

for (model in models_to_run) {
  
  model_dir <- file.path(base_dir, model)
  if (!dir.exists(model_dir)) {
    message("!! Model directory not found, skipping: ", model_dir)
    next
  }
  
  message("=== ", model, " ===")
  
  result <- tryCatch({
    
    found <- discover_files(model_dir)
    if (nrow(found) == 0) {
      message("  No .nc files found under ", model_dir)
      return(NULL)
    }
    
    # Files whose variable folder we couldn't match to expected_variables
    unmatched <- found %>% filter(is.na(variable))
    if (nrow(unmatched) > 0) {
      message("  !! ", nrow(unmatched), " file(s) did not match any expected ",
              "variable folder -- check these folder(s):")
      print(unique(dirname(unmatched$path)))
    }
    
    # Expected variables with no files at all
    found_vars <- unique(na.omit(found$variable))
    missing_vars <- setdiff(expected_variables, found_vars)
    if (length(missing_vars) > 0) {
      message("  !! Expected variable(s) with NO files found: ",
              paste(missing_vars, collapse = ", "))
    }
    
    # Frequency folders that aren't "day" (shouldn't happen per the download script)
    freqs <- unique(na.omit(found$frequency))
    if (length(setdiff(freqs, "day")) > 0) {
      message("  !! Non-daily frequency folder(s) present: ",
              paste(setdiff(freqs, "day"), collapse = ", "))
    }
    
    message("  ", nrow(found), " files across ", length(found_vars),
            " variable(s), ", length(unique(na.omit(found$scenario))), " scenario(s)")
    
    # Read metadata per file (the slow part -- progress every 50 files)
    meta_list <- vector("list", nrow(found))
    for (i in seq_len(nrow(found))) {
      if (i %% 50 == 0) message("    ...", i, "/", nrow(found))
      meta_list[[i]] <- read_file_meta(found$path[i])
    }
    meta_df <- bind_rows(lapply(meta_list, as_tibble))
    
    nominal <- parse_nominal_range(found$filename)
    
    combined <- bind_cols(found, meta_df) %>%
      mutate(
        model = model,
        nominal_start = to_date_vec(nominal$start_raw, end = FALSE),
        nominal_end   = to_date_vec(nominal$end_raw,   end = TRUE),
        first_date = as.Date(first_time),
        last_date  = as.Date(last_time),
        # The classic "spills into the next chunk" idiosyncrasy: does the
        # file's ACTUAL last date sit after its NOMINAL end date?
        spillover_suspected = !is.na(nominal_end) & !is.na(last_date) &
          (last_date > nominal_end)
      )
    
    combined
    
  }, error = function(e) {
    message("  !! Error processing model ", model, ": ", conditionMessage(e))
    NULL
  })
  
  if (!is.null(result)) file_level[[model]] <- result
}

file_level_df <- bind_rows(file_level)

if (nrow(file_level_df) == 0) {
  stop("No files were successfully processed for any model -- check base_dir and models.")
}

readr::write_csv(file_level_df, paste0(diag_outpath, "/PolarRes26_diagnostics_file_level_", models_to_run, ".csv"))
message("\nWrote file-level diagnostics: ",
        file.path(diag_outpath, models_to_run, "PolarRes26_diagnostics_file_level.csv"))

# ---- 4. Roll-up per (model, variable, scenario) --------------------------------

variable_summary <- file_level_df %>%
  filter(!is.na(variable)) %>%
  group_by(model, variable, scenario) %>%
  summarise(
    n_files               = n(),
    n_distinct_units      = n_distinct(na.omit(units)),
    units_found           = paste(unique(na.omit(units)), collapse = " | "),
    n_distinct_calendars  = n_distinct(na.omit(calendar)),
    calendars_found       = paste(unique(na.omit(calendar)), collapse = " | "),
    typical_step_days     = median(step_days_median, na.rm = TRUE),
    all_steps_daily       = all(round(step_days_median, 4) == 1, na.rm = TRUE),
    any_non_uniform_step  = any(!step_uniform, na.rm = TRUE),
    n_spillover_suspected = sum(spillover_suspected, na.rm = TRUE),
    earliest_date         = min(first_date, na.rm = TRUE),
    latest_date            = max(last_date, na.rm = TRUE),
    n_files_with_read_error = sum(!is.na(error)),
    .groups = "drop"
  ) %>%
  mutate(
    expected_units = unname(expected_units[variable]),
    units_match_expected = str_trim(tolower(units_found)) == str_trim(tolower(expected_units))
  )


readr::write_csv(variable_summary, paste0(diag_outpath, "/PolarRes26_diagnostics_variable_summary_", models_to_run, ".csv"))
message("Wrote variable-level summary: ",
        file.path(diag_outpath, models_to_run, "PolarRes26_diagnostics_variable_summary.csv"))

# ---- 5. Date-coverage check: gaps & duplicates ----------------------------------
# For each (model, variable, scenario), build the set of dates actually
# covered across ALL its files (each file's first->last date, inclusive --
# a reasonable approximation given a daily step), compare to a full daily
# sequence, and report gaps + duplicated dates (same kind of issue as the
# old missing rsus months).

check_coverage <- function(df) {
  ranges <- df %>% filter(!is.na(first_date) & !is.na(last_date))
  if (nrow(ranges) == 0) return(tibble())
  
  all_dates <- unlist(purrr::map2(ranges$first_date, ranges$last_date,
                                  ~ seq(.x, .y, by = "day")))
  all_dates <- as.Date(all_dates, origin = "1970-01-01")
  
  full_seq <- seq(min(all_dates), max(all_dates), by = "day")
  missing  <- as.Date(setdiff(full_seq, all_dates), origin = "1970-01-01")
  dupes    <- as.Date(unique(all_dates[duplicated(all_dates)]), origin = "1970-01-01")
  
  tibble(
    n_missing_dates          = length(missing),
    missing_dates_sample     = paste(head(missing, 10), collapse = ", "),
    n_duplicated_dates       = length(dupes),
    duplicated_dates_sample  = paste(head(dupes, 10), collapse = ", ")
  )
}

date_gaps <- file_level_df %>%
  filter(!is.na(variable)) %>%
  group_by(model, variable, scenario) %>%
  group_modify(~ check_coverage(.x)) %>%
  ungroup()

readr::write_csv(date_gaps, paste0(diag_outpath, "/PolarRes26_diagnostics_date_gaps_", models_to_run,".csv"))
message("Wrote date-gap report: ",
        file.path(diag_outpath, models_to_run, "PolarRes26_diagnostics_date_gaps.csv"))

# ---- 5b. Units for the same variable across ALL models -------------------------
# Useful if you'll eventually want to combine/compare across models -- flags
# if e.g. HCLIM's pr and RACMO's pr aren't stored in the same units.

units_by_variable <- file_level_df %>%
  filter(!is.na(variable)) %>%
  group_by(variable) %>%
  summarise(
    units_across_all_models = paste(unique(na.omit(units)), collapse = " | "),
    n_distinct_across_models = n_distinct(na.omit(units)),
    .groups = "drop"
  )

readr::write_csv(units_by_variable, paste0(diag_outpath, "/PolarRes26_diagnostics_units_by_variable_", models_to_run,".csv"))

# ---- 6. Console warnings worth a human's attention ------------------------------

message("\n================ SUMMARY OF THINGS TO CHECK ================\n")

bad_units <- variable_summary %>% filter(!units_match_expected | n_distinct_units > 1)
if (nrow(bad_units) > 0) {
  message("UNITS DIFFERENT FROM WHAT THE OLD SCRIPTS ASSUMED (or inconsistent ",
          "within a variable/scenario) -- check before assuming old conversions apply:")
  print(as.data.frame(bad_units %>%
                        select(model, variable, scenario, units_found, expected_units)))
}

bad_cal <- variable_summary %>%
  filter(n_distinct_calendars > 1 |
           !calendars_found %in% c("", "standard", "gregorian", "proleptic_gregorian"))
if (nrow(bad_cal) > 0) {
  message("\nNON-STANDARD OR INCONSISTENT CALENDARS -- this affects day-counting ",
          "(degree days, day-of-year math):")
  print(as.data.frame(bad_cal %>% select(model, variable, scenario, calendars_found)))
}

not_daily <- variable_summary %>% filter(!all_steps_daily | any_non_uniform_step)
if (nrow(not_daily) > 0) {
  message("\nNOT CONFIRMED DAILY / NON-UNIFORM TIME STEPS:")
  print(as.data.frame(not_daily %>%
                        select(model, variable, scenario, typical_step_days, any_non_uniform_step)))
}

spillover <- variable_summary %>% filter(n_spillover_suspected > 0)
if (nrow(spillover) > 0) {
  message("\nSPILLOVER SUSPECTED (a file's actual last date is after its filename's ",
          "nominal end date) -- the old 'drop the last layer' fix may still be needed:")
  print(as.data.frame(spillover %>% select(model, variable, scenario, n_spillover_suspected)))
}

gappy <- date_gaps %>% filter(n_missing_dates > 0)
if (nrow(gappy) > 0) {
  message("\nMISSING DATES WITHIN COVERAGE RANGE (potential corrupt/absent files, ",
          "same kind of issue as the old missing rsus months):")
  print(as.data.frame(gappy %>%
                        select(model, variable, scenario, n_missing_dates, missing_dates_sample)))
}

dupey <- date_gaps %>% filter(n_duplicated_dates > 0)
if (nrow(dupey) > 0) {
  message("\nDUPLICATED DATES ACROSS FILES (chunks may overlap):")
  print(as.data.frame(dupey %>%
                        select(model, variable, scenario, n_duplicated_dates, duplicated_dates_sample)))
}

cross_model <- units_by_variable %>% filter(n_distinct_across_models > 1)
if (nrow(cross_model) > 0) {
  message("\nUNITS DIFFER FOR THE SAME VARIABLE ACROSS MODELS:")
  print(as.data.frame(cross_model))
}

errored <- file_level_df %>% filter(!is.na(error))
if (nrow(errored) > 0) {
  message("\nFILES THAT COULD NOT BE READ AT ALL:")
  print(as.data.frame(errored %>% select(model, variable, path, error)))
}

message("\n===============================================================\n")
message("Diagnostics complete. Review the CSVs in: ", diag_outpath)
message("Once you've checked these outputs (especially units, calendars, ",
        "spillover, and missing dates), we can write Script 2 to summarise ",
        "into climatologies using years_hist/years_mid/years_future.")