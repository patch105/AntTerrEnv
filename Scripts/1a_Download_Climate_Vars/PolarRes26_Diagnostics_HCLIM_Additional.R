# ==============================================================================
# PolarRes26 -- DIAGNOSTICS (SUBSET): HCLIM siconca / rsus / snc / snm
# ==============================================================================
# A separate, simplified diagnostics pass for a small subset of variables
# (siconca, rsus, snc, snm) that came from a different source/route than the
# main PolarRes26 download: flat monthly-chunked files, not organised into
# per-variable folders and not one folder per model/GCM, originally from
# http://prudence.dmi.dk/data/temp/JAT/Charlotte/ (its "midfut" and "CESM"
# subfolders, plus loose files directly in the base folder). Filenames look
# like:
#   rsus_fp_ANT11_ANT11_ALADIN43_v1_CESM2_r11i1p1f1_ssp370_day_204001010000-204002010000.nc
#   rsus_fp_ANT11_ANT11_ALADIN43_v1_MPI_r1i1p1f1_historical_day_198601010000-198602010000.nc
#   snm_sfx_ANT11_ANT11_ALADIN43_v1_CESM2_r11i1p1f1_ssp370_day_209611010000-209612010000.nc
#
# This is the same general shape of data -- and potentially the same known
# idiosyncrasy -- described earlier for the old PolarRes .tif pipeline: one
# file per calendar month, where the filename's end-date token is the FIRST
# DAY OF THE NEXT MONTH (00:00) rather than the last day of the file's own
# month. That's exactly the "monthly raster includes the name of the first
# day of the next month" problem that pipeline had to work around by
# dropping each file's last layer. This script checks, per file, whether
# that's still happening here, by comparing each file's actual layer count
# to the number of calendar days in its own month (taken from the START
# token, which is unambiguous):
#   nlyr == days_in_month       -> clean, no extra layer
#   nlyr == days_in_month + 1   -> classic extra-layer idiosyncrasy
#   anything else               -> something else is wrong, worth a look
#
# GCM, scenario, and variable aren't in folder names here -- they're parsed
# straight out of each filename instead. Everything not matching one of the
# 4 variables of interest is ignored. subset_base_dir is scanned
# recursively, so this naturally picks up loose files directly in it plus
# everything under any subfolder (midfut, CESM, or otherwise) in one pass --
# no need to list subfolders separately.
#
# Run one GCM at a time (job-array style, matching the main diagnostics
# script) -- args[1] picks a row from `subset_gcms`, e.g. 1 = CESM2,
# 2 = MPI. With no argument, all GCMs found are done in one run instead.
#
# Outputs (written to Data/PolarRes26/Diagnostics/, alongside the main
# diagnostics outputs, with "_subset_<GCM>" in every name -- GCM uses the
# same naming convention as everywhere else in this project, i.e. "MPI"
# in the filenames themselves becomes "MPI_ESM1" in these output names):
#   - PolarRes26_diagnostics_subset_file_level_<GCM>.csv
#   - PolarRes26_diagnostics_subset_variable_summary_<GCM>.csv
#   - PolarRes26_diagnostics_subset_month_gaps_<GCM>.csv
# Plus a console summary per GCM.
# ==============================================================================

# ---- 0. Setup ------------------------------------------------------------------

lib_loc <- paste(getwd(), "/r_lib_new", sep = "")
# .libPaths(lib_loc)  # uncomment if r_lib isn't already on the search path

library(terra)
library(ncdf4, lib.loc = lib_loc)
library(here)
library(lubridate)
library(stringr)
library(dplyr, lib.loc = lib_loc)
library(purrr, lib.loc = lib_loc)
library(readr)
library(tibble)

has_ncdf4 <- requireNamespace("ncdf4", quietly = TRUE)

# ---- 1. Configuration -----------------------------------------------------------

# EDIT ME: wherever the prudence.dmi.dk/.../Charlotte files have been
# mirrored to locally.
subset_base_dir <- "/mnt/hpccs01/home/patterc2/n11222026/AntarcticFutureHabitat/Data/PolarRes26"

subset_variables <- c("siconca", "rsus", "snc", "snm")

expected_units_subset <- c(
  siconca = "%",
  rsus    = "W m-2",
  snc     = "%",
  snm     = "kg m-2 s-1"
)

# GCM tokens as they appear in the filenames themselves, mapped to the
# naming convention used everywhere else in this project (HCLIM_CESM2 /
# HCLIM_MPI_ESM1 etc.) for labelling the output files.
subset_gcms <- c("CESM2", "MPI")
gcm_output_label <- c(CESM2 = "CESM2", MPI = "MPI_ESM1")

# Same output location as the main diagnostics script.
diag_outpath <- here("Data/PolarRes26/Diagnostics")
dir.create(diag_outpath, recursive = TRUE, showWarnings = FALSE)

# Run one GCM at a time (job-array style, matching the main diagnostics
# script), or all GCMs found in one go if no argument is given.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && !is.na(suppressWarnings(as.integer(args[1])))) {
  gcms_to_run <- subset_gcms[as.integer(args[1])]
} else {
  gcms_to_run <- subset_gcms
}
message("Running subset diagnostics for GCM(s): ", paste(gcms_to_run, collapse = ", "))

# ---- 2. Helper functions ---------------------------------------------------------

# 2a. Pull variable / GCM / scenario / nominal month straight out of the
#     filename -- there's no folder structure to lean on here.
parse_filename <- function(filenames) {
  variable <- map_chr(filenames, function(fn) {
    hit <- subset_variables[str_starts(fn, paste0(subset_variables, "_"))]
    if (length(hit) == 1) hit else NA_character_
  })
  gcm      <- str_extract(filenames, "CESM2|MPI")
  scenario <- str_extract(filenames, "ssp[0-9]+|historical")
  
  # Trailing "..._<8-digit-date><4-digit-time>-<8-digit-date><4-digit-time>.nc".
  # HHMM is always 0000 in every example seen so far -- only the date part
  # is used below. If that assumption ever breaks (a non-zero HHMM turns
  # up), this regex will just fail to match that file, which will show up
  # as an NA nominal_start/nominal_end_token rather than silently misparsing it.
  m <- str_match(filenames, "_(\\d{8})\\d{4}-(\\d{8})\\d{4}\\.nc$")
  nominal_start      <- ymd(m[, 2])
  nominal_end_token  <- ymd(m[, 3])  # literal date in the token (first of NEXT month)
  
  tibble(variable, gcm, scenario, nominal_start, nominal_end_token)
}

# 2b. Read the metadata we need from one file, without loading pixel data
#     (same approach as the main diagnostics script).
read_file_meta <- function(path) {
  out <- list(
    nlyr = NA_integer_, first_time = as.POSIXct(NA), last_time = as.POSIXct(NA),
    step_days_median = NA_real_, step_uniform = NA, calendar = NA_character_,
    units = NA_character_, main_var = NA_character_, error = NA_character_
  )
  
  main_var <- NA_character_
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

# 2c. Missing AND duplicate months within one (variable, scenario)'s
#     observed span. Duplicates matter here because a GCM's files can be
#     split across more than one folder (e.g. a CESM subfolder as well as
#     midfut) -- if the same month exists in two different files, this
#     catches it, since a missing-months check alone wouldn't.
check_month_coverage <- function(df) {
  months_all <- floor_date(df$nominal_start, "month")
  months_present <- unique(months_all)
  
  full_seq <- seq(min(months_present), max(months_present), by = "month")
  missing  <- setdiff(full_seq, months_present)
  
  dupe_months <- unique(months_all[duplicated(months_all)])
  
  tibble(
    n_missing_months = length(missing),
    missing_months_sample = paste(head(as.Date(missing, origin = "1970-01-01"), 12), collapse = ", "),
    n_duplicate_months = length(dupe_months),
    duplicate_months_sample = paste(head(as.Date(dupe_months, origin = "1970-01-01"), 12), collapse = ", ")
  )
}

# ---- 3. Discover every matching file (once -- filtered per GCM below) ----------

if (!dir.exists(subset_base_dir)) {
  stop("Folder not found: ", subset_base_dir)
}

all_files <- list.files(subset_base_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
if (length(all_files) == 0) {
  stop("No .nc files found under: ", subset_base_dir)
}

parsed <- parse_filename(basename(all_files))
found  <- bind_cols(tibble(path = all_files, filename = basename(all_files)), parsed)

unmatched <- found %>% filter(is.na(variable))
if (nrow(unmatched) > 0) {
  message(nrow(unmatched), " file(s) didn't match one of siconca/rsus/snc/snm -- ignored. Examples:")
  print(head(unmatched$filename, 10))
}
found <- found %>% filter(!is.na(variable))
if (nrow(found) == 0) stop("None of the files found matched siconca/rsus/snc/snm.")

# ---- 4. Run diagnostics for each requested GCM, writing its own set of CSVs ----

for (gcm in gcms_to_run) {
  
  gcm_label <- unname(gcm_output_label[gcm])
  message("\n=== ", gcm_label, " ===")
  
  found_gcm <- found %>% filter(gcm == !!gcm)
  if (nrow(found_gcm) == 0) {
    message("No matching files found for GCM: ", gcm)
    next
  }
  
  message(nrow(found_gcm), " matching file(s), ", length(unique(found_gcm$variable)), " variable(s).")
  
  meta_list <- vector("list", nrow(found_gcm))
  for (i in seq_len(nrow(found_gcm))) {
    if (i %% 50 == 0) message("  ...", i, "/", nrow(found_gcm))
    meta_list[[i]] <- read_file_meta(found_gcm$path[i])
  }
  meta_df <- bind_rows(lapply(meta_list, as_tibble))
  
  file_level_df <- bind_cols(found_gcm, meta_df) %>%
    mutate(
      first_date = as.Date(first_time),
      last_date  = as.Date(last_time),
      expected_days_in_month  = lubridate::days_in_month(nominal_start),
      extra_layer_suspected   = !is.na(nlyr) & nlyr == expected_days_in_month + 1,
      layer_count_unexplained = !is.na(nlyr) &
        !(nlyr %in% c(expected_days_in_month, expected_days_in_month + 1))
    )
  
  readr::write_csv(file_level_df,
                   file.path(diag_outpath, sprintf("PolarRes26_diagnostics_subset_file_level_%s.csv", gcm_label)))
  message("  wrote PolarRes26_diagnostics_subset_file_level_", gcm_label, ".csv")
  
  # ---- roll-up per (variable, scenario) ----
  variable_summary <- file_level_df %>%
    group_by(variable, scenario) %>%
    summarise(
      n_files                   = n(),
      n_distinct_units          = n_distinct(na.omit(units)),
      units_found               = paste(unique(na.omit(units)), collapse = " | "),
      n_distinct_calendars      = n_distinct(na.omit(calendar)),
      calendars_found           = paste(unique(na.omit(calendar)), collapse = " | "),
      n_extra_layer_suspected   = sum(extra_layer_suspected, na.rm = TRUE),
      n_layer_count_unexplained = sum(layer_count_unexplained, na.rm = TRUE),
      any_non_uniform_step      = any(!step_uniform, na.rm = TRUE),
      earliest_month            = min(nominal_start, na.rm = TRUE),
      latest_month              = max(nominal_start, na.rm = TRUE),
      n_files_with_read_error   = sum(!is.na(error)),
      .groups = "drop"
    ) %>%
    mutate(
      expected_units = unname(expected_units_subset[variable]),
      units_match_expected = str_trim(tolower(units_found)) == str_trim(tolower(expected_units))
    )
  
  readr::write_csv(variable_summary,
                   file.path(diag_outpath, sprintf("PolarRes26_diagnostics_subset_variable_summary_%s.csv", gcm_label)))
  message("  wrote PolarRes26_diagnostics_subset_variable_summary_", gcm_label, ".csv")
  
  # ---- missing / duplicate months per (variable, scenario) ----
  month_gaps <- file_level_df %>%
    group_by(variable, scenario) %>%
    group_modify(~ check_month_coverage(.x)) %>%
    ungroup()
  
  readr::write_csv(month_gaps,
                   file.path(diag_outpath, sprintf("PolarRes26_diagnostics_subset_month_gaps_%s.csv", gcm_label)))
  message("  wrote PolarRes26_diagnostics_subset_month_gaps_", gcm_label, ".csv")
  
  # ---- console summary for this GCM ----
  message("\n---- ", gcm_label, ": things to check ----")
  
  bad_units <- variable_summary %>% filter(!units_match_expected | n_distinct_units > 1)
  if (nrow(bad_units) > 0) {
    message("UNITS DIFFERENT FROM EXPECTED (or inconsistent):")
    print(as.data.frame(bad_units %>% select(variable, scenario, units_found, expected_units)))
  }
  
  bad_cal <- variable_summary %>%
    filter(n_distinct_calendars > 1 |
             !calendars_found %in% c("", "standard", "gregorian", "proleptic_gregorian"))
  if (nrow(bad_cal) > 0) {
    message("\nNON-STANDARD OR INCONSISTENT CALENDARS:")
    print(as.data.frame(bad_cal %>% select(variable, scenario, calendars_found)))
  }
  
  extra_layer <- variable_summary %>% filter(n_extra_layer_suspected > 0)
  if (nrow(extra_layer) > 0) {
    message("\nCLASSIC 'EXTRA LAYER' IDIOSYNCRASY STILL PRESENT (nlyr = days_in_month + 1) -- ",
            "these files likely still need their last layer dropped before use:")
    print(as.data.frame(extra_layer %>% select(variable, scenario, n_files, n_extra_layer_suspected)))
  }
  
  unexplained <- variable_summary %>% filter(n_layer_count_unexplained > 0)
  if (nrow(unexplained) > 0) {
    message("\nUNEXPLAINED LAYER COUNTS (neither days_in_month nor days_in_month+1) -- ",
            "check these files individually, this isn't the known idiosyncrasy:")
    print(as.data.frame(unexplained %>% select(variable, scenario, n_files, n_layer_count_unexplained)))
  }
  
  not_uniform <- variable_summary %>% filter(any_non_uniform_step)
  if (nrow(not_uniform) > 0) {
    message("\nNON-UNIFORM TIME STEPS WITHIN AT LEAST ONE FILE:")
    print(as.data.frame(not_uniform %>% select(variable, scenario)))
  }
  
  gappy <- month_gaps %>% filter(n_missing_months > 0)
  if (nrow(gappy) > 0) {
    message("\nMISSING MONTHS (the old 'missing raster, fill with blanks' issue) -- ",
            "gaps within the observed date range:")
    print(as.data.frame(gappy %>% select(variable, scenario, n_missing_months, missing_months_sample)))
  }
  
  dupe_months <- month_gaps %>% filter(n_duplicate_months > 0)
  if (nrow(dupe_months) > 0) {
    message("\nDUPLICATE MONTHS (the same calendar month appears in more than one file) -- ",
            "worth checking for stale/overlapping files, since this GCM's data may be ",
            "split across more than one folder:")
    print(as.data.frame(dupe_months %>% select(variable, scenario, n_duplicate_months, duplicate_months_sample)))
  }
  
  errored <- file_level_df %>% filter(!is.na(error))
  if (nrow(errored) > 0) {
    message("\nFILES THAT COULD NOT BE READ AT ALL:")
    print(as.data.frame(errored %>% select(variable, path, error)))
  }
}

message("\n=====================================================================\n")
message("Subset diagnostics complete. Review the CSVs (with '_subset_<GCM>' in the name) in: ", diag_outpath)