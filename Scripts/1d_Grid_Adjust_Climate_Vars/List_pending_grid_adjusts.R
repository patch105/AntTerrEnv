# ==============================================================================
# PolarRes26 -- list_pending_jobs.R
# ------------------------------------------------------------------------------
# Rebuilds the exact same job table and expected-output paths as Script 3
# (sections 1-4), WITHOUT reading any raster data. Its only purpose is to
# tell you which job_index values in Script 3 would NOT be skipped, so you
# can submit only those to your scheduler instead of all ~3000+.
#
# This duplicates Script 3's config and path-construction logic on purpose
# (rather than sourcing Script 3 itself), because Script 3 does raster work
# and quit()s as soon as job_index is missing/invalid -- it's not written to
# be sourced. If you ever change the output-path logic, filenames, or the
# is_sea_ice/is_future detection in Script 3, mirror the change here too.
#
# Usage:
#   Rscript list_pending_jobs.R
#     -> prints the full job table with a `pending` column, a summary count,
#        and (on the last line) a comma-separated list of pending job_index
#        values ready to paste into e.g. `sbatch --array=<that list>`.
#
#   Rscript list_pending_jobs.R --indices-only
#     -> prints ONLY the comma-separated list of pending indices (nothing
#        else), so it can be captured directly in a shell variable, e.g.:
#          PENDING=$(Rscript list_pending_jobs.R --indices-only)
#          sbatch --array=$PENDING run_script3.sh
# ==============================================================================

suppressMessages({
  lib_loc <- paste(getwd(), "/r_lib_new", sep = "")
  library(dplyr, lib.loc = lib_loc)
  library(purrr, lib.loc = lib_loc)
  library(here)
})

indices_only <- "--indices-only" %in% commandArgs(trailingOnly = TRUE)

# ---- Config (must match Script 3 sections 1/2/4 exactly) -----------------------

models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5", "RACMO_CESM2",
            "RACMO_MPI_ESM1", "RACMO_ERA5", "MetUM_ERA5")

input_base  <- here("Data/Environmental_predictors/PolarRes26")
output_base <- here("Data/Environmental_predictors/PolarRes26/Regridded")

sea_ice_buffer_km <- c(2, 5, 10, 50, 100, 200)

# ---- Build job table (identical to Script 3 section 2) -------------------------

job_table <- map_dfr(models, function(model) {
  model_dir <- file.path(input_base, model)
  if (!dir.exists(model_dir)) {
    return(tibble())
  }
  files <- list.files(model_dir, pattern = "\\.tif$", full.names = TRUE)
  tibble(model = model, input_path = files, filename = basename(files))
})

if (nrow(job_table) == 0) {
  stop("No .tif files found under ", input_base, " -- has Script 2 been run yet?")
}

# ---- For every job, compute expected outputs and check existence (Script 3 section 4) --

job_table <- job_table %>%
  mutate(job_index = row_number()) %>%
  rowwise() %>%
  mutate(
    out_dir       = file.path(output_base, model),
    stem          = tools::file_path_sans_ext(filename),
    is_sea_ice    = grepl("Sea_Ice_Concentration", filename, fixed = TRUE),
    is_future     = grepl("_FUTURE_", filename, fixed = TRUE),
    expected_outputs = list({
      coastline_path     <- file.path(out_dir, paste0(stem, "_COASTLINE.tif"))
      icefree_path       <- file.path(out_dir, paste0(stem, "_ICEFREE.tif"))
      icefree_futurepath <- file.path(out_dir, paste0(stem, "_ICEFREE_FUTURE.tif"))
      concentration_path <- file.path(out_dir, paste0(stem, "_CONCENTRATION.tif"))
      buffer_paths       <- file.path(out_dir, paste0(stem, "_BUFFER_", sea_ice_buffer_km, "km.tif"))
      
      if (is_sea_ice) {
        c(concentration_path, buffer_paths)
      } else {
        domain_mask_outputs <- c(coastline_path, icefree_path)
        if (is_future) domain_mask_outputs <- c(domain_mask_outputs, icefree_futurepath)
        domain_mask_outputs
      }
    }),
    all_outputs_exist = all(file.exists(unlist(expected_outputs))),
    pending = !all_outputs_exist
  ) %>%
  ungroup()

pending_indices <- job_table$job_index[job_table$pending]

if (indices_only) {
  cat(paste(pending_indices, collapse = ","), "\n")
  quit(save = "no", status = 0)
}

print(job_table %>% select(job_index, model, filename, is_sea_ice, is_future, pending), n = Inf)

message("\n", nrow(job_table), " total job(s); ",
        length(pending_indices), " pending (not yet fully processed); ",
        nrow(job_table) - length(pending_indices), " already complete (would be skipped).")

message("\nPending job_index values:")
message(paste(pending_indices, collapse = ","))

# ---- Manifest for downstream scripts --------------------------------------------
# Records exactly which Script-3 input files were pending at the moment this
# was run, so a later script can process "only what was just grid-adjusted"
# instead of re-scanning all of output_base (which would also pick up
# unrelated pre-existing outputs from earlier runs).
manifest <- job_table %>%
  filter(pending) %>%
  mutate(stem = tools::file_path_sans_ext(filename)) %>%
  select(job_index, model, filename, stem, is_sea_ice, is_future)

manifest_path <- here("Data/Environmental_predictors/PolarRes26/pending_manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)
message("\nWrote manifest of ", nrow(manifest), " pending job(s) to: ", manifest_path)