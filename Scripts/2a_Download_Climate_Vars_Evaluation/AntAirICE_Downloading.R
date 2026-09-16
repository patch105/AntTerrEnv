
# HPC version
lib_loc <- paste(getwd(),"/r_lib",sep="")

library(terra)
library(glue)
library(gdalraster, lib.loc = lib_loc)
library(vapour, lib.loc = lib_loc)

# extract the arguments provided in the command line
args <- commandArgs(trailingOnly = TRUE)
# The first argument is now the job index
job_index <- as.integer(args[1])

# Stagger jobs to avoid hammering the server simultaneously
Sys.sleep(runif(1, min = job_index * 2, max = job_index * 2 + 30))

# URL for file location
u <- "https://download.pangaea.de/dataset/954750/files"

# Build the full file list from the known naming pattern
# Files run from 2003 to 2014 (adjust end year as needed), quarters 1-4
years <- 2003:2014
quarters <- 1:4
# files <- as.vector(outer(years, quarters, function(y, q) glue("AntAir_ICE_{y}_{q}.zip")))
# Correct ordering: 2003_1, 2003_2, 2003_3, 2003_4, 2004_1, ...
grid <- expand.grid(quarter = 1:4, year = 2003:2014)  # expand.grid fills first argument fastest
files <- glue("AntAir_ICE_{grid$year}_{grid$quarter}.zip")

# Use job_index to select the file for this job (1-based indexing for SLURM arrays)
i <- job_index
file_name <- files[i]

dirpath <- "/mnt/hpccs01/home/n11222026/AntarcticFutureHabitat/Data/AntAirICE"
file_url <- glue("{u}/{file_name}")
localfile <- glue("{dirpath}/{file_name}")

# --- Check 1: zip already exists ---
if (file.exists(localfile)) {
  message(glue("Job {job_index}: Zip exists already — {file_name}"))
} else {
  message(glue("Job {job_index}: Downloading {file_name}"))
  
  max_attempts <- 5
  for (attempt in 1:max_attempts) {
    result <- tryCatch({
      curl::curl_download(file_url, localfile)
      "success"
    }, error = function(e) {
      message(glue("Attempt {attempt} failed: {e$message}"))
      "error"
    })
    
    if (result == "success") {
      message(glue("Download complete: {localfile}"))
      break
    } else if (attempt < max_attempts) {
      wait <- 30 * attempt
      message(glue("Retrying in {wait}s..."))
      Sys.sleep(wait)
    } else {
      stop(glue("Failed to download {file_name} after {max_attempts} attempts"))
    }
  }
}

# --- Check 2: already unzipped ---
# Since files are extracted flat (junkpaths), check for any extracted file
# whose name contains the year_quarter stem of this zip (e.g. "2003_1")
zip_stem <- gsub("AntAir_ICE_|\\.zip", "", file_name)  # e.g. "2003_1"
extracted_files <- list.files(dirpath, pattern = zip_stem, full.names = FALSE)
extracted_files <- extracted_files[!grepl("\\.zip$", extracted_files)]  # exclude the zip itself

if (length(extracted_files) > 0) {
  message(glue("Job {job_index}: Unzipped already — found {length(extracted_files)} file(s) matching '{zip_stem}'"))
} else {
  message(glue("Unzipping {file_name} into {dirpath}"))
  exit_code <- system2("unzip", args = c("-o", "-j", shQuote(localfile), "-d", shQuote(dirpath)))
  
  if (exit_code == 0) {
    message(glue("Unzip complete for {file_name}"))
  } else {
    stop(glue("Unzip failed for {file_name} with exit code {exit_code}"))
  }
}
