lib_loc <- paste(getwd(),"/r_lib",sep="")
library(glue)

dirpath <- "/mnt/hpccs01/home/patterc2/n11222026/AntarcticFutureHabitat/Data/AntAirICE"

grid <- expand.grid(quarter = 1:4, year = 2003:2014)
files <- glue("AntAir_ICE_{grid$year}_{grid$quarter}.zip")

validate_extraction <- function(localfile, dirpath) {
  listing <- system2("unzip", args = c("-l", shQuote(localfile)), stdout = TRUE)
  data_lines <- grep("^\\s*[0-9]+\\s+[0-9-]{10}\\s+[0-9:]{5}\\s+", listing, value = TRUE)
  if (length(data_lines) == 0) return(FALSE)
  sizes <- as.numeric(sub("^\\s*([0-9]+).*$", "\\1", data_lines))
  names <- basename(sub("^\\s*[0-9]+\\s+[0-9-]{10}\\s+[0-9:]{5}\\s+", "", data_lines))
  ok <- logical(length(names))
  for (i in seq_along(names)) {
    fpath <- file.path(dirpath, names[i])
    ok[i] <- file.exists(fpath) && file.size(fpath) == sizes[i]
  }
  if (!all(ok)) message(glue("  Bad/missing: {paste(names[!ok], collapse=', ')}"))
  all(ok)
}

for (file_name in files) {
  localfile <- file.path(dirpath, file_name)
  if (!file.exists(localfile)) {
    message(glue("MISSING ZIP: {file_name} — skipping (needs re-download)"))
    next
  }
  
  if (validate_extraction(localfile, dirpath)) {
    message(glue("OK: {file_name}"))
    next
  }
  
  message(glue("REPAIRING: {file_name}"))
  exit_code <- system2("unzip", args = c("-o", "-j", shQuote(localfile), "-d", shQuote(dirpath)))
  
  if (validate_extraction(localfile, dirpath)) {
    message(glue("  fixed: {file_name}"))
  } else {
    message(glue("  STILL BAD after re-extract: {file_name} — check zip integrity (unzip -t)"))
  }
}