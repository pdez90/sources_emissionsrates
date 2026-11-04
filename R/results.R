# key_numbers.json: headline numbers written by each step and checked by scripts/99_vet.R
KEYFILE <- file.path(RESULTS, "key_numbers.json")
record <- function(section, ...) {
  d <- if (file.exists(KEYFILE)) jsonlite::read_json(KEYFILE, simplifyVector = TRUE) else list()
  kv <- list(...); if (is.null(d[[section]])) d[[section]] <- list()
  for (k in names(kv)) d[[section]][[k]] <- kv[[k]]
  jsonlite::write_json(d, KEYFILE, auto_unbox = TRUE, digits = 6, pretty = TRUE)
}
load_keys <- function() if (file.exists(KEYFILE)) jsonlite::read_json(KEYFILE, simplifyVector = TRUE) else list()
# write a table with POSIXct columns as MST "YYYY-MM-DD HH:MM:SS" text
wcsv <- function(d, file) {
  d <- as.data.table(copy(d))
  for (k in names(d)) if (inherits(d[[k]], "POSIXct")) set(d, j = k, value = format(d[[k]], "%Y-%m-%d %H:%M:%S", tz = TZ_MST))
  fwrite(d, file)
}
