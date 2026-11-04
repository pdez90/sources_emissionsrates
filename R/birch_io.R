# CDPHE MT GOAT Adams County Birch Street (AQDx v3) parsing.
# Zips hold weekly 1-s CSVs (CR line endings) per parameter; timestamps are fixed MST (-07:00).
# Qualifiers: MD = below MDL (value kept, incl. negatives); AQ = zero/bump QC; BD/BR/AN = invalid; validity_code 9 = invalid.
# Wind speed unit code 012 = mph.  H2S clock fix: stamps before BIRCH_CLOCK_CUT are 6 h early (6.0-h gap at the cut;
# before it, H2S lags Birch's own wind by ~6 h and its diurnal cycle is shifted).
BIRCH_PARAMS <- c("hydrogen_sulfide", "wind_speed", "wind_direction")
birch_quarters <- function() sort(sub(".*_(20\\d\\dQ\\d)\\.zip$", "\\1", list.files(RAW_BIRCH, "^MTGOAT_AdamsCountyBirchStreetBeta_20\\d\\dQ\\d\\.zip$")))
birch_parse <- function(q, param) {
  out <- file.path(INTERIM, "birch", paste0(param, "_", q, ".rds")); if (file.exists(out)) return(readRDS(out))
  z <- file.path(RAW_BIRCH, paste0("MTGOAT_AdamsCountyBirchStreetBeta_", q, ".zip"))
  mem <- sort(grep(paste0("^weekly_files_", param, "/.*\\.csv$"), unzip(z, list = TRUE)$Name, value = TRUE))
  d <- rbindlist(lapply(mem, function(m) fread(cmd = sprintf("unzip -p '%s' '%s' | tr '\\r' '\\n'", z, m),
                 select = c("datetime", "parameter_value", "validity_code"), colClasses = list(character = "datetime"))))
  stopifnot(all(substr(d$datetime, 20, 25) == "-07:00"))
  r <- data.table(t = as.numeric(as.POSIXct(substr(d$datetime, 1, 19), format = "%Y-%m-%dT%H:%M:%S", tz = TZ_MST)),
                  v = suppressWarnings(as.numeric(d$parameter_value)), valid = d$validity_code)
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE); saveRDS(r, out); r
}
birch_series <- function(q, param, fix_clock = TRUE) {
  d <- copy(birch_parse(q, param)); d[valid != 1, v := NA]
  if (fix_clock && param == "hydrogen_sulfide") d[t < as.numeric(BIRCH_CLOCK_CUT), t := t + 6*3600]
  d[, .(v = mean(v, na.rm = TRUE)), by = t][order(t)][is.nan(v), v := NA]
}
birch_one_minute <- function(q, fix_clock = TRUE) {
  h <- birch_series(q, "hydrogen_sulfide", fix_clock); ws <- birch_series(q, "wind_speed"); wd <- birch_series(q, "wind_direction")
  w <- merge(ws, wd, by = "t", suffixes = c("_ws", "_wd")); w[, `:=`(ms = v_ws*0.44704)]
  w[, `:=`(u = -ms*sin(v_wd*pi/180), v = -ms*cos(v_wd*pi/180))]
  h[, b := ceiling(t/60)*60]; w[, b := ceiling(t/60)*60]
  H <- h[, .(h2s = mean(v, na.rm = TRUE), n_h2s = sum(!is.na(v))), by = b]
  W <- w[, .(u = mean(u, na.rm = TRUE), v = mean(v, na.rm = TRUE), ws = mean(ms, na.rm = TRUE), n_wind = sum(!is.na(u))), by = b]
  D <- merge(H, W, by = "b", all = TRUE)
  D[is.na(n_h2s) | n_h2s < 30, h2s := NA]; D[is.na(n_wind) | n_wind < 30, `:=`(u = NA, v = NA, ws = NA)]
  D[is.na(n_h2s), n_h2s := 0L]; D[is.na(n_wind), n_wind := 0L]
  for (k in c("h2s", "u", "v", "ws")) set(D, which(is.nan(D[[k]])), k, NA)
  D[, quarter := q]; D
}
birch_build <- function() {
  D <- rbindlist(lapply(birch_quarters(), birch_one_minute)); D <- D[!duplicated(b)][order(b)]
  out <- D[, .(t_mst = format(as.POSIXct(b, origin = "1970-01-01", tz = TZ_MST), "%Y-%m-%d %H:%M:%S"),
               h2s = round(h2s, 4), n_h2s, u = round(u, 4), v = round(v, 4), ws = round(ws, 4), n_wind, quarter)]
  fwrite(out, file.path(PROCESSED, "birch_1min.csv.gz")); out
}
birch_load_1min <- function() {
  d <- fread(file.path(PROCESSED, "birch_1min.csv.gz"), colClasses = list(character = "t_mst"))
  d[, t := as.POSIXct(t_mst, tz = TZ_MST)][]
}
# 5-min means (label = interval end), requiring >= 3 valid minutes
birch_five_min <- function(b1 = birch_load_1min()) {
  b1[, b := bin_end(t, 300)]
  d <- b1[, .(h2s = mean(h2s, na.rm = TRUE), nh = sum(!is.na(h2s)), u = mean(u, na.rm = TRUE), v = mean(v, na.rm = TRUE), nw = sum(!is.na(u))), by = b]
  d[nh < 3, h2s := NA][nw < 3, `:=`(u = NA, v = NA)]
  for (k in c("h2s", "u", "v")) set(d, which(is.nan(d[[k]])), k, NA)
  setnames(d, "b", "t"); d[order(t), .(t, h2s, u, v)]
}
