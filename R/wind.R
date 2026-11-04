# Wind: NOAA ISD (DEN, 725650-03017) parsing, CCND vane QA, and the regional transport wind used for inversion.
# CCND anemometers read ~0.4x Birch speed and several vanes drift/jump by tens of degrees between quarters
# (results/wind_qa_vs_den.csv). CM8 is stable vs DEN (-3 to -17 deg) and vs Birch (+6 deg), so the regional wind is
# Birch where available, otherwise CM8 with its speed scaled to Birch (CM8_SPEED_SCALE).
den_build <- function() {
  L <- unlist(lapply(sort(list.files(RAW_DEN, "^725650-03017-", full.names = TRUE)), readLines))
  L <- L[substr(L, 42, 46) == "FM-15"]
  d <- data.table(t_utc = as.POSIXct(paste(substr(L, 16, 23), substr(L, 24, 27)), format = "%Y%m%d %H%M", tz = "UTC"),
                  wd = as.integer(substr(L, 61, 63)), wsr = as.integer(substr(L, 66, 69)))
  d <- d[wd != 999 & wsr != 9999][!duplicated(t_utc)][, ws := wsr/10][, wsr := NULL]
  d[, t_mst := as.POSIXct(ceiling(as.numeric(t_utc)/3600)*3600, origin = "1970-01-01", tz = TZ_MST)]   # :53 obs -> hour ending
  d <- d[!duplicated(t_mst)]
  fwrite(d[, .(t_mst = format(t_mst, "%Y-%m-%d %H:%M:%S"), wd, ws)], file.path(PROCESSED, "den_hourly.csv.gz")); d
}
den_load <- function() {
  d <- fread(file.path(PROCESSED, "den_hourly.csv.gz"), colClasses = list(character = "t_mst"))
  d[, t := as.POSIXct(t_mst, tz = TZ_MST)][, `:=`(u = -ws*sin(wd*pi/180), v = -ws*cos(wd*pi/180))][]
}
# trailing-1h (12 x 5 min) vector-mean regional wind on grid `tt`: Birch where available (if given) else CM8 x scale
regional_wind <- function(P, tt, birch5 = NULL, smooth = 12) {
  i <- match(as.numeric(tt), as.numeric(P$t))
  U <- roll_mean(P$u[i, "CM8"]*CM8_SPEED_SCALE, smooth, smooth %/% 2); V <- roll_mean(P$v[i, "CM8"]*CM8_SPEED_SCALE, smooth, smooth %/% 2)
  if (!is.null(birch5)) {
    j <- match(as.numeric(tt), as.numeric(birch5$t))
    bu <- roll_mean(birch5$u[j], smooth, smooth %/% 2); bv <- roll_mean(birch5$v[j], smooth, smooth %/% 2)
    U <- ifelse(is.na(bu), U, bu); V <- ifelse(is.na(bv), V, bv)
  }
  list(U = U, V = V)
}
# per quarter x site: median direction offset (deg) and speed ratio against a reference hourly wind (u, v m/s)
wind_qa <- function(P, ref, min_ref = 3, min_site = 0.1) {
  hb <- bin_end(P$t, 3600); rows <- list()
  for (qq in setdiff(unique(P$quarter), "")) {
    m <- P$quarter == qq; r <- list(quarter = qq)
    for (s in SITE_LIST) {
      d <- data.table(b = hb[m], u = P$u[m, s], v = P$v[m, s])[, .(u = mean(u, na.rm = TRUE), v = mean(v, na.rm = TRUE)), by = b]
      j <- merge(d, ref[, .(b = t, ur = u, vr = v)], by = "b")[is.finite(u) & is.finite(ur)]
      ok <- sqrt(j$ur^2 + j$vr^2) > min_ref & sqrt(j$u^2 + j$v^2) > min_site
      if (sum(ok) < 100) { r[[paste0(s, "_dir")]] <- NA; r[[paste0(s, "_spd")]] <- NA; next }
      r[[paste0(s, "_dir")]] <- median(angdiff(wind_from(j$u, j$v), wind_from(j$ur, j$vr))[ok])
      r[[paste0(s, "_spd")]] <- median((sqrt(j$u^2 + j$v^2)/sqrt(j$ur^2 + j$vr^2))[ok])
    }
    rows[[qq]] <- as.data.table(r)
  }
  rbindlist(rows, fill = TRUE)
}
