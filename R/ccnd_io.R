# CCND quarterly workbook extraction and harmonization to a 5-min MST panel.
#
# Published CCND gas values are trailing 1-h means of 5-min data, stamped every 5 min in local time (MST/MDT).
# Wind is 1-min (mph). Known file issues fixed here:
#   * Q4 2023 wind: 'Date' and 'Time' columns disagree -> use 'Date' (full timestamp)
#   * Q1 2024 wind: stamped in UTC -> convert to America/Denver
#   * Q1 2025 wind: AM/PM lost (hours 12-23 stored as 00-11) -> add 12 h after each negative jump within a day
SHEETS <- c(h2s = "H2S-1hr", ws = "WS", wd = "WD", no2 = "NO2", so2 = "SO2", co = "CO", tvoc = "tVOC",
            benz = "Benzene-1hr", tol = "Toluene-1hr", ebz = "Ethylbenzene-1hr", xyl = "Xylenes-1hr")
GASES <- c("no2", "so2", "co", "tvoc", "benz", "tol", "ebz", "xyl")
tag_of <- function(f) gsub(" ", "_", sub("\\.xlsx$", "", basename(f)))
quarter_of <- function(tag) { m <- regmatches(tag, regexec("Q(\\d)_(\\d{4})", tag))[[1]]; c(q = as.integer(m[2]), y = as.integer(m[3])) }

SHEET_RX <- c(h2s = "^h2s-?(1-?h(ou)?r)?$", ws = "^ws$", wd = "^wd$", no2 = "^no2", so2 = "^so2", co = "^co$",
              tvoc = "^t?vocs?$", benz = "^benzene-?1-?h(ou)?r", tol = "^toluene-?1-?h(ou)?r",
              ebz = "^ethylbenzene-?1-?h(ou)?r", xyl = "^xylenes?-?1-?h(ou)?r")
find_sheet <- function(names, key) {
  nn <- tolower(gsub(" ", "", names)); hit <- names[grepl(SHEET_RX[[key]], nn)]
  if (length(hit)) hit[1] else NA_character_
}
# one cell list-column -> character timestamps (Excel datetimes, serial numbers, or text)
cell_str <- function(col) vapply(col, function(v) {
  if (length(v) == 0 || (length(v) == 1 && is.na(v))) return(NA_character_)
  if (inherits(v, "POSIXct")) return(format(v, "%Y-%m-%d %H:%M:%S", tz = "UTC"))
  if (is.numeric(v)) return(format(as.POSIXct(round((v - 25569) * 86400), origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d %H:%M:%S"))
  as.character(v)
}, character(1), USE.NAMES = FALSE)
cell_num <- function(col) suppressWarnings(vapply(col, function(v) {
  if (length(v) == 0) return(NA_real_); if (is.numeric(v)) return(as.numeric(v)); as.numeric(as.character(v)) }, numeric(1), USE.NAMES = FALSE))
parse_dt <- function(x, tz = "UTC") {
  out <- as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = tz)
  for (fm in c("%Y-%m-%d %H:%M", "%Y-%m-%d")) { i <- is.na(out) & !is.na(x); if (any(i)) out[i] <- as.POSIXct(x[i], format = fm, tz = tz) }
  out
}
# timestamps exactly as the published sheet lays them out (local wall clock, returned as UTC-tagged POSIXct)
sheet_times <- function(raw) {
  nm <- names(raw)
  if ("DateTime" %in% nm) return(list(t = parse_dt(cell_str(raw[["DateTime"]])), date_full = NULL))
  d <- cell_str(raw[[1]])
  if (length(nm) < 2 || !(nm[2] %in% c("Time", "None", "...2"))) return(list(t = parse_dt(d), date_full = parse_dt(d)))
  tt <- regmatches(cell_str(raw[[2]]), regexpr("\\d{1,2}:\\d{2}(:\\d{2})?$", cell_str(raw[[2]])))
  tv <- rep(NA_character_, length(d)); ok <- grepl("\\d{1,2}:\\d{2}(:\\d{2})?$", cell_str(raw[[2]])); tv[ok] <- tt
  list(t = parse_dt(ifelse(is.na(tv) | is.na(d), NA, paste(substr(d, 1, 10), tv))), date_full = parse_dt(d))
}
# dump the needed sheets of one workbook to gzipped CSV (cached; a .done file records completion)
extract_workbook <- function(f, out = file.path(INTERIM, "ccnd_sheets")) {
  dir.create(out, showWarnings = FALSE, recursive = TRUE); tag <- tag_of(f); done <- file.path(out, paste0(tag, ".done"))
  if (file.exists(done)) return(invisible(tag))
  sn <- readxl::excel_sheets(f); got <- character()
  for (key in names(SHEETS)) {
    s <- find_sheet(sn, key); if (is.na(s)) next
    got <- c(got, paste0(key, ":", s)); dest <- file.path(out, paste0(tag, "_", key, ".csv.gz")); if (file.exists(dest)) next
    raw <- suppressMessages(suppressWarnings(readxl::read_excel(f, s, col_types = "list", .name_repair = "unique")))
    tm <- sheet_times(raw)
    site_cols <- grep("CM(\\d+)", names(raw), value = TRUE)
    num <- as.integer(sub(".*CM(\\d+).*", "\\1", site_cols)); site_cols <- site_cols[num >= 1 & num <= 10]
    dt <- data.table(t_local = format(tm$t, "%Y-%m-%d %H:%M:%S"),
                     date_full = if (is.null(tm$date_full)) NA_character_ else format(tm$date_full, "%Y-%m-%d %H:%M:%S"))
    for (c in site_cols) { nmc <- paste0("CM", sub(".*CM(\\d+).*", "\\1", c)); if (!nmc %in% names(dt)) set(dt, j = nmc, value = cell_num(raw[[c]])) }
    fwrite(dt, dest)
  }
  writeLines(got, done); invisible(tag)
}
load_sheet <- function(tag, key, src = file.path(INTERIM, "ccnd_sheets")) {
  f <- file.path(src, paste0(tag, "_", key, ".csv.gz")); if (!file.exists(f)) return(NULL)
  d <- fread(f, colClasses = list(character = c("t_local", "date_full")))
  d[, t_local := parse_dt(t_local)][, date_full := parse_dt(date_full)]; d[]
}
# naive America/Denver wall clock (UTC-tagged POSIXct) -> fixed MST POSIXct; ambiguous / nonexistent hours -> NA
local_to_mst <- function(t) {
  s <- format(t, "%Y-%m-%d %H:%M:%S", tz = "UTC")
  c6 <- t + 6*3600; c7 <- t + 7*3600                                  # candidate UTC instants (MDT, MST)
  ok6 <- format(c6, "%Y-%m-%d %H:%M:%S", tz = "America/Denver") == s
  ok7 <- format(c7, "%Y-%m-%d %H:%M:%S", tz = "America/Denver") == s
  u <- rep(as.POSIXct(NA, tz = "UTC"), length(t)); one <- xor(ok6, ok7) & !is.na(ok6) & !is.na(ok7)
  u[one & ok6] <- c6[one & ok6]; u[one & ok7] <- c7[one & ok7]
  attr(u, "tzone") <- TZ_MST; u
}
fix12 <- function(t) {
  d <- c(0, diff(as.numeric(t))); add <- numeric(length(t)); cur <- 0
  for (i in seq_along(t)[-1]) { if (!is.na(d[i]) && d[i] < -40000) cur <- 12 else if (!is.na(d[i]) && d[i] > 80000) cur <- 0; add[i] <- cur }
  t + add*3600
}
# harmonize one quarter (cached): list(gases = list(key -> data.table t_mst + sites), wind = data.table site,t5,u,v,ws, range)
quarter_frames <- function(tag) {
  cache <- file.path(INTERIM, "ccnd_quarters", paste0(tag, ".rds")); if (file.exists(cache)) return(readRDS(cache))
  qy <- quarter_of(tag); s <- as.POSIXct(sprintf("%d-%02d-01", qy["y"], 3*qy["q"] - 2), tz = "UTC")
  e <- as.POSIXct(seq(as.Date(s), by = "3 months", length.out = 2)[2], tz = "UTC")
  G <- list()
  for (key in setdiff(names(SHEETS), c("ws", "wd"))) {
    d <- load_sheet(tag, key); if (is.null(d)) next
    d <- d[!is.na(t_local) & t_local >= s & t_local < e]; d <- d[!duplicated(t_local)]
    d[, t_mst := local_to_mst(t_local)]; d <- d[!is.na(t_mst)][!duplicated(t_mst)]
    G[[key]] <- d[, c("t_mst", intersect(SITE_LIST, names(d))), with = FALSE]
  }
  W <- NULL; ws <- load_sheet(tag, "ws"); wd <- load_sheet(tag, "wd")
  if (!is.null(ws) && !is.null(wd)) {
    if (startsWith(tag, "CCND_Q4_2023")) { ws[, t_local := date_full]; wd[, t_local := date_full] }
    if (startsWith(tag, "CCND_Q1_2024")) for (x in list(ws, wd)) x[, t_local := as.POSIXct(format(t_local, tz = "America/Denver", usetz = FALSE), tz = "UTC")]
    if (startsWith(tag, "CCND_Q1_2025")) { ws[, t_local := fix12(t_local)]; wd[, t_local := fix12(t_local)] }
    # timestamps repeated >2 times are unusable (Q3 2021: from 2021-09-04 22:02 stamps collapse to midnight);
    # pairs (Q1 2025 after the AM/PM fix) are kept and averaged, as in the original analysis
    dropdup <- function(x) { x <- x[!is.na(t_local)]; x[, n_ := .N, by = t_local]; x <- x[n_ <= 2]; x[, n_ := NULL]; x }
    ws <- dropdup(ws); wd <- dropdup(wd)
    sc <- intersect(SITE_LIST, intersect(names(ws), names(wd)))
    a <- melt(ws[!is.na(t_local), c("t_local", sc), with = FALSE], id.vars = "t_local", variable.name = "site", value.name = "ws")
    b <- melt(wd[!is.na(t_local), c("t_local", sc), with = FALSE], id.vars = "t_local", variable.name = "site", value.name = "wd")
    m <- merge(a, b, by = c("t_local", "site"), allow.cartesian = TRUE)
    m <- m[ws >= 0 & wd >= 0 & wd <= 360 & t_local >= s & t_local < e]
    m[, t := local_to_mst(t_local)]; m <- m[!is.na(t)]
    m[, `:=`(u = -ws*0.44704*sin(wd*pi/180), v = -ws*0.44704*cos(wd*pi/180), wsm = ws*0.44704)]
    m[, t5 := as.POSIXct(ceiling(as.numeric(t)/300)*300, origin = "1970-01-01", tz = TZ_MST)]
    W <- m[, .(u = mean(u), v = mean(v), ws = mean(wsm)), by = .(site, t5)]
  }
  r <- list(gases = G, wind = W, range = c(local_to_mst(s), local_to_mst(e)))
  dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE); saveRDS(r, cache); r
}
GRID <- function(t0 = "2021-08-16 00:00:00", t1 = "2026-01-01 00:00:00")
  seq(as.POSIXct(t0, tz = TZ_MST), as.POSIXct(t1, tz = TZ_MST), by = 300)
build_panel <- function(tags) {
  grid <- GRID(); n <- length(grid); quarter <- rep("", n); acc <- list(); W <- list()
  for (tag in tags) {
    r <- quarter_frames(tag)
    for (k in names(r$gases)) acc[[k]] <- c(acc[[k]], list(r$gases[[k]]))
    if (!is.null(r$wind)) W[[tag]] <- r$wind
    qy <- quarter_of(tag); quarter[grid >= r$range[1] & grid < r$range[2]] <- sprintf("Q%d_%d", qy["q"], qy["y"])
  }
  P <- list(t = grid, quarter = quarter)
  to_mat <- function(d) {
    d <- d[!duplicated(t_mst)]; i <- match(as.numeric(d$t_mst), as.numeric(grid))
    M <- matrix(NA_real_, n, 10, dimnames = list(NULL, SITE_LIST)); ok <- !is.na(i)
    for (s in intersect(SITE_LIST, names(d))) M[i[ok], s] <- d[[s]][ok]
    M
  }
  for (k in names(acc)) P[[k]] <- to_mat(rbindlist(acc[[k]], fill = TRUE))
  W <- rbindlist(W); W <- W[!duplicated(W[, .(t5, site)])]
  i <- match(as.numeric(W$t5), as.numeric(grid)); j <- match(as.character(W$site), SITE_LIST); ok <- !is.na(i) & !is.na(j)
  for (k in c("u", "v", "ws")) { M <- matrix(NA_real_, n, 10, dimnames = list(NULL, SITE_LIST)); M[cbind(i[ok], j[ok])] <- W[[k]][ok]; P[[k]] <- M }
  P
}
fmt_t <- function(t) format(t, "%Y-%m-%d %H:%M:%S", tz = TZ_MST)
save_panel <- function(P) {
  wide <- function(keys) { d <- data.table(t_mst = fmt_t(P$t))
    for (k in keys) if (!is.null(P[[k]])) { M <- signif(P[[k]], 6); colnames(M) <- paste0(k, "_", SITE_LIST); d <- cbind(d, as.data.table(M)) }
    d }
  h <- wide("h2s"); h[, quarter := P$quarter]; setcolorder(h, c("t_mst", "quarter"))
  fwrite(h, file.path(PROCESSED, "ccnd_5min_h2s.csv.gz"))
  fwrite(wide(c("u", "v", "ws")), file.path(PROCESSED, "ccnd_5min_wind.csv.gz"))
  fwrite(wide(GASES), file.path(PROCESSED, "ccnd_5min_gases.csv.gz"))
}
load_panel <- function(gases = FALSE) {
  rd <- function(f) fread(file.path(PROCESSED, f), colClasses = list(character = "t_mst"))
  h <- rd("ccnd_5min_h2s.csv.gz"); P <- list(t = as.POSIXct(h$t_mst, tz = TZ_MST), quarter = h$quarter)
  getm <- function(d, k) { M <- as.matrix(d[, paste0(k, "_", SITE_LIST), with = FALSE]); colnames(M) <- SITE_LIST; storage.mode(M) <- "double"; M }
  P$h2s <- getm(h, "h2s"); w <- rd("ccnd_5min_wind.csv.gz"); for (k in c("u", "v", "ws")) P[[k]] <- getm(w, k)
  if (gases) { g <- rd("ccnd_5min_gases.csv.gz"); for (k in GASES) if (paste0(k, "_CM1") %in% names(g)) P[[k]] <- getm(g, k) }
  P
}
