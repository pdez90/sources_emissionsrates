# Step 30: emission estimates for Hite and the refinery complex (Suncor + P66).
# Part A: Q4 2025 with Birch wind and event masks (CCND + Birch).  Part B: 2024-2025 with CM8-scaled regional wind (CCND).
.repo <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "..")); source(file.path(.repo, "R", "config.R")); src_files()
P <- load_panel(gases = TRUE); b1 <- birch_load_1min(); b5 <- birch_five_min(b1)
EV <- fread(file.path(RESULTS, "q4_2025_events.csv"), colClasses = list(character = c("t0", "t1"))); SR <- fread(file.path(RESULTS, "q4_2025_event_sites.csv"), colClasses = list(character = "tpeak"))
for (k in c("t0", "t1")) EV[[k]] <- as.POSIXct(EV[[k]], tz = TZ_MST); SR[, tpeak := as.POSIXct(tpeak, tz = TZ_MST)]
rows <- list(); add <- function(...) { r <- list(...); rows[[length(rows) + 1]] <<- as.data.table(r)
  cat(paste(names(r), sapply(r, function(v) if (is.numeric(v)) sprintf("%.1f", v) else v), sep = "=", collapse = " "), "\n") }
# ---------------- Part A: Q4 2025, Birch wind ----------------
pq <- panel_q4(P, b5); tt <- pq$t; W <- regional_wind(P, tt, birch5 = b5); U <- W$U; V <- W$V; sp <- sqrt(U^2 + V^2)
hr <- mst_hour(tt); day <- day_of(tt); enh <- enhancement(pq$X$h2s); sites <- setdiff(colnames(enh), "CM6")
event_mask <- function(classes) { m <- rep(TRUE, length(tt)); for (r in which(EV$final %in% classes)) m[tt >= EV$t0[r] & tt <= EV$t1[r]] <- FALSE; m }
mH <- event_mask(c("Refinery", "other/local")); mR <- event_mask(c("Hite", "other/local"))
trans <- !is.na(sp) & sp >= TRANSPORT_MIN_WS; allh <- !is.na(sp)
geoms <- list(split = HITE_SPLIT, main = list(HiteMain = 1), north = list(HiteNorth = 1), point = list(HitePoint = 1))
for (g in names(geoms)) {
  FP <- footprints(sites, U, V, hr, parts = list(Hite = geoms[[g]], Suncor = list(Suncor = 1), P66 = list(P66 = 1)))
  for (regime in c("transport", "all")) {
    base <- (if (regime == "transport") trans else allh) & mH
    r <- ratio_by_site(enh, FP, "Hite", base, day); q9 <- r$q[setdiff(names(r$q), "CM9")]; ci <- boot_site_median(r$daily[setdiff(names(r$daily), "CM9")])
    add(part = "A_Q4_2025_birchwind", target = "Hite", geometry = g, regime = regime, site_median = median(q9), ci_lo = ci[1], ci_hi = ci[3],
        spread_log = site_spread(r$q), birch = unname(r$q["BIRCH"]), cm9 = unname(r$q["CM9"]))
    if (g == "split") fwrite(data.table(site = names(r$q), Q_tpy = r$q), file.path(RESULTS, sprintf("q4_2025_hite_by_site_%s.csv", regime)))
  }
  if (g == "split") {
    for (gn in c(0.64, 0.83, 1.0)) { r <- ratio_by_site(enh, FP, "Hite", trans & mH, day, gain = gn)
      add(part = "A_Q4_2025_birchwind", target = "Hite", geometry = "split", regime = paste0("transport_gain", gn), site_median = median(r$q[setdiff(names(r$q), "CM9")])) }
    for (dm in c(1, 1e9)) { r <- ratio_by_site(enh, FP, "Hite", trans & mH, day, dom = dm)     # weaker / no dominance filter
      add(part = "A_Q4_2025_birchwind", target = "Hite", geometry = "split", regime = paste0("transport_dom", format(dm)), site_median = median(r$q[setdiff(names(r$q), "CM9")])) }
    for (tg in c("Refinery", "Suncor")) { r <- ratio_by_site(enh, FP, tg, trans & mR, day); ci <- boot_site_median(r$daily[setdiff(names(r$daily), "CM9")])
      add(part = "A_Q4_2025_birchwind", target = tg, geometry = "point", regime = "transport", site_median = median(r$q[setdiff(names(r$q), "CM9")]), ci_lo = ci[1], ci_hi = ci[3]) }
  }
}
# coverage: Hite-downwind transport time that falls inside detected Hite events
wf <- wind_from(U, V); dw <- rep(FALSE, length(tt))
for (s in sites) for (src in c("HiteMain", "HiteNorth")) { bd <- bearing(ALL_SITES[[s]], SOURCES[[src]][1:2]); if (bd["dist"] <= 3500) dw <- dw | (!is.na(wf) & abs(angdiff(wf, bd["bearing"])) <= 20) }
dw <- dw & trans; inev <- !event_mask("Hite"); cov <- sum(dw & inev)/sum(dw)
cat(sprintf("Hite-downwind transport time %.1f%%; inside detected Hite events %.0f%%\n", 100*mean(dw), 100*cov))
record("invert_q4_2025", hite_downwind_frac = mean(dw), hite_event_coverage = cov)
for (tg in c("Hite", "Refinery")) { pc <- plume_conditioned(EV, SR, b1, target = tg); wcsv(pc$Q, file.path(RESULTS, sprintf("q4_2025_plume_conditioned_%s.csv", tg)))
  cat("plume-conditioned", tg, ":", paste(names(pc$s), round(unlist(pc$s), 1), sep = "=", collapse = " "), "\n")
  do.call(record, c(list("invert_q4_2025"), setNames(pc$s, paste0("plumecond_", tg, "_", names(pc$s))))) }
# CM9 directional diagnostic
F9 <- footprint("CM9", U, V, hr, HITE_SPLIT); m <- trans & is.finite(enh[, "CM9"])
fwrite(data.table(wd_bin = (wf[m] %/% 15)*15, obs = enh[m, "CM9"]/CCND_GAIN, pred35 = F9[m]*35/TPY_PER_GS)[, .(obs_ppb = mean(obs), pred_ppb_at35tpy = mean(pred35), n = .N), keyby = wd_bin],
       file.path(RESULTS, "q4_2025_cm9_by_direction.csv"))
# ---------------- Part B: 2024-2025, CM8-scaled regional wind ----------------
t2 <- seq(as.POSIXct("2024-01-01 00:05:00", tz = TZ_MST), as.POSIXct("2026-01-01 00:00:00", tz = TZ_MST), by = 300)
W2 <- regional_wind(P, t2); U2 <- W2$U; V2 <- W2$V; sp2 <- sqrt(U2^2 + V2^2); hr2 <- mst_hour(t2); day2 <- day_of(t2)
enh2 <- enhancement(P$h2s[match(as.numeric(t2), as.numeric(P$t)), ]); sites2 <- setdiff(SITE_LIST, c("CM6", "CM9"))
FP2 <- footprints(sites2, U2, V2, hr2)
cm8q1 <- list(CM8 = !(t2 >= as.POSIXct("2025-01-01", tz = TZ_MST) & t2 < as.POSIXct("2025-04-01", tz = TZ_MST)))   # CM8 Q1 2025 excluded
periods <- list("2024" = t2 < as.POSIXct("2025-01-01", tz = TZ_MST), "2025" = t2 >= as.POSIXct("2025-01-01", tz = TZ_MST), "Q4_2025" = t2 >= as.POSIXct("2025-10-01", tz = TZ_MST))
for (pn in names(periods)) for (regime in c("transport", "all")) {
  base <- periods[[pn]] & !is.na(sp2) & (if (regime == "transport") sp2 >= TRANSPORT_MIN_WS else TRUE)
  for (tg in c("Hite", "Refinery", "Suncor")) {
    r <- ratio_by_site(enh2, FP2, tg, base, day2, site_mask = cm8q1); ci <- boot_site_median(r$daily, nb = 300)
    add(part = "B_cm8wind", target = tg, geometry = if (tg == "Hite") "split" else "point", regime = regime, period = pn, site_median = median(r$q), ci_lo = ci[1], ci_hi = ci[3], site_min = min(r$q), site_max = max(r$q))
    if (pn == "2025" && regime == "transport") fwrite(data.table(site = names(r$q), Q_tpy = r$q), file.path(RESULTS, sprintf("2025_%s_by_site_transport.csv", tg)))
  }
  nn <- binned_nnls(enh2, FP2, base, wind_from(U2, V2), sp2, hr2, site_mask = cm8q1)
  add(part = "B_cm8wind", target = "NNLS", regime = regime, period = pn, hite = nn[["Hite"]], suncor = nn[["Suncor"]], p66 = nn[["P66"]])
}
R <- rbindlist(rows, fill = TRUE); fwrite(R, file.path(RESULTS, "estimates.csv"))
pick <- function(...) { k <- list(...); r <- R; for (n in names(k)) r <- r[r[[n]] == k[[n]]]; r$site_median[1] }
record("invert_q4_2025", hite_split_transport = pick(part = "A_Q4_2025_birchwind", target = "Hite", geometry = "split", regime = "transport"),
       hite_split_all = pick(part = "A_Q4_2025_birchwind", target = "Hite", geometry = "split", regime = "all"),
       spread_split = R[part == "A_Q4_2025_birchwind" & geometry == "split" & regime == "transport", spread_log][1],
       spread_north = R[part == "A_Q4_2025_birchwind" & geometry == "north" & regime == "transport", spread_log][1],
       refinery_transport = pick(part = "A_Q4_2025_birchwind", target = "Refinery", regime = "transport"))
record("invert_2024_2025", hite_2025_transport = pick(part = "B_cm8wind", target = "Hite", regime = "transport", period = "2025"),
       hite_2025_all = pick(part = "B_cm8wind", target = "Hite", regime = "all", period = "2025"),
       refinery_2025_transport = pick(part = "B_cm8wind", target = "Refinery", regime = "transport", period = "2025"),
       hite_2024_transport = pick(part = "B_cm8wind", target = "Hite", regime = "transport", period = "2024"),
       hite_q4_cm8wind_transport = pick(part = "B_cm8wind", target = "Hite", regime = "transport", period = "Q4_2025"))
cat("wrote", file.path(RESULTS, "estimates.csv"), "\n")
