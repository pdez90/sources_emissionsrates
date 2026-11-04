# Step 99: vet the run. Checks data integrity, physics/unit consistency, the reference validations, internal consistency
# of the estimates, and that the headline numbers match tests/expected.json within tolerance. Exits 1 if any check fails.
.repo <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "..")); source(file.path(.repo, "R", "config.R")); src_files()
K <- load_keys(); EXP <- jsonlite::read_json(file.path(REPO, "tests", "expected.json")); fails <- character()
check <- function(name, ok, detail = "") { ok <- isTRUE(ok); cat(if (ok) "PASS " else "FAIL ", name, if (nzchar(detail)) paste0("  [", detail, "]"), "\n", sep = "")
  if (!ok) fails <<- c(fails, name) }
# 1 data integrity
P <- load_panel(); cov <- colMeans(is.finite(P$h2s))
check("CCND panel shape 460513 x 10", identical(dim(P$h2s), c(460513L, 10L)), paste(dim(P$h2s), collapse = " x "))
check("CCND H2S coverage >= 75% at 9+ sites", sum(cov >= 0.75) >= 9, paste(round(cov, 2), collapse = " "))
check("CCND H2S non-negative (published floor)", min(P$h2s, na.rm = TRUE) >= 0, format(min(P$h2s, na.rm = TRUE)))
b <- birch_load_1min()
check("Birch 1-min spans 2025Q4-2026Q2", min(b$t) <= as.POSIXct("2025-10-01 00:05", tz = TZ_MST) && max(b$t) >= as.POSIXct("2026-06-30", tz = TZ_MST), paste(format(range(b$t)), collapse = " .. "))
check("Birch keeps sub-zero values (MD qualifier)", any(b$h2s < 0, na.rm = TRUE))
check("Birch wind speed plausible (median 0.5-3 m/s)", median(b$ws, na.rm = TRUE) > 0.5 && median(b$ws, na.rm = TRUE) < 3, sprintf("%.2f", median(b$ws, na.rm = TRUE)))
# 2 physics / units
for (s in c("CM8", "BIRCH", "CM2")) {
  bd <- bearing(SOURCES$HitePoint[1:2], ALL_SITES[[s]]); th <- bd["bearing"]*pi/180
  g <- footprint(s, 2.5*sin(th), 2.5*cos(th), 12, list(HitePoint = 1)); f <- centerline_factor(bd["dist"], 2.5, 13, H = 3)
  check(paste("centerline consistency", s), abs(g*f - 1) < 0.02, sprintf("footprint*factor=%.3f", g*f))
}
check("upwind footprint is zero", footprint("CM8", 0, -2, 12, list(HitePoint = 1)) == 0)
# 3 reference validations
v <- K$validation
if (!is.null(v$clock_lag_raw_min)) {
  check("Birch clock fix: raw H2S lags own wind by ~6 h", v$clock_lag_raw_min >= 300 && v$clock_lag_raw_min <= 390, v$clock_lag_raw_min)
  check("Birch clock fix: fixed lag within 1 h", abs(v$clock_lag_fixed_min) <= 60, v$clock_lag_fixed_min)
} else cat("SKIP Birch clock diagnostic (needs the raw MT GOAT zips)\n")
check("CM4 vs Birch timing: >= 75% of 3-day windows within 30 min", v$lag_frac_within_30min >= 0.75, round(v$lag_frac_within_30min, 2))
check("CM4 gain (Deming) 0.6-0.85", v$cm4_deming_slope >= 0.6 && v$cm4_deming_slope <= 0.85, round(v$cm4_deming_slope, 3))
check("CM8 direction within 10 deg of Birch", abs(v$cm8_dir_offset_vs_birch) <= 10, round(v$cm8_dir_offset_vs_birch, 1))
check("CM8 direction stable vs DEN (-20..0 deg every quarter)", v$cm8_dir_vs_den_range[1] >= -20 && v$cm8_dir_vs_den_range[2] <= 0, paste(round(v$cm8_dir_vs_den_range, 1), collapse = ".."))
# 4 internal consistency
q4 <- K$invert_q4_2025; yr <- K$invert_2024_2025
check("plume-conditioned Hite > zeros-retained (selection bias)", q4$plumecond_Hite_median > q4$hite_split_transport)
check("area split fits better than north-only (site spread)", q4$spread_split < q4$spread_north, sprintf("%.2f vs %.2f", q4$spread_split, q4$spread_north))
check("CM8-scaled wind reproduces Birch-wind Q4 Hite within 30%", abs(yr$hite_q4_cm8wind_transport/q4$hite_split_transport - 1) < 0.3, sprintf("%.1f vs %.1f", yr$hite_q4_cm8wind_transport, q4$hite_split_transport))
# 5 regression against expected values
for (sec in names(EXP)) for (k in names(EXP[[sec]])) {
  e <- EXP[[sec]][[k]]; got <- K[[sec]][[k]]
  check(sprintf("expected %s.%s = %s +/- %.0f%%", sec, k, format(e[[1]]), 100*e[[2]]), !is.null(got) && abs(got - e[[1]]) <= e[[2]]*abs(e[[1]]) + 1e-9, paste("got", format(got)))
}
cat(sprintf("\n%d failure(s)\n", length(fails))); quit(status = if (length(fails)) 1 else 0)
