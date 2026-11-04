# Step 10: reference checks. Writes results/validation_*.csv, results/wind_qa_vs_den.csv and key_numbers.json[validation].
.repo <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "..")); source(file.path(.repo, "R", "config.R")); src_files()
P <- load_panel(); b5 <- birch_five_min(); b1h <- birch_1h_ccnd_style(b5); cm4 <- P$h2s[, "CM4"]
lt <- lag_table(b5$t, b1h, P$t, cm4); wcsv(lt, file.path(RESULTS, "validation_lag_birch_vs_cm4.csv"))
cat("CM4-vs-Birch best lag by 3-day window (min):", lt$lag_min, "\n")
G <- gain_fit(b5$t, b1h, P$t, cm4); print(round(unlist(G), 3))
do.call(record, c(list("validation"), setNames(G, paste0("cm4_", names(G)))))
record("validation", lag_median_min = median(lt$lag_min, na.rm = TRUE), lag_frac_within_30min = mean(abs(lt$lag_min) <= 30, na.rm = TRUE))
W <- wind_compare(P, b5); wcsv(W, file.path(RESULTS, "validation_wind_birch_vs_ccnd.csv")); print(W)
record("validation", cm8_dir_offset_vs_birch = W[site == "CM8", dir_offset_med], cm8_speed_ratio_vs_birch = W[site == "CM8", speed_ratio])
Q <- wind_qa(P, den_load()); wcsv(Q, file.path(RESULTS, "wind_qa_vs_den.csv"))
q <- Q[grepl("2024|2025", quarter)]; print(q[, c("quarter", grep("_dir$", names(q), value = TRUE)), with = FALSE], digits = 2)
record("validation", cm8_dir_vs_den_range = range(q$CM8_dir, na.rm = TRUE))
if (file.exists(file.path(INTERIM, "birch", "hydrogen_sulfide_2025Q4.rds"))) {      # clock diagnostic needs the raw-parse cache
  raw <- birch_series("2025Q4", "hydrogen_sulfide", fix_clock = FALSE)
  raw[, b := ceiling(t/300)*300]; r5 <- raw[, .(h = mean(v, na.rm = TRUE)), by = b]
  r1h <- roll_mean(pmax(r5$h[match(as.numeric(b5$t), r5$b)], 0), 12, 12)
  cc <- clock_check(b5$t, r1h, b1h, b5); cat("Birch H2S vs own S-sector wind lag (min, r): raw", cc$raw, " fixed", cc$fixed, "\n")
  record("validation", clock_lag_raw_min = cc$raw[1], clock_lag_fixed_min = cc$fixed[1])
}
