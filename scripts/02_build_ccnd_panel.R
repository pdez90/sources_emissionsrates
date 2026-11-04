# Step 2: harmonize all extracted quarters -> data/processed/ccnd_5min_{h2s,wind,gases}.csv.gz (MST, 5-min).
.repo <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "..")); source(file.path(.repo, "R", "config.R")); src_files()
tags <- sort(sub("\\.done$", "", list.files(file.path(INTERIM, "ccnd_sheets"), "\\.done$")))
for (tag in tags) if (!file.exists(file.path(INTERIM, "ccnd_quarters", paste0(tag, ".rds")))) {
  t0 <- Sys.time(); quarter_frames(tag); cat(sprintf("%s: %.0fs\n", tag, as.numeric(Sys.time() - t0, units = "secs"))) }
P <- build_panel(tags); save_panel(P)
cat("panel", dim(P$h2s), " H2S coverage by site:", round(colMeans(is.finite(P$h2s)), 2), "\n")
