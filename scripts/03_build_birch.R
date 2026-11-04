# Step 3: parse the MT GOAT 1-s zips -> data/processed/birch_1min.csv.gz (cached per quarter/parameter).
.repo <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "..")); source(file.path(.repo, "R", "config.R")); src_files()
qs <- birch_quarters(); if (!length(qs)) stop("No MT GOAT zips in ", RAW_BIRCH, " (set EMIS_RAW_BIRCH)")
for (q in qs) for (p in BIRCH_PARAMS) if (!file.exists(file.path(INTERIM, "birch", paste0(p, "_", q, ".rds")))) {
  t0 <- Sys.time(); birch_parse(q, p); cat(sprintf("%s %s: %.0fs\n", q, p, as.numeric(Sys.time() - t0, units = "secs"))) }
d <- birch_build(); cat("birch 1-min rows", nrow(d), " h2s valid", round(mean(!is.na(d$h2s)), 2), " wind valid", round(mean(!is.na(d$u)), 2), "\n")
