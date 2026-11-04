# Step 4: parse NOAA ISD DEN hourly wind -> data/processed/den_hourly.csv.gz
.repo <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "..")); source(file.path(.repo, "R", "config.R")); src_files()
d <- den_build(); cat("DEN hourly", nrow(d), format(min(d$t_mst)), "..", format(max(d$t_mst)), "\n")
