# Step 1: dump the needed sheets of every CCND workbook to data/interim/ccnd_sheets (cached, resumable).
.repo <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "..")); source(file.path(.repo, "R", "config.R")); src_files()
files <- sort(list.files(RAW_CCND, "\\.xlsx$", full.names = TRUE))
if (!length(files)) stop("No CCND workbooks in ", RAW_CCND, " (set EMIS_RAW_CCND)")
for (f in files) { t0 <- Sys.time(); extract_workbook(f); dt <- as.numeric(Sys.time() - t0, units = "secs"); if (dt > 1) cat(sprintf("%s: %.0fs\n", tag_of(f), dt)) }
cat("extraction step done\n")
