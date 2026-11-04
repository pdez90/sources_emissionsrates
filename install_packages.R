# One-time setup: install the R packages the pipeline uses.
pkgs <- c("data.table", "R.utils", "readxl", "jsonlite", "caTools", "nnls")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
cat("R", R.version$major, ".", R.version$minor, " packages OK:", paste(pkgs, collapse = ", "), "\n")
